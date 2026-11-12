"""Six-state TVSA model coupled to a lumped bed energy balance."""

from math import exp

import numpy as np

from gas_balance import (
    effective_aw,
    flows_from_pressure,
    pressure_chain,
    select_p_co2,
    update_y,
)
from isotherms import co2_capacity, co2_equilibrium, h2o_equilibrium
from kinetics import co2_slow_residual, h2o_bound_residual, pool_fraction, rates
from parameters import equipment
from thermal import (
    adsorption_heat_rate,
    cooler_controller,
    cooler_power,
    dynamic_energy_balance,
    eff_heat_capacity,
    gas_flow_heat_rate,
    heat_loss_rate,
    heater_controller,
    heater_power,
)


def clamp(x, lo, hi):
    return min(max(x, lo), hi)


def flow_ads(ctrl, design, opt):
    """Return the inlet volumetric flow [m3/s] using the selected area convention."""
    return ctrl["u_feed_m_s"] * (
        design["A_contact_m2"]
        if opt["flow_mode"].lower() == "acontact"
        and opt["legacy"]["keepQ_with_Acontact"]
        else design["A_flow_m2"]
    )


def forward_cycle(drv, amb, ctrl, design, phys, opt, x_init=None):
    """Advance one cycle using the bed energy balance for temperature.

    State columns: CO2 fast/slow and H2O free/bound loadings [mol/kg],
    gas CO2 mole fraction [-], and bed temperature [K]. Heat diagnostics
    are in W; x_end can initialise the next cycle."""
    dt = opt["dt_s"]
    t = drv["t_s"]
    sid_all = drv["stage_id"]
    n = len(t)
    R = phys["eq"]["R"]
    # State ordering [qCO2_fast, qCO2_slow, qH2O_free, qH2O_bound, yCO2_g, T_bed]
    x = np.zeros((n, 6), dtype=float)
    if x_init is None:
        tk = drv["TK_cmd"][0]
        pk = drv["P_Pa"][0]
        ph, pdry, aw = pressure_chain(drv["RH_cmd"][0], tk, tk, pk, opt)
        qh = h2o_equilibrium(tk, aw, phys["iso"])
        fb = pool_fraction(aw, int(sid_all[0]), opt["h2o_2pool"]["fbound"])
        qco2 = co2_equilibrium(tk, amb["yCO2_air"] * pdry, qh, phys["iso"])
        aw_eff = effective_aw(
            int(sid_all[0]), tk, aw, (1 - fb) * qh, fb * qh, phys["iso"], opt
        )
        fs = pool_fraction(aw_eff, int(sid_all[0]), opt["co2_2pool"]["fslow"])
        x[0] = [(1 - fs) * qco2, fs * qco2, (1 - fb) * qh, fb * qh, amb["yCO2_air"], tk]
    else:
        x[0] = np.asarray(x_init, dtype=float).reshape(6)
        x[0, :4] = np.maximum(x[0, :4], 0.0)
        x[0, 4] = clamp(x[0, 4], opt["yco2"]["y_floor"], opt["yco2"]["y_max"])
        x[0, 5] = max(x[0, 5], 200.0)
    # Diagnostic arrays in MATLAB state terms.
    variables = [
        "qCO2eq",
        "pCO2_inv",
        "pCO2_drive",
        "p_h2o",
        "pdry",
        "aw_gas",
        "aw_eff",
        "q_in",
        "q_out",
        "y_in",
        "n_dry",
        "nd_res",
        "kco2_fast",
    ]
    variables += [
        "Q_heater_W",
        "Q_ads_heat_W",
        "Q_cooler_W",
        "Q_flow_W",
        "Q_loss_W",
        "Q_ads_CO2_W",
        "Q_ads_H2O_W",
    ]
    z = {name: np.zeros(n) for name in variables}
    qads_prev = x[0, 0] + x[0, 1]
    peak_h2o_bound = x[0, 3]
    Qin_filt = 0.0
    Qout_filt = 0.0
    n_dry_filt_prev = float("nan")
    pdrv_filt_prev = float("nan")
    aw_filt_prev = float("nan")
    Q_ads = flow_ads(ctrl, design, opt)
    u_h_prev = 0  # Heater initially OFF
    u_c_prev = 0
    for k in range(n - 1):
        sid = int(sid_all[k])
        sid_prev = int(sid_all[max(k - 1, 0)])
        leaving_ads = sid != 6 and sid_prev == 6
        if leaving_ads:
            qads_prev = x[k, 0] + x[k, 1]
        stage_changed = int(sid_all[k + 1]) != sid
        TK_cmd = float(drv["TK_cmd"][k])
        rh = float(drv["RH_cmd"][k])
        pk = float(drv["P_Pa"][k])
        tk = x[k, 5]
        # Stage names index the pressure and dry-gas filter settings.
        names = ("vac", "heat", "des", "cool", "press", "ads")

        ph, pdry, aw = pressure_chain(rh, TK_cmd, tk, pk, opt)

        # CO2-effective water activity based on water state and GAB inverse.
        aw_eff = effective_aw(sid, tk, aw, x[k, 2], x[k, 3], phys["iso"], opt)
        if opt["stab"]["enable"] and opt["stab"]["aw_lpf_enable"]:
            if not np.isfinite(aw_filt_prev):
                aw_filt_prev = aw_eff
            if stage_changed and opt["stab"]["aw_stagejump_reset"]:
                aw_filt_prev = aw_eff
            else:
                aa = 1 - exp(-dt / max(opt["stab"]["aw_tau_s"], 1e-6))
                aw_filt_prev += aa * (aw_eff - aw_filt_prev)
            aw_eff = clamp(aw_filt_prev, 0.0, 0.999)
        p_inv = clamp(x[k, 4], opt["yco2"]["y_floor"], opt["yco2"]["y_max"]) * pdry
        p_in = amb["yCO2_air"] * pdry
        pdrv = select_p_co2(sid, p_inv, p_in, opt)
        if opt["stab"]["enable"] and opt["stab"]["pdrv_lpf_enable"]:
            if not np.isfinite(pdrv_filt_prev):
                pdrv_filt_prev = pdrv
            if stage_changed and opt["stab"]["pdrv_stagejump_reset"]:
                pdrv_filt_prev = pdrv
            else:
                aa = 1 - exp(-dt / max(opt["stab"]["pdrv_tau_s"], 1e-6))
                pdrv_filt_prev += aa * (pdrv - pdrv_filt_prev)
            pdrv = max(pdrv_filt_prev, 0.0)
        qh_eq = h2o_equilibrium(tk, aw, phys["iso"], opt["aw"]["max"])
        fb = pool_fraction(aw, sid, opt["h2o_2pool"]["fbound"])
        qhf_eq = (1 - fb) * qh_eq
        qhb_eq = fb * qh_eq
        qco2_eq = co2_equilibrium(tk, pdrv, x[k, 2] + x[k, 3], phys["iso"])
        fs = pool_fraction(aw_eff, sid, opt["co2_2pool"]["fslow"])
        qcf_eq = (1 - fs) * qco2_eq
        qcs_eq = fs * qco2_eq
        kc, kh = rates(tk, sid, pdrv, aw_eff, phys["kin"], design)
        ks = clamp(
            opt["co2_2pool"]["kslow"]["ratio"] * kc,
            opt["co2_2pool"]["kslow"]["kmin"],
            opt["co2_2pool"]["kslow"]["kmax"],
        )
        kb = clamp(
            opt["h2o_2pool"]["kslow"]["ratio"] * kh,
            opt["h2o_2pool"]["kslow"]["kmin"],
            opt["h2o_2pool"]["kslow"]["kmax"],
        )
        kc_eff = kc
        ks_eff = ks
        if sid == 6:
            sat = max(
                1 - min(max(x[k, 0] + x[k, 1], 0.0) / max(qco2_eq, 1e-12), 1.0),
                opt["co2_sat"]["floor"],
            )
            fac = sat ** opt["co2_sat"]["alpha"]
            kc_eff *= fac
            ks_eff *= fac
        qcf = qcf_eq + (x[k, 0] - qcf_eq) * exp(-kc_eff * dt)
        qcs = qcs_eq + (x[k, 1] - qcs_eq) * exp(-ks_eff * dt)
        qhf = qhf_eq + (x[k, 2] - qhf_eq) * exp(-kh * dt)
        qhb = qhb_eq + (x[k, 3] - qhb_eq) * exp(-kb * dt)
        qcs = co2_slow_residual(
            qcs, qads_prev, tk, aw_eff, sid, opt["co2_2pool"]["resid"]
        )
        qhb = h2o_bound_residual(
            qhb, peak_h2o_bound, tk, ctrl, sid, opt["h2o_2pool"]["resid"]
        )
        qcf = max(qcf, 0.0)
        qcs = max(qcs, 0.0)
        qhf = clamp(qhf, 0.0, 80.0)
        qhb = clamp(qhb, 0.0, 80.0)
        qco2tot_next = qcf + qcs
        qinf = co2_capacity(tk, phys["iso"])
        if qco2tot_next > qinf:
            qcf = max(qcf - (qco2tot_next - qinf), 0.0)

        # Actual loading changes after kinetic corrections
        dq_co2_dt = ((qcf + qcs) - (x[k, 0] + x[k, 1])) / dt

        dq_h2o_dt = ((qhf + qhb) - (x[k, 2] + x[k, 3])) / dt

        # Adsorption/desorption heat [W]
        Q_ads_heat_W = adsorption_heat_rate(dq_co2_dt, dq_h2o_dt, design, phys)

        Q_ads_CO2_W = design["msorb_kg"] * phys["eq"]["dHads_CO2"] * dq_co2_dt

        Q_ads_H2O_W = design["msorb_kg"] * phys["eq"]["dHads_H2O"] * dq_h2o_dt

        z["Q_ads_CO2_W"][k] = Q_ads_CO2_W
        z["Q_ads_H2O_W"][k] = Q_ads_H2O_W

        # Gas concentrations at current timestep [mol/m³]
        c_total = pk / (R * tk)

        y_co2 = clamp(p_inv / pk, 0.0, 0.95)
        y_h2o = clamp(ph / pk, 0.0, 0.98)
        y_air = clamp(1.0 - y_co2 - y_h2o, 0.0, 1.0)

        c_co2 = y_co2 * c_total
        c_h2o = y_h2o * c_total
        c_air = y_air * c_total

        # Effective heat capacity [J/K]
        C_eff = eff_heat_capacity(c_co2, c_h2o, c_air, design, phys)
        # Heater hysteresis control
        T_regen_K = ctrl["T_regenC"] + 273.15
        T_h_K = equipment["T_h_C"] + 273.15

        u_h = heater_controller(
            tk, T_regen_K, sid, u_h_prev, hysteresis=equipment["hysteresis_h_K"]
        )

        # Actual heater thermal power [W]
        Q_heater_W = heater_power(
            tk, T_h_K, sid, equipment["UA_h_W_K"], equipment["Q_h_max_W"], u_h
        )

        # Store controller state for next timestep
        u_h_prev = u_h

        # Cooler target and utility temperatures [K]
        T_cool_target_K = amb["T_ambC"] + 5.0 + 273.15
        T_c_K = equipment["T_c_C"] + 273.15

        # Cooler hysteresis control
        u_c = cooler_controller(
            tk, T_cool_target_K, sid, u_c_prev, hysteresis=equipment["hysteresis_c_K"]
        )

        # Actual cooling thermal power [W]
        Q_cooler_W = cooler_power(
            tk, T_c_K, sid, equipment["UA_c_W_K"], equipment["Q_c_max_W"], u_c
        )

        # Update controller state
        u_c_prev = u_c

        # Heat loss to surroundings [W]
        T_amb_K = amb["T_ambC"] + 273.15

        Q_loss_W = heat_loss_rate(tk, T_amb_K)

        # Predict actual inlet flow during adsorption [m3/s]
        if sid == 6:
            Qlim_flow = opt["ptraj"]["Qmax_mult_Vg_dt"] * design["Vg_m3"] / dt

            qin_raw_heat = clamp(Q_ads, 0.0, Qlim_flow)

            if opt["ptraj"]["enable_smoothing"]:
                alpha_flow = 1.0 - exp(-dt / max(opt["ptraj"]["tau_flow_ramp_s"], 1e-9))

                qin_for_heat = Qin_filt + alpha_flow * (qin_raw_heat - Qin_filt)

            else:
                qin_for_heat = qin_raw_heat

        else:
            qin_for_heat = 0.0

        # Gas-flow sensible heat transfer [W]
        Q_flow_W = gas_flow_heat_rate(
            sid=sid,
            Q_in_m3_s=qin_for_heat,
            P_dry_Pa=pdry,
            P_h2o_Pa=ph,
            T_bed_K=tk,
            T_in_K=T_amb_K,
            y_co2_dry=amb["yCO2_air"],
            eq=phys["eq"],
        )

        # Dynamic thermal energy balance [K/s]
        dT_dt = dynamic_energy_balance(
            C_eff, Q_heater_W, Q_loss_W, Q_ads_heat_W, Q_cooler_W, Q_flow_W
        )

        # Update bed temperature [K]
        x[k + 1, 5] = tk + dt * dT_dt

        # Record heat transfer rates [W]
        z["Q_heater_W"][k] = Q_heater_W
        z["Q_ads_heat_W"][k] = Q_ads_heat_W
        z["Q_cooler_W"][k] = Q_cooler_W
        z["Q_flow_W"][k] = Q_flow_W
        z["Q_loss_W"][k] = Q_loss_W

        _, pdry1, _ = pressure_chain(
            float(drv["RH_cmd"][k + 1]),
            float(drv["TK_cmd"][k + 1]),
            x[k + 1, 5],
            float(drv["P_Pa"][k + 1]),
            opt,
        )
        nd = max(pdry * design["Vg_m3"] / (R * max(tk, 1.0)), 1e-12)
        nd_next = max(pdry1 * design["Vg_m3"] / (R * max(x[k + 1, 5], 1.0)), 1e-12)
        z["n_dry"][k] = nd
        # pressure inventory filtering; reset uses old n_dry then next under stage boundary
        tau_nd = opt["ptraj"]["tau_ndry_s"][names[int(sid_all[k + 1]) - 1]]
        if not np.isfinite(n_dry_filt_prev):
            n_dry_filt_prev = nd
        if (
            opt["ptraj"]["enable_ndry_filter"]
            and opt["ptraj"]["ndry_reset_on_stage_change"]
            and stage_changed
        ):
            nd_filt_next = nd_next
            n_dry_filt_prev = nd
        elif opt["ptraj"]["enable_ndry_filter"] and tau_nd > 0:
            nd_filt_next = (
                exp(-dt / tau_nd) * n_dry_filt_prev + (1 - exp(-dt / tau_nd)) * nd_next
            )
        else:
            nd_filt_next = nd_next
        ndot_target = (nd_filt_next - n_dry_filt_prev) / dt
        n_dry_filt_prev = nd_filt_next

        nd_sorb = design["msorb_kg"] * ((qcf + qcs) - (x[k, 0] + x[k, 1]))
        nd_sorb_dot = nd_sorb / dt
        # Figure 2 option enable_ndot_rate_limit=1.
        Qlim = opt["ptraj"]["Qmax_mult_Vg_dt"] * design["Vg_m3"] / dt
        c_in = max(pdry / (R * max(tk, 1.0)), 1e-12)
        c_out = max(nd / max(design["Vg_m3"], 1e-12), 1e-12)
        if opt["ptraj"]["enable_ndot_rate_limit"] and sid != 6:
            rate_limit = Qlim * (c_in if sid == 5 else c_out)
            req = ndot_target + nd_sorb_dot
            req = clamp(req, 0, rate_limit) if sid == 5 else clamp(req, -rate_limit, 0)
            ndot_target = req - nd_sorb_dot
        qin_raw, qout_raw, yin = flows_from_pressure(
            sid,
            Q_ads,
            amb["yCO2_air"],
            nd,
            ndot_target,
            nd_sorb_dot,
            pdry,
            tk,
            R,
            design["Vg_m3"],
        )
        qin_raw = clamp(qin_raw, 0, Qlim)
        qout_raw = clamp(qout_raw, 0, Qlim)
        if opt["ptraj"]["enable_smoothing"]:
            if (
                stage_changed
                and opt["ptraj"]["reset_flow_filters_on_stage_change"]
                and sid != 6
                and int(sid_all[k + 1]) != 6
            ):
                Qin_filt = Qout_filt = 0.0
            af = 1 - exp(-dt / max(opt["ptraj"]["tau_flow_ramp_s"], 1e-9))
            Qin_filt += af * (qin_raw - Qin_filt)
            Qout_filt += af * (qout_raw - Qout_filt)
            qin = Qin_filt
            qout = Qout_filt
        else:
            qin = qin_raw
            qout = qout_raw
        # Verify predicted inlet flow matches actual flow
        if sid == 6:
            assert np.isclose(
                qin_for_heat, qin, rtol=1e-10, atol=1e-10
            ), f"Inlet flow mismatch at timestep {k}"
        yc = update_y(
            x[k, 4],
            nd,
            nd_next,
            qin,
            qout,
            yin,
            pdry,
            tk,
            design["Vg_m3"],
            R,
            nd_sorb,
            dt,
            opt["yco2"],
        )
        x[k + 1, :5] = [qcf, qcs, qhf, qhb, yc]

        peak_h2o_bound = max(peak_h2o_bound, qhb)
        z["qCO2eq"][k] = qco2_eq
        z["pCO2_inv"][k] = p_inv
        z["pCO2_drive"][k] = pdrv
        z["p_h2o"][k] = ph
        z["pdry"][k] = pdry
        z["aw_gas"][k] = aw
        z["aw_eff"][k] = aw_eff
        z["q_in"][k] = qin
        z["q_out"][k] = qout
        z["y_in"][k] = yin
        z["kco2_fast"][k] = kc
        if opt["ptraj"]["log_ndry_residual"]:
            c_out = max(nd / max(design["Vg_m3"], 1e-12), 1e-12)
            z["nd_res"][k] = ndot_target - (qin * c_in - qout * c_out - nd_sorb_dot)
    # Last point diagnostics
    for key in ("qCO2eq", "kco2_fast", "nd_res"):
        z[key][-1] = z[key][-2]
    p_h, p_d, aw = pressure_chain(
        float(drv["RH_cmd"][-1]),
        float(drv["TK_cmd"][-1]),
        x[-1, 5],
        float(drv["P_Pa"][-1]),
        opt,
    )
    z["pdry"][-1] = p_d
    z["p_h2o"][-1] = p_h
    z["aw_gas"][-1] = aw
    z["aw_eff"][-1] = effective_aw(
        int(sid_all[-1]), x[-1, 5], aw, x[-1, 2], x[-1, 3], phys["iso"], opt
    )
    z["pCO2_inv"][-1] = x[-1, 4] * p_d
    z["pCO2_drive"][-1] = select_p_co2(
        int(sid_all[-1]), z["pCO2_inv"][-1], amb["yCO2_air"] * p_d, opt
    )
    z["n_dry"][-1] = p_d * design["Vg_m3"] / (R * x[-1, 5])
    qco2 = x[:, 0] + x[:, 1]
    qh2o = x[:, 2] + x[:, 3]
    ia = np.where(sid_all == 6)[0]
    ides = np.where(sid_all == 3)[0]
    qads = qco2[ia[-1]] if len(ia) else qco2[-1]
    qdes = qco2[ides[-1]] if len(ides) else qco2[0]
    delta = max(qads - qdes, 0.0)
    n_work = design["msorb_kg"] * delta
    k0 = ia[0]
    ninit = max(
        (
            x[k0, 4]
            - (
                opt["yco2"]["y_floor"]
                if opt.get("ads_supply", {}).get("keep_y_floor_inventory", True)
                else 0.0
            )
        )
        * z["n_dry"][k0],
        0.0,
    )
    p_in_ads = amb["yCO2_air"] * z["pdry"][ia]
    nin = float(np.trapezoid(Q_ads * p_in_ads / (R * x[ia, 5]), t[ia]))
    ncap = min(n_work, ninit + nin)
    cyc = {
        "qCO2_ads_end": qads,
        "qCO2_des_end": qdes,
        "dqCO2_raw": delta,
        "nCO2_working_raw": n_work,
        "nCO2_supply_max": ninit + nin,
        "nCO2_captured_eff": ncap,
        "cap_ratio": ncap / max(n_work, 1e-12),
    }
    prof = {
        "time_s": t,
        "qCO2_total": qco2,
        "qH2O_total": qh2o,
        "qCO2_fast": x[:, 0],
        "qCO2_slow": x[:, 1],
        "qH2O_free": x[:, 2],
        "qH2O_bound": x[:, 3],
        "yCO2g": x[:, 4],
        "TK_bed": x[:, 5],
        "TC_bed": x[:, 5] - 273.15,
        "TK_cmd": drv["TK_cmd"],
        "RH_cmd": drv["RH_cmd"],
    }
    prof.update(z)
    return {"drv": drv, "prof": prof, "x_end": x[-1].copy(), "states": x, "cyc": cyc}


def fig2_simulations(ctrl, design, phys, opt, make_driver, extract_driver, ambient_fn):
    """Exact Figure 2 memory protocol: precondition, moderate seed,
    then one adapt cycle + one plotted cycle for each ambient."""
    amb_base = ambient_fn(20.0, 0.5)
    drv_base = make_driver(amb_base, ctrl, opt)
    drv_ads = extract_driver(drv_base, opt)
    clean = np.array(
        [0.0, 0.0, 0.0, 0.0, amb_base["yCO2_air"], amb_base["T_ambC"] + 273.15]
    )
    pre = forward_cycle(drv_ads, amb_base, ctrl, design, phys, opt, clean)
    x0 = pre["x_end"].copy()
    x0[4] = amb_base["yCO2_air"]
    x0[5] = amb_base["T_ambC"] + 273.15
    baseline = forward_cycle(drv_base, amb_base, ctrl, design, phys, opt, x0)
    memory = baseline["x_end"].copy()
    memory[4] = amb_base["yCO2_air"]
    memory[5] = amb_base["T_ambC"] + 273.15
    results = []
    for name, T, RH in [
        ("Cold and dry", 4.0, 0.2),
        ("Moderate", 20.0, 0.5),
        ("Hot and humid", 40.0, 0.9),
    ]:
        amb = ambient_fn(T, RH)
        drv = make_driver(amb, ctrl, opt)
        x0 = memory.copy()
        x0[4] = amb["yCO2_air"]
        x0[5] = T + 273.15
        adapt = forward_cycle(drv, amb, ctrl, design, phys, opt, x0)
        xcase = adapt["x_end"].copy()
        xcase[4] = amb["yCO2_air"]
        xcase[5] = T + 273.15
        sim = forward_cycle(drv, amb, ctrl, design, phys, opt, xcase)
        results.append((name, amb, adapt, sim))
    return baseline, results
