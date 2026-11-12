"""Gas pressure/flow closure, core MATLAB lines 1782-1938, 2033-2057."""

from isotherms import h2o_equilibrium, psat_h2o


def clamp(x, a, b):
    return min(max(x, a), b)


def pressure_chain(rh, t_cmd, t_bed, p_total, opt):
    rh = clamp(rh, 0, 0.999)
    p_total = max(p_total, 200.0)
    if opt["aw"]["chain_mode"] == "bed_consistent":
        p_h2o = clamp(rh * psat_h2o(t_bed), 0, 0.98 * p_total)
    else:
        p_h2o = clamp(rh * psat_h2o(t_cmd), 0, 0.98 * p_total)
    p_dry = max(p_total - p_h2o, 1.0)
    a_w = clamp(p_h2o / max(psat_h2o(t_bed), 1.0), 0.0, opt["aw"]["max"])
    return p_h2o, p_dry, a_w


def invert_gab(T, q, iso, aw_max=0.999):
    q = max(q, 0.0)
    if q <= h2o_equilibrium(T, 0, iso, aw_max):
        return 0.0
    if q >= h2o_equilibrium(T, aw_max, iso, aw_max):
        return aw_max
    lo = 0.0
    hi = aw_max
    for _ in range(18):
        mid = 0.5 * (lo + hi)
        if h2o_equilibrium(T, mid, iso, aw_max) < q:
            lo = mid
        else:
            hi = mid
    return clamp(0.5 * (lo + hi), 0.0, aw_max)


def effective_aw(sid, T, aw_g, q_h2o_free, q_h2o_bound, iso, opt):
    maxaw = opt["aw"]["max"]
    aw_g = clamp(aw_g, 0.0, maxaw)
    q_free = max(q_h2o_free, 0.0)
    q_bound = max(q_h2o_bound, 0.0)
    qtot = q_free + q_bound
    aw_gab = invert_gab(T, qtot, iso, maxaw)
    mode = opt["aw"]["mode_for_co2_effects"].lower()
    if mode == "use_aw_gas":
        aw = aw_g
    else:
        if opt["h2o_2pool"]["enable"] and opt["aw"].get("co2_boundmask", {}).get(
            "enable", False
        ):
            mask = opt["aw"].get("co2_boundmask", {})
            free_frac = (q_free / max(qtot + 1e-9, 1e-9)) ** max(
                mask.get("power", 1.0), 0.25
            )
            aw = aw_g + free_frac * (aw_gab - aw_g)
        else:
            aw = aw_gab
    aw = clamp(aw, 0.0, maxaw)
    if sid == 6 and opt.get("co2_awfloor", {}).get("enable", False):
        prm = opt["co2_awfloor"]
        if aw_g < prm.get("apply_below", 0.25):
            aw = clamp(max(aw, prm.get("value", 0.0)), 0, maxaw)
    return aw


def select_p_co2(sid, p_inv, p_in, opt):
    cfg = opt["pco2"]
    mode = cfg["ads_mode"] if sid == 6 else cfg["regen_mode"]
    if mode == "inlet":
        p = p_in
    elif mode == "blend":
        p = (
            cfg["blend_alpha_inventory"] * p_inv
            + (1 - cfg["blend_alpha_inventory"]) * p_in
        )
    else:
        p = p_inv
    if sid == 6:
        p = max(p, cfg.get("ads_p_floor_frac", 1.0) * p_in)
    return max(p, 0.0)


def flows_from_pressure(
    sid, q_ads, y_air, n_dry, ndot_target, dn_sorb_dot, pdry, T, R, vg
):
    # Reproduced from MATLAB flows_from_ptraj at 1880-1909.
    c_in = max(pdry / (R * max(T, 1.0)), 1e-12)
    c_out = max(n_dry / max(vg, 1e-12), 1e-12)
    if sid == 6:
        qin = q_ads
        qout = q_ads
        yin = y_air
    elif sid == 5:
        qin = max((ndot_target + dn_sorb_dot) / c_in, 0.0)
        qout = 0.0
        yin = y_air
    else:
        qin = 0.0
        qout = max(-(ndot_target + dn_sorb_dot) / c_out, 0.0)
        yin = 0.0
    return qin, qout, yin


def update_y(y0, n_dry, n_next, qin, qout, yin, pdry, T, vg, R, dn_co2, dt, cfg):
    # Reproduced from MATLAB update_yCO2g_from_flows, see source comments.
    y0 = clamp(y0, cfg["y_floor"], cfg["y_max"])
    c_in_co2 = yin * pdry / (R * max(T, 1.0))
    nco2_prev = y0 * max(n_dry, 1e-12)
    rhs = nco2_prev + qin * c_in_co2 * dt - dn_co2
    nco2_mb = rhs / (1 + qout * dt / max(vg, 1e-12))
    nco2_mb = max(nco2_mb, cfg["y_floor"] * max(n_dry, 1e-12))
    ycalc = clamp(nco2_mb / max(n_next, 1e-12), cfg["y_floor"], cfg["y_max"])
    if cfg["dy_cap"]["enable"]:
        ycalc = clamp(ycalc, y0 - cfg["dy_cap"]["dy_max"], y0 + cfg["dy_cap"]["dy_max"])
    return clamp(ycalc, cfg["y_floor"], cfg["y_max"])
