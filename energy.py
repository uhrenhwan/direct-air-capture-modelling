"""Figure 2 enabled energy accounting from MATLAB lines 2312-2671.

No optional correction layers enabled in the Figure 2 script. The only
active dry-condition adjustment is vacuum aw penalty.
"""

from math import log

import numpy as np

from cycle_driver import ambient_pressure_bar
from isotherms import psat_h2o
from simulator import flow_ads


def clamp(x, lo, hi):
    return min(max(x, lo), hi)


def energy_metrics(sim, amb, ctrl, design, phys, opt):
    """Calculate the reference Figure 2 electricity accounting and annual productivity.

    This accounting belongs to the prescribed-temperature baseline.
    Thermal-extension heater SEC is calculated separately in optimisation.py."""
    drv = sim["drv"]
    p = sim["prof"]
    cyc = sim["cyc"]
    eq = phys["eq"]
    R = eq["R"]
    dt = opt["dt_s"]
    t = drv["t_s"]
    tcycle = t[-1]
    cycles_y = 365 * 86400.0 / tcycle
    nco2 = max(cyc["nCO2_captured_eff"], 0.0)
    kg = nco2 * eq["MW_CO2"]
    tonne = max(kg / 1000.0, 1e-12)
    Prod = max(kg * cycles_y / 1000.0, 0.0)
    results = dict(
        Prod_tCO2_collector_yr=Prod, Prod_tCO2_m3bed_yr=Prod / design["Vbed_m3"]
    )
    if kg < 1e-12:
        results.update(Elec_total_MWh_tCO2=float("nan"))
        return results
    idx = np.where((drv["stage_id"] == 2) | (drv["stage_id"] == 3))[0]
    ti = t[idx]
    T = p["TK_bed"][idx]
    dT = np.gradient(T, dt)
    dqco2 = np.gradient(p["qCO2_total"][idx], dt)
    dqh2o = np.gradient(p["qH2O_total"][idx], dt)
    ptot = drv["P_Pa"][idx]
    ctot = ptot / (R * np.maximum(T, 1.0))
    yH = np.clip(p["p_h2o"][idx] / np.maximum(ptot, 1.0), 0, 0.98)
    yC = np.clip(p["pCO2_inv"][idx] / np.maximum(ptot, 1.0), 0, 0.95)
    yA = np.clip(1 - yH - yC, 0, 1)
    gasCp = (
        design["eps"]
        * ctot
        * (yC * eq["Cp_CO2_g"] + yH * eq["Cp_H2O_g"] + yA * eq["Cp_air_g"])
    )
    solidCp = phys["bed"]["rb"] * eq["Cp_s"]
    sens = (gasCp + solidCp) * dT
    rxn = -phys["bed"]["rb"] * (eq["dHads_CO2"] * dqco2 + eq["dHads_H2O"] * dqh2o)
    Qbed = design["Vbed_m3"] * float(np.trapezoid(np.maximum(sens + rxn, 0.0), ti))
    E_bed_th = Qbed / 3.6e9 / tonne
    nsteam = eq["r_steam_CO2"] * nco2
    msteam = nsteam * eq["MW_H2O"]
    E_steam_sens_th = (
        msteam * eq["Cp_liq"] * max(eq["TsatC"] - amb["T_ambC"], 0.0) / 3.6e9 / tonne
    )
    E_steam_gen_th = msteam * eq["hvap"] / 3.6e9 / tonne
    # In original raw PPI, absent phys.ctrl -> T heating from max bed T.
    TH = max(T) + eq["dTmin_K"]
    TC = max(amb["T_ambC"] + 273.15 - eq["dTmin_K"], 273.15)
    TC = min(TC, TH - 1e-3)
    COP = clamp(eq["eta2nd"] * TH / max(TH - TC, 1e-6), 1.2, 6.0)
    E_bed = E_bed_th / COP
    E_steam_sens = E_steam_sens_th / COP
    E_steam_gen = E_steam_gen_th / COP
    # dry-gas vacuum compression, discrete sum just like MATLAB.
    patm = ambient_pressure_bar(amb["altitude_m"]) * 1e5
    ph_amb = clamp(amb["RH_in"] * psat_h2o(amb["T_ambC"] + 273.15), 0.0, 0.98 * patm)
    patm_dry = max(patm - ph_amb, 1.0)
    pressure = p["pdry"][:-1]
    tk = p["TK_bed"][:-1]
    mask = (pressure < patm_dry) & (p["q_out"][:-1] > 0)
    n_dot = (pressure[mask] / (R * np.maximum(tk[mask], 250.0))) * p["q_out"][:-1][mask]
    w_mol = (
        R
        * np.maximum(tk[mask], 250.0)
        * np.log(patm_dry / np.maximum(pressure[mask], 1.0))
    )
    W_vac = float(np.sum(n_dot * w_mol * dt))
    E_vac = W_vac / max(eq["eta_comp"], 1e-6) / 3.6e9 / tonne
    aw_reg = np.clip(
        p["p_h2o"][idx] / np.maximum([psat_h2o(float(v)) for v in T], 1.0), 0.0, 0.999
    )
    aw_min = float(np.min(aw_reg))
    if aw_min < opt["energy"]["vac_aw_ref"]:
        E_vac *= 1 + opt["energy"]["vac_aw_gain"] * (
            opt["energy"]["vac_aw_ref"] - aw_min
        )
    W_prod = nco2 * R * eq["T_comp_K"] * log(eq["P_CO2_out_bar"] * 1e5 / max(patm, 1.0))
    E_prod = W_prod / max(eq["eta_comp"], 1e-6) / 3.6e9 / tonne
    # Ergun contactor pressure drop; MATLAB uses A_flow even though Qads = u*Acontact.
    Q_ads = flow_ads(ctrl, design, opt)
    adsidx = np.where(drv["stage_id"] == 6)[0]
    k0 = adsidx[0]
    Ts = p["TK_bed"][k0]
    Pdry = max(p["pdry"][k0], 1.0)
    mu0 = 1.716e-5
    S = 111.0
    mu = mu0 * (Ts / 273.15) ** 1.5 * (273.15 + S) / (Ts + S)
    mu = max(mu, 8e-6)
    rho = max(Pdry * eq["MW_air"] / (R * max(Ts, 1.0)), 0.2)
    u = Q_ads / max(design["A_flow_m2"], 1e-12)
    eps = design["eps"]
    dp = design["dp_m"]
    term1 = 150 * (1 - eps) ** 2 / eps**3 * mu * u / dp**2
    term2 = 1.75 * (1 - eps) / eps**3 * rho * u * u / dp
    dP = (
        design["N_series_pass"] * (term1 + term2) * design["L_flow_m"]
        + design["dP_misc_Pa"]
    )
    t_ads = np.sum(drv["stage_id"][:-1] == 6) * dt
    E_fan = Q_ads * t_ads * dP / max(eq["eta_fan"], 1e-6) / 3.6e9 / tonne
    eheat = E_bed + E_steam_sens + E_steam_gen
    etot = eheat + E_vac + E_prod + E_fan
    results.update(
        COP=COP,
        Elec_total_MWh_tCO2=etot,
        breakdown={
            "E_bed": E_bed,
            "E_steam_sens": E_steam_sens,
            "E_steam_gen": E_steam_gen,
            "E_heating": eheat,
            "E_comp_vac": E_vac,
            "E_comp_prod": E_prod,
            "E_comp": E_vac + E_prod,
            "E_fan": E_fan,
            "E_parasitic": 0.0,
            "E_idle": 0.0,
        },
        audit={"aw_min_reg": aw_min, "nCO2_captured_mol": nco2, "mCO2_kg": kg},
    )
    return results
