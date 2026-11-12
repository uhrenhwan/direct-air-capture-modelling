"""LDF rates and two-pool equilibrium fractions, MATLAB core lines 2062-2166, 2221-2282."""

from math import exp

from isotherms import R


def clamp(x, lo, hi):
    return min(max(x, lo), hi)


def pool_fraction(activity, stage, p):
    a = "ads" if stage == 6 else "regen"
    return clamp(
        p["f0_" + a] + p["k_aw_" + ("ads" if stage == 6 else "reg")] * activity,
        p["f_min"],
        p["f_max"],
    )


def rates(T, stage, p_co2, aw_eff, kin, design):
    kc = kin["kCO2_ref"] * exp(-kin["Ea_CO2"] / R * (1 / T - 1 / kin["Tref"]))
    kh = kin["kH2O_ref"] * exp(-kin["Ea_H2O"] / R * (1 / T - 1 / kin["Tref"]))
    kc *= kin["mCO2"][stage - 1]
    kh *= kin["mH2O"][stage - 1]
    kc *= 1 + kin["beta"] * (p_co2 / kin["p_ref_Pa"]) ** kin["n_exp"]
    if stage == 6:
        kc *= kin.get("kCO2_ads_mult", 1.0)
    if stage == 6:
        gate = kin.get("gate_lowaw", dict(aw_hi=0.22, aw_lo=0.14, fmax=4.0, power=2.0))
        if gate["aw_hi"] > gate["aw_lo"] and gate["fmax"] > 1:
            u = clamp(
                (gate["aw_hi"] - aw_eff) / (gate["aw_hi"] - gate["aw_lo"]), 0.0, 1.0
            )
            kc *= 1 + (gate["fmax"] - 1) * u ** gate["power"]
    # MATLAB caps inside kinetics_LDF first, THEN multiplies sqrt(A_mult) again.
    kc = clamp(kc, 1e-8, 0.05) * (design["A_mult"] ** 0.5)
    kh = clamp(kh, 1e-8, 0.10) * (design["A_mult"] ** 0.5)
    return kc, kh


def co2_slow_residual(q_slow, q_ads_peak, T, aw, stage, p):
    if stage not in p["apply_in_stage_ids"]:
        return q_slow
    frac = clamp(
        p["f_res"]
        + p["k_wet"] * max(aw, 0.0) ** p["n_wet"]
        + p["k_T"] * max(p["Tref_K"] / max(T, 1.0) - 1, 0.0),
        p["f_min"],
        p["f_max"],
    )
    qmin = max(p["q_abs_min"], frac * max(q_ads_peak, 1e-12))
    return qmin if q_slow < qmin * (1 - p["tol_frac"]) else q_slow


def h2o_bound_residual(q_bound, q_peak, T, ctrl, stage, p):
    if stage not in p["apply_in_stage_ids"]:
        return q_bound
    t_reg = ctrl["time_heat_s"] + ctrl["time_des_s"]
    frac = clamp(
        p["f0"]
        + p["k_T"] * max(p["Tref_K"] / max(T, 1.0) - 1, 0.0)
        + p["k_time"] * max(p["tref_s"] / max(t_reg, 1.0) - 1, 0.0),
        0.0,
        p["f_max"],
    )
    qmin = max(p["q_abs_min"], frac * max(q_peak, 1e-12))
    return qmin if q_bound < qmin * (1 - p["tol_frac"]) else q_bound
