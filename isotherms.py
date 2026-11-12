"""Water GAB and coupled CO2-water equilibrium relations from the reference model."""

from math import exp

R = 8.314462618  # MATLAB universal gas constant [J/(mol K)]

# saturation vapour pressure of water


def psat_h2o(T_K):
    T_C = T_K - 273.15
    # Arden Buck equation (empirical correlation)
    p_sat_hPa = 6.1121 * exp((18.678 - T_C / 234.5) * (T_C / (257.14 + T_C)))
    p_sat_Pa = max(p_sat_hPa * 100, 1.0)

    return p_sat_Pa


def h2o_vapour_pressure(T_K, RH, P):
    p_sat = psat_h2o(T_K)
    p_h2o = RH * p_sat
    p_h2o = min(max(p_h2o, 0), 0.98 * P)
    return p_h2o


def water_activity(T_K, RH, P, a_w_max=0.999):
    p_h2o = h2o_vapour_pressure(T_K, RH, P)
    p_sat = psat_h2o(T_K)

    a_w = p_h2o / p_sat
    a_w = min(max(a_w, 0.0), a_w_max)

    return a_w


def h2o_equilibrium(T_K, a_w, iso, a_w_max=0.999):
    R = 8.314462618  # J/mol/K; match MATLAB exactly
    a_w = min(max(a_w, 0), a_w_max)

    gab = iso["gab"]

    # eqn from Young et al. The impact of binary water-CO2 isotherm models on the optimal performance of sorbent-based direct air capture processes. 2021
    # experimentally derived equations
    # E1: first molecular layer adsorption heat; E2m9: 2nd-9th molecular layer adsorption heat; E10p: 10th and higher molecular layer adsorption heat
    # C,D,F,G: experimentally fitted constants; empirical correlation

    E10p = -44.38 * T_K + 57220
    E1 = gab["C"] - exp(gab["D"] * T_K)
    E2m9 = gab["F"] + gab["G"] * T_K

    C_gab = exp((E1 - E10p) / (R * T_K))
    K_gab = exp((E2m9 - E10p) / (R * T_K))

    denominator = (1 - K_gab * a_w) * (1 + (C_gab - 1) * K_gab * a_w)
    denominator = max(denominator, 1e-12)

    q_h2o = gab["qm"] * (K_gab * C_gab * a_w) / denominator

    return max(q_h2o, 0.0)


#######################################################
# ------CO2 adsorption equilibrium model------#


# maximum co2 adsorption capacity
def co2_capacity(T_K, iso):
    mech = iso["mech"]

    q_inf = mech["S_cap"] * mech["qinf0"] * exp(mech["chi"] * (1 - T_K / mech["T0"]))

    return q_inf


# H2O CO2 adsorption site blocking effect


def co2_water_effects(q_h2o, iso):
    mech = iso["mech"]

    q_h2o_eff = max(q_h2o, 1e-9)

    w = exp(-mech["Aexp"] / q_h2o_eff)

    arg_blk = (mech["k_blk"] * q_h2o_eff) ** mech["n_blk"]

    f_blk = mech["fblk_max"] * (1 - exp(-arg_blk))

    f_blk = min(max(f_blk, 0.0), mech["fblk_max"])

    phi_avail = mech["phi_max"] - f_blk

    phi = mech["phi_dry"] + (phi_avail - mech["phi_dry"]) * w

    phi = min(max(phi, 0.0), mech["phi_max"])

    return w, f_blk, phi


def co2_affinity(T_K, w, iso):
    mech = iso["mech"]
    R = 8.314462618

    dH_avg = (1 - w) * mech["DH_dry"] + w * mech["DH_wet"]

    b_T = mech["b0"] * exp(dH_avg / (R * T_K))

    return dH_avg, b_T


def co2_heterogeneity(T_K, iso):
    mech = iso["mech"]

    tau_T = mech["tau0"] + mech["alpha"] * (1 - mech["T0"] / T_K)

    tau_T = max(tau_T, 1e-6)

    return tau_T


def co2_equilibrium(T_K, p_co2, q_h2o, iso):
    mech = iso["mech"]

    # Step 1: Maximum CO2 adsorption capacity
    q_inf = co2_capacity(T_K, iso)

    # Step 2: Water-induced site blocking
    w, f_blk, phi = co2_water_effects(q_h2o, iso)

    # Step 3: Temperature-dependent adsorption affinity
    dH_avg, b_T = co2_affinity(T_K, w, iso)

    # Step 4: Adsorption heterogeneity
    tau_T = co2_heterogeneity(T_K, iso)

    # Step 5: Dimensionless adsorption term
    bp = max(b_T * max(p_co2, 0.0), 0.0)

    denominator = (1 + bp**tau_T) ** (1 / tau_T)

    # Step 6: Equilibrium CO2 loading
    q_co2_eq = (phi / mech["phi_dry"]) * q_inf * (bp / max(denominator, 1e-12))

    # Physical bounds
    q_co2_eq = min(max(q_co2_eq, 0.0), q_inf)

    return q_co2_eq


# Functional checks can be run as a script, without printing during imports.
if __name__ == "__main__":
    from parameters import phys

    T_K, RH, P, p_co2 = 293.15, 0.5, 101325, 40.0
    iso = phys["iso"]
    a_w = water_activity(T_K, RH, P)
    q_h2o_eq = h2o_equilibrium(T_K, a_w, iso)
    q_co2_eq = co2_equilibrium(T_K, p_co2, q_h2o_eq, iso)
    print("a_w =", a_w)
    print("q_h2o_eq =", q_h2o_eq, "mol/kg")
    print("q_co2_eq =", q_co2_eq, "mol/kg")
