"""Bed heat capacity, sorption heat and capacity-limited heater/cooler control."""

from parameters import design, phys

R = phys["eq"]["R"]  # J/(mol K)

# Hardware heat capacity (assumed)
C_hardware = 0.0  # J/K


def gas_concentrations(T_bed, P, p_co2, p_h2o):
    """Return CO2, water and air concentrations [mol/m3] from pressures [Pa] and K."""
    c_tot = P / (R * T_bed)

    y_co2 = min(max(p_co2 / max(P, 1.0), 0.0), 0.95)
    y_h2o = min(max(p_h2o / max(P, 1.0), 0.0), 0.98)
    y_air = min(max(1.0 - y_co2 - y_h2o, 0.0), 1.0)

    c_co2 = y_co2 * c_tot
    c_h2o = y_h2o * c_tot
    c_air = y_air * c_tot

    return c_co2, c_h2o, c_air


def eff_heat_capacity(c_co2, c_h2o, c_air, design_params, phys_params):
    """Return the combined sorbent, void-gas and hardware heat capacity [J/K]."""
    eq = phys_params["eq"]

    C_sorbent = design_params["msorb_kg"] * eq["Cp_s"]

    C_gas = (
        design_params["eps"]
        * design_params["Vbed_m3"]
        * (c_co2 * eq["Cp_CO2_g"] + c_h2o * eq["Cp_H2O_g"] + c_air * eq["Cp_air_g"])
    )

    C_eff = C_sorbent + C_gas + C_hardware

    return C_eff


def adsorption_heat_rate(dq_co2_dt, dq_h2o_dt, design_params=None, phys_params=None):
    # Compatibility with the original simulator.py
    """Return sorption heat [W]; positive loading rates release heat into the bed."""
    d = design if design_params is None else design_params
    p = phys if phys_params is None else phys_params

    Q_ads_heat_W = d["msorb_kg"] * (
        p["eq"]["dHads_CO2"] * dq_co2_dt + p["eq"]["dHads_H2O"] * dq_h2o_dt
    )

    return Q_ads_heat_W


# Overall heat loss coefficient [W/K]
UA = 0.0  # Initial adiabatic assumption


def heat_loss_rate(T_bed, T_amb):
    Q_loss = UA * (T_bed - T_amb)

    return Q_loss  # W


def heater_controller(T_bed, T_regen, stage_id, u_h_prev, hysteresis=1.0):
    # Heater operates only during heating and desorption
    """Switch heating on below the lower threshold and retain state in the deadband."""
    if stage_id not in (2, 3):
        return 0

    # Hysteresis On-Off Control
    if T_bed <= T_regen - hysteresis:
        u_h = 1

    elif T_bed >= T_regen:
        u_h = 0

    else:
        u_h = u_h_prev

    return u_h


def heater_power(T_bed, T_h, stage_id, UA_h, Q_h_max, u_h):
    if stage_id not in (2, 3):
        return 0.0

    Q_transfer = UA_h * max(T_h - T_bed, 0.0)

    Q_heater = u_h * min(Q_h_max, Q_transfer)

    return Q_heater


def cooler_controller(T_bed, T_target, stage_id, u_c_prev, hysteresis=1.0):
    """Switch cooling on above the upper threshold and retain state in the deadband."""
    if stage_id != 4:
        return 0

    if T_bed >= T_target + hysteresis:
        return 1

    elif T_bed <= T_target:
        return 0

    else:
        return u_c_prev


def cooler_power(T_bed, T_c, stage_id, UA_c, Q_c_max, u_c):
    if stage_id != 4:
        return 0.0

    Q_transfer = UA_c * max(T_bed - T_c, 0.0)

    Q_cooler = u_c * min(Q_c_max, Q_transfer)

    return Q_cooler


def dynamic_energy_balance(C_eff, Q_heater, Q_loss, Q_ads, Q_cooler, Q_flow=0.0):
    """Return dT/dt [K/s]; cooling and heat loss are positive outward duties."""
    dT_dt = (Q_heater - Q_loss + Q_ads - Q_cooler + Q_flow) / C_eff

    return dT_dt


def gas_flow_heat_rate(
    sid, Q_in_m3_s, P_dry_Pa, P_h2o_Pa, T_bed_K, T_in_K, y_co2_dry, eq
):
    """
    Sensible heat transfer by inlet/outlet gas [W].
    Applied during adsorption (Stage 6).

    Assumption: T_out = T_bed.
    Q_in follows the original simulator's bed-state volume convention.
    """

    if sid != 6:
        return 0.0

    R = eq["R"]

    # Inlet dry-gas molar flow [mol/s]
    n_dry = Q_in_m3_s * P_dry_Pa / (R * T_bed_K)

    # Inlet water-vapour molar flow [mol/s]
    n_h2o = n_dry * P_h2o_Pa / max(P_dry_Pa, 1e-12)

    # Dry gas mixture heat capacity [J/mol/K]
    Cp_dry = y_co2_dry * eq["Cp_CO2_g"] + (1.0 - y_co2_dry) * eq["Cp_air_g"]

    # Total gas heat-capacity flow rate [W/K]
    C_flow = n_dry * Cp_dry + n_h2o * eq["Cp_H2O_g"]

    return C_flow * (T_in_K - T_bed_K)
