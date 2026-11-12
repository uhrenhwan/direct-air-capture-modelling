"""Estimate sensible-only heating and cooling duties and utility conductances."""

import numpy as np

from parameters import amb_base, ctrl, design, phys


def main():
    # Existing DAC parameters
    m_s = design["msorb_kg"]
    Cp_s = phys["eq"]["Cp_s"]

    T_initial = amb_base["T_ambC"]  # °C
    T_target = ctrl["T_regenC"]  # °C
    t_heat = ctrl["time_heat_s"]  # s

    # Preliminary design assumptions
    T_hot = 140.0  # °C, illustrative heating utility
    C_eff = m_s * Cp_s  # J/K, sensible-only estimate

    # Sensible heating requirement
    Q_sensible = C_eff * (T_target - T_initial)

    # Required heat transfer conductance
    UA_h = (C_eff / t_heat) * np.log((T_hot - T_initial) / (T_hot - T_target))

    # Heating power requirements
    Q_h_initial = UA_h * (T_hot - T_initial)
    Q_h_average = Q_sensible / t_heat

    print(f"Sensible heating duty: {Q_sensible / 1e6:.2f} MJ")
    print(f"Required UA_h: {UA_h:.2f} W/K")
    print(f"Initial heating power: {Q_h_initial / 1e3:.2f} kW")
    print(f"Average heating power: {Q_h_average / 1e3:.2f} kW")

    # Preliminary Cooler Sizing

    T_cool_initial = ctrl["T_regenC"]  # 100°C
    T_cool_target = amb_base["T_ambC"] + 5  # 25°C
    T_cool_medium = 15.0  # °C, assumed
    t_cool = ctrl["time_cool_s"]  # 400 s

    # Check temperature feasibility
    if not T_cool_medium < T_cool_target < T_cool_initial:
        raise ValueError("Invalid cooling temperature conditions")

    # Sensible cooling requirement [J]
    Q_c_sensible = C_eff * (T_cool_initial - T_cool_target)

    # Required heat transfer conductance [W/K]
    UA_c = (C_eff / t_cool) * np.log(
        (T_cool_initial - T_cool_medium) / (T_cool_target - T_cool_medium)
    )

    # Cooling power [W]
    Q_c_initial = UA_c * (T_cool_initial - T_cool_medium)
    Q_c_average = Q_c_sensible / t_cool

    print(f"Sensible cooling duty: {Q_c_sensible / 1e6:.2f} MJ")
    print(f"Required UA_c: {UA_c:.2f} W/K")
    print(f"Initial cooling power: {Q_c_initial / 1e3:.2f} kW")
    print(f"Average cooling power: {Q_c_average / 1e3:.2f} kW")


if __name__ == "__main__":
    main()
