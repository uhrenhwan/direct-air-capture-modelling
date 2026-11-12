"""Run five thermal cycles, carrying the final bed state into the next cycle."""

import numpy as np

from cycle_driver import make_driver
from parameters import ambient, baseline_parameters
from simulator_thermal import forward_cycle


def main():
    ctrl, design, phys, opt = baseline_parameters()

    amb = ambient(20.0, 0.50)
    drv = make_driver(amb, ctrl, opt)

    x_init = None

    for cycle in range(1, 6):
        result = forward_cycle(drv, amb, ctrl, design, phys, opt, x_init=x_init)

        # Carry final state into the next cycle
        x_init = result["x_end"].copy()

        T_final = result["prof"]["TC_bed"][-1]
        n_CO2 = result["cyc"]["nCO2_captured_eff"]

        print(
            f"Cycle {cycle}: "
            f"T_final = {T_final:.2f} °C, "
            f"CO2 captured = {n_CO2:.2f} mol"
        )

    # Use the final cycle
    prof = result["prof"]
    dt = opt["dt_s"]

    # CO2 captured [kg/cycle]
    n_co2 = result["cyc"]["nCO2_captured_eff"]
    m_co2_kg = n_co2 * phys["eq"]["MW_CO2"]

    # Cycle duration [h]
    cycle_time_h = (drv["t_s"][-1] - drv["t_s"][0]) / 3600

    # Thermal energy [kWh_th/cycle]
    E_heater = np.sum(prof["Q_heater_W"][:-1]) * dt / 3.6e6
    E_cooler = np.sum(prof["Q_cooler_W"][:-1]) * dt / 3.6e6

    # Performance metrics
    productivity = m_co2_kg / cycle_time_h
    SEC_thermal = E_heater / m_co2_kg  # Specific energy consumption

    print(f"CO2 captured: {m_co2_kg:.3f} kg/cycle")
    print(f"Productivity: {productivity:.3f} kg/h")
    print(f"Heater energy: {E_heater:.3f} kWh_th/cycle")
    print(f"Cooling duty: {E_cooler:.3f} kWh_th/cycle")
    print(f"Thermal SEC: {SEC_thermal:.3f} kWh_th/kgCO2")


if __name__ == "__main__":
    main()
