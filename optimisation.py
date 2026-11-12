"""Evaluate fifth-cycle heater energy and CO2 productivity for an operating point."""

import numpy as np

from cycle_driver import make_driver
from parameters import ambient, baseline_parameters
from simulator_thermal import forward_cycle


def evaluate(T_regenC, t_des_s, dt_s=5.0, T_ambC=20.0, RH_in=0.50):
    """Return fifth-cycle performance at the requested operating and ambient conditions.

    T_regenC and T_ambC are in Celsius, t_des_s and dt_s in seconds,
    and RH_in is a fraction. Results contain productivity [kg CO2/h],
    SEC [kWh thermal/kg CO2], E_heater [kWh] and m_co2 [kg/cycle].
    Five cycles are used consistently; convergence is not checked."""
    ctrl, design, phys, opt = baseline_parameters()

    ctrl["T_regenC"] = T_regenC
    ctrl["time_des_s"] = t_des_s
    opt["dt_s"] = dt_s

    amb = ambient(T_ambC, RH_in)
    drv = make_driver(amb, ctrl, opt)

    # Run repeated cycles
    x_init = None

    for _ in range(5):
        result = forward_cycle(drv, amb, ctrl, design, phys, opt, x_init=x_init)

        x_init = result["x_end"].copy()

    prof = result["prof"]

    n_co2 = result["cyc"]["nCO2_captured_eff"]
    m_co2 = n_co2 * phys["eq"]["MW_CO2"]

    cycle_time_h = (drv["t_s"][-1] - drv["t_s"][0]) / 3600

    E_heater = np.sum(prof["Q_heater_W"][:-1]) * opt["dt_s"] / 3.6e6

    productivity = m_co2 / cycle_time_h
    SEC = E_heater / m_co2 if m_co2 > 0 else float("inf")

    return {
        "productivity": productivity,
        "SEC": SEC,
        "E_heater": E_heater,
        "m_co2": m_co2,
    }


if __name__ == "__main__":
    result = evaluate(100, 20000)

    print(f"Productivity: {result['productivity']:.3f} kg/h")
    print(f"Thermal SEC: {result['SEC']:.3f} kWh_th/kgCO2")
