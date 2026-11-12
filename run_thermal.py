"""Inspect one thermal cycle and plot heat transfer during adsorption."""

import matplotlib.pyplot as plt
import numpy as np

from cycle_driver import STAGE_NAMES, make_driver
from parameters import ambient, baseline_parameters
from simulator_thermal import forward_cycle


def main():
    ctrl, design, phys, opt = baseline_parameters()
    amb = ambient(20.0, 0.50)
    driver = make_driver(amb, ctrl, opt)
    result = forward_cycle(driver, amb, ctrl, design, phys, opt)
    profile = result["prof"]
    temperature = profile["TC_bed"]

    print(
        f"Bed temperature range: {temperature.min():.2f} to {temperature.max():.2f} C"
    )
    print(f"Final temperature: {temperature[-1]:.2f} C")
    print("All temperatures finite:", np.isfinite(temperature).all())
    for stage_id, name in enumerate(STAGE_NAMES, start=1):
        indices = np.flatnonzero(driver["stage_id"] == stage_id)
        if indices.size:
            print(
                f"{name}: {temperature[indices[0]]:.2f} to {temperature[indices[-1]]:.2f} C"
            )

    # Heat rates describe intervals, so exclude the terminal sample.
    adsorption = driver["stage_id"][:-1] == 6
    time = profile["time_s"][:-1][adsorption]
    hours = (time - time[0]) / 3600
    co2_heat = profile["Q_ads_CO2_W"][:-1][adsorption]
    water_heat = profile["Q_ads_H2O_W"][:-1][adsorption]
    sorption_heat = profile["Q_ads_heat_W"][:-1][adsorption]
    flow_heat = profile["Q_flow_W"][:-1][adsorption]

    fig, ax = plt.subplots(figsize=(9, 5))
    for values, label, style in [
        (co2_heat, "CO2 adsorption heat", "-"),
        (water_heat, "Water adsorption heat", "-"),
        (sorption_heat, "Total adsorption heat", "-"),
        (flow_heat, "Gas-flow sensible heat", "-"),
        (sorption_heat + flow_heat, "Net adsorption and gas-flow heat", "--"),
    ]:
        ax.plot(hours, values / 1000, style, label=label)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set(
        xlabel="Adsorption time [h]",
        ylabel="Heat transfer rate [kW]",
        title="Thermal energy balance during adsorption",
    )
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    error = np.max(np.abs(sorption_heat - co2_heat - water_heat))
    print(f"Maximum sorption heat decomposition error: {error:.6e} W")
    plt.show()


if __name__ == "__main__":
    main()
