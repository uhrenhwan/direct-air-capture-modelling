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

    peak_idx = np.argmax(sorption_heat)

    print(
        f"\nPeak adsorption heat: "
        f"{sorption_heat[peak_idx]/1000:.2f} kW "
        f"at {hours[peak_idx]*3600:.1f} s"
    )

    print("\nInitial adsorption heat rates:")
    for i in range(min(10, len(hours))):
        print(
            f"t = {hours[i]*3600:7.1f} s | "
            f"CO2 = {co2_heat[i]/1000:8.2f} kW | "
            f"H2O = {water_heat[i]/1000:8.2f} kW"
        )

    # Plot initial and long-term adsorption heat behaviour
    from pathlib import Path

    time_s = time - time[0]

    heat_terms = [
        (co2_heat, "CO2 adsorption heat"),
        (water_heat, "Water adsorption heat"),
        (sorption_heat, "Total adsorption heat"),
        (flow_heat, "Gas-flow sensible heat"),
        (sorption_heat + flow_heat, "Adsorption + gas-flow heat"),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    initial = time_s <= 120
    remaining = time_s >= 120

    for values, label in heat_terms:
        ax1.plot(time_s[initial], values[initial] / 1000, label=label)
        ax2.plot(time_s[remaining] / 3600, values[remaining] / 1000, label=label)

    ax1.set_title("Initial adsorption transient (0–120 s)")
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Heat transfer rate [kW]")

    ax2.set_title("Subsequent adsorption behaviour")
    ax2.set_xlabel("Adsorption time [h]")
    ax2.set_ylabel("Heat transfer rate [kW]")

    for ax in (ax1, ax2):
        ax.axhline(0, color="black", linewidth=0.8)
        ax.grid(alpha=0.3)

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3)

    fig.suptitle("TVSA Adsorption: Dynamic Thermal Energy Balance")
    fig.tight_layout(rect=[0, 0.12, 1, 0.93])

    output_dir = Path("outputs/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_dir / "adsorption_thermal_balance.png", dpi=300, bbox_inches="tight"
    )

    plt.show()


if __name__ == "__main__":
    main()
