"""Plot one thermal cycle, its reference temperature command and equipment duties."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cycle_driver import STAGE_NAMES, make_driver
from parameters import ambient, baseline_parameters
from simulator_thermal import forward_cycle

# Run baseline thermal model


def main():
    ctrl, design, phys, opt = baseline_parameters()
    amb = ambient(20.0, 0.50)

    driver = make_driver(amb, ctrl, opt)
    result = forward_cycle(driver, amb, ctrl, design, phys, opt)

    profile = result["prof"]

    time_h = profile["time_s"] / 3600

    T_bed = profile["TC_bed"]
    T_command = driver["TC_cmd"]

    Q_heater = profile["Q_heater_W"] / 1000
    Q_cooler = profile["Q_cooler_W"] / 1000

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 8), sharex=True, gridspec_kw={"height_ratios": [1.3, 1]}
    )

    ax1.plot(time_h, T_bed, label="Simulated bed temperature", linewidth=2)

    ax1.plot(
        time_h,
        T_command,
        "--",
        label="Prescribed temperature trajectory",
        linewidth=1.7,
    )

    ax1.set_ylabel("Temperature [°C]")
    ax1.set_title("Dynamic Bed Temperature Response")
    ax1.legend()

    ax2.plot(time_h, Q_heater, label="Heater duty", linewidth=1.5)

    ax2.plot(time_h, Q_cooler, label="Cooler duty", linewidth=1.5)

    ax2.set_xlabel("Cycle time [h]")
    ax2.set_ylabel("Thermal duty [kW]")
    ax2.set_title("Heating and Cooling Equipment Duty")
    ax2.legend()

    # Indicate the six TVSA cycle stages
    edges = driver["edges_s"] / 3600

    for ax in (ax1, ax2):
        for i, name in enumerate(STAGE_NAMES):
            start, end = edges[i], edges[i + 1]

            if i % 2 == 0:
                ax.axvspan(start, end, color="grey", alpha=0.08)

            if i > 0:
                ax.axvline(start, color="grey", linestyle=":", linewidth=0.9)

        ax.grid(alpha=0.25)

    fig.suptitle("TVSA Cycle: Dynamic Thermal Model", fontsize=16, y=0.98)

    fig.text(
        0.5,
        0.92,
        "Stage sequence: Vacuum → Heating → Desorption → "
        "Cooling → Pressurisation → Adsorption",
        ha="center",
        fontsize=10,
    )

    fig.tight_layout(rect=[0, 0, 1, 0.84])

    output_dir = Path(__file__).resolve().parent / "outputs" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(output_dir / "tvsa_cycle_temperature.png", dpi=300, bbox_inches="tight")

    print("Figure saved successfully.")
    print(f"Maximum bed temperature: {np.max(T_bed):.2f} °C")
    print(f"Minimum bed temperature: {np.min(T_bed):.2f} °C")

    plt.show()


if __name__ == "__main__":
    main()
