"""Run one-variable sensitivity sweeps and save heater SEC curves and model results."""

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

from optimisation import evaluate

DT = 1.0  # Simulation timestep [s]

# Baseline operating and ambient conditions
BASELINE = {
    "T_regenC": 100.0,  # Regeneration temperature [°C]
    "t_des_s": 20000.0,  # Desorption duration [s]
    "T_ambC": 20.0,  # Ambient temperature [°C]
    "RH_in": 0.50,  # Ambient relative humidity [-]
}

# Sensitivity ranges
sweeps = [
    {
        "title": "SEC vs Regeneration Temperature",
        "field": "T_regenC",
        "values": [90, 95, 100, 105, 110],
        "xlabel": "Regeneration temperature [°C]",
        "scale": 1.0,
    },
    {
        "title": "SEC vs Desorption Duration",
        "field": "t_des_s",
        "values": [16000, 18000, 20000, 22000, 24000],
        "xlabel": "Desorption duration [h]",
        "scale": 1 / 3600,
    },
    {
        "title": "SEC vs Ambient Temperature",
        "field": "T_ambC",
        "values": [4, 10, 20, 30, 40],
        "xlabel": "Ambient temperature [°C]",
        "scale": 1.0,
    },
    {
        "title": "SEC vs Ambient Relative Humidity",
        "field": "RH_in",
        "values": [0.20, 0.35, 0.50, 0.70, 0.90],
        "xlabel": "Ambient relative humidity [%]",
        "scale": 100.0,
    },
]


def run_model(settings):
    """Evaluate the existing DAC thermal model."""

    output = evaluate(
        settings["T_regenC"],
        settings["t_des_s"],
        dt_s=DT,
        T_ambC=settings["T_ambC"],
        RH_in=settings["RH_in"],
    )

    if not (math.isfinite(output["SEC"]) and math.isfinite(output["productivity"])):
        raise ValueError(f"Non-finite model result at {settings}")

    return output


def main():
    print("Running baseline simulation...")
    baseline = run_model(BASELINE)

    print(f"Baseline SEC: {baseline['SEC']:.4f} kWh_th/kg CO2")
    print(f"Baseline productivity: {baseline['productivity']:.4f} kg/h")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)

    fig.suptitle("TVSA Sensitivity Analysis", fontsize=16, fontweight="bold")

    all_results = []

    for ax, sweep in zip(axes.flat, sweeps):
        field = sweep["field"]
        scale = sweep["scale"]

        x_values = []
        sec_values = []

        print(f"\n--- {sweep['title']} ---")

        for value in sweep["values"]:
            # Start from baseline and change ONE variable only
            settings = BASELINE.copy()
            settings[field] = float(value)

            output = run_model(settings)

            SEC = output["SEC"]
            productivity = output["productivity"]

            x_values.append(value * scale)
            sec_values.append(SEC)

            all_results.append(
                {
                    "analysis": sweep["title"],
                    "T_regen_C": settings["T_regenC"],
                    "t_des_s": settings["t_des_s"],
                    "T_amb_C": settings["T_ambC"],
                    "RH_pct": settings["RH_in"] * 100,
                    "SEC_kWh_th_per_kgCO2": SEC,
                    "productivity_kg_h": productivity,
                }
            )

            print(
                f"{field} = {value:.2f} | "
                f"SEC = {SEC:.4f} | "
                f"Productivity = {productivity:.4f}"
            )

        ax.plot(
            x_values,
            sec_values,
            marker="o",
            markersize=7,
            linewidth=2,
            label="Model prediction",
        )

        # Highlight baseline operating condition
        baseline_x = BASELINE[field] * scale

        ax.scatter(
            baseline_x,
            baseline["SEC"],
            marker="*",
            s=200,
            color="black",
            zorder=5,
            label="Baseline",
        )

        ax.set_title(sweep["title"], fontsize=12, pad=12)

        ax.set_xlabel(sweep["xlabel"], fontsize=11)

        ax.set_ylabel("Heater SEC [kWh$_{th}$/kg CO$_2$]", fontsize=11)

        ax.grid(alpha=0.25)
        ax.legend(fontsize=9)

        ax.tick_params(axis="both", labelsize=10)

    root = Path(__file__).resolve().parent

    figure_dir = root / "outputs" / "figures"
    data_dir = root / "outputs" / "data"

    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    figure_path = figure_dir / "tvsa_sensitivity_analysis.png"

    fig.savefig(figure_path, dpi=300, bbox_inches="tight")

    csv_path = data_dir / "tvsa_sensitivity_analysis.csv"

    with open(csv_path, "w", newline="", encoding="utf-8") as f:

        writer = csv.DictWriter(f, fieldnames=list(all_results[0].keys()))

        writer.writeheader()
        writer.writerows(all_results)

    print("\n===================================")
    print("SENSITIVITY ANALYSIS COMPLETE")
    print("===================================")
    print(f"Figure saved: {figure_path}")
    print(f"Data saved: {csv_path}")

    plt.show()


if __name__ == "__main__":
    main()
