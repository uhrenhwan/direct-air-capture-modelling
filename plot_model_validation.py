"""Compare baseline MATLAB and Python trajectories and report numerical agreement."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from validate_matlab_parity import validate

ROOT = Path(__file__).resolve().parent
REFERENCE = ROOT / "validation_reference"
DATA = ROOT / "outputs" / "data"
FIGURES = ROOT / "outputs" / "figures"


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)

    CASES = [
        ("1_cold_and_dry", "Cold and dry\n4°C, 20% RH"),
        ("2_moderate", "Moderate\n20°C, 50% RH"),
        ("3_hot_and_humid", "Hot and humid\n40°C, 90% RH"),
    ]

    VARIABLES = [
        ("q_co2_total_mol_kg", "CO₂ loading [mol/kg]"),
        ("T_bed_K", "Bed temperature [°C]"),
    ]

    fig, axes = plt.subplots(
        2, 3, figsize=(15, 7.5), sharex="col", constrained_layout=True
    )

    fig.suptitle(
        "TVSA Baseline Comparison: MATLAB vs Python", fontsize=16, fontweight="bold"
    )

    records = []

    for j, (case, title) in enumerate(CASES):
        matlab_path = REFERENCE / f"matlab_fig2_{case}.csv"
        python_path = DATA / f"fig2_{case}.csv"

        for path in (matlab_path, python_path):
            if not path.exists():
                raise FileNotFoundError(
                    f"Missing data: {path}\n"
                    "Run 'python run_baseline.py' to regenerate "
                    "Python results if necessary."
                )

        matlab = pd.read_csv(matlab_path)
        python = pd.read_csv(python_path)

        # Check the simulations use identical time grids
        if not np.array_equal(matlab["time_s"].to_numpy(), python["time_s"].to_numpy()):
            raise ValueError(f"Time grids do not match for {case}")

        if not np.array_equal(
            matlab["stage_id"].to_numpy(), python["stage_id"].to_numpy()
        ):
            raise ValueError(f"Stage sequences do not match for {case}")

        time_h = python["time_s"].to_numpy() / 3600

        for i, (variable, ylabel) in enumerate(VARIABLES):
            ax = axes[i, j]

            y_matlab = matlab[variable].to_numpy(dtype=float)
            y_python = python[variable].to_numpy(dtype=float)

            # Convert K to °C for presentation
            if variable == "T_bed_K":
                y_matlab = y_matlab - 273.15
                y_python = y_python - 273.15

            error = y_python - y_matlab

            rmse = np.sqrt(np.mean(error**2))
            max_error = np.max(np.abs(error))

            ax.plot(time_h, y_matlab, linewidth=2, label="MATLAB reference")

            ax.plot(
                time_h,
                y_python,
                linestyle="--",
                linewidth=1.5,
                marker="o",
                markevery=300,
                markersize=2.5,
                label="Python implementation",
            )

            ax.text(
                0.97,
                0.05,
                f"RMSE = {rmse:.2e}\n" f"Max error = {max_error:.2e}",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=9,
                bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "lightgrey"},
            )

            ax.grid(alpha=0.25)

            if j == 0:
                ax.set_ylabel(ylabel)

            if i == 0:
                ax.set_title(title, fontsize=12)

            if i == 1:
                ax.set_xlabel("Cycle time [h]")

            records.append(
                {
                    "case": case,
                    "variable": variable,
                    "RMSE": rmse,
                    "max_absolute_error": max_error,
                }
            )

    axes[0, 0].legend(fontsize=9)

    figure_path = FIGURES / "matlab_python_validation.png"
    csv_path = DATA / "matlab_python_validation_plot_metrics.csv"

    fig.savefig(figure_path, dpi=300, bbox_inches="tight")

    pd.DataFrame(records).to_csv(csv_path, index=False)

    print("\nPlot comparison results:")
    print(pd.DataFrame(records).to_string(index=False))

    print("\nRunning full MATLAB-Python numerical validation...")

    validate(tolerance=1e-9)

    print("\nValidation figure saved to:")
    print(figure_path)

    print("\nPlot metrics saved to:")
    print(csv_path)

    plt.show()


if __name__ == "__main__":
    main()
