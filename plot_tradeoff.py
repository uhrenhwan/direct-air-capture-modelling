"""Plot heater SEC against CO2 productivity across a small operating grid."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt

from optimisation import evaluate

DT = 1.0  # Simulation timestep [s]


def main():
    temperatures = [95.0, 100.0, 105.0]  # °C
    desorption_times = [18000.0, 20000.0, 23076.65]  # s

    baseline_T = 100.0
    baseline_t = 20000.0

    opt_T = 95.0
    opt_t = 23076.65

    baseline = evaluate(baseline_T, baseline_t, dt_s=DT)

    P_min = 0.90 * baseline["productivity"]

    results = []

    for T in temperatures:
        for t in desorption_times:
            print(f"Evaluating T = {T:.1f} °C, " f"t_des = {t:.2f} s", flush=True)

            output = evaluate(T, t, dt_s=DT)

            P = output["productivity"]
            SEC = output["SEC"]

            results.append(
                {
                    "T_regen_C": T,
                    "t_des_s": t,
                    "productivity": P,
                    "SEC": SEC,
                    "feasible": P >= P_min,
                }
            )

            print(
                f"  SEC = {SEC:.4f} kWh_th/kg CO2 | "
                f"P = {P:.4f} kg/h | "
                f"Feasible = {P >= P_min}",
                flush=True,
            )

    fig, ax = plt.subplots(figsize=(9, 6))

    for T in temperatures:
        subset = [row for row in results if row["T_regen_C"] == T]

        P = [row["productivity"] for row in subset]
        SEC = [row["SEC"] for row in subset]

        ax.plot(
            P,
            SEC,
            marker="o",
            linewidth=1.8,
            markersize=7,
            label=f"T_regen = {T:.0f} °C",
        )

    # Minimum productivity constraint
    ax.axvline(
        P_min,
        color="grey",
        linestyle="--",
        linewidth=1.5,
        label="90% productivity constraint",
    )

    # Highlight baseline
    ax.scatter(
        baseline["productivity"],
        baseline["SEC"],
        marker="*",
        s=250,
        color="black",
        zorder=5,
        label="Baseline",
    )

    # Highlight optimised operating point
    optimised = evaluate(opt_T, opt_t, dt_s=DT)

    ax.scatter(
        optimised["productivity"],
        optimised["SEC"],
        marker="D",
        s=90,
        color="red",
        zorder=6,
        label="Recorded candidate",
    )

    ax.annotate(
        "Baseline",
        (baseline["productivity"], baseline["SEC"]),
        xytext=(10, 10),
        textcoords="offset points",
    )

    ax.annotate(
        "Recorded candidate",
        (optimised["productivity"], optimised["SEC"]),
        xytext=(15, 15),
        textcoords="offset points",
        ha="left",
        va="bottom",
    )

    ax.set_xlabel("CO₂ productivity [kg/h]")
    ax.set_ylabel("Heater SEC [kWh$_{th}$/kg CO₂]")

    ax.set_title("TVSA Optimisation: Energy–Productivity Trade-off")

    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)

    fig.tight_layout()

    root = Path(__file__).resolve().parent

    figure_dir = root / "outputs" / "figures"
    data_dir = root / "outputs" / "data"

    figure_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        figure_dir / "sec_productivity_tradeoff.png", dpi=300, bbox_inches="tight"
    )

    with open(data_dir / "sec_productivity_tradeoff.csv", "w", newline="") as f:

        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))

        writer.writeheader()
        writer.writerows(results)

    print("\nFigure and CSV saved successfully.")

    plt.show()


if __name__ == "__main__":
    main()
