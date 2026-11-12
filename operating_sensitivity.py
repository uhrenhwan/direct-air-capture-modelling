"""Explore regeneration temperature and desorption time near the candidate optimum."""

import matplotlib.pyplot as plt
import pandas as pd

from optimisation import evaluate


def main():
    # Operating conditions near the SciPy optimum
    temperatures = [93.0, 95.0, 97.0]
    desorption_times = [22000, 23076.65, 24000]

    dt = 1.0

    # Baseline productivity constraint
    baseline = evaluate(100, 20000, dt_s=dt)
    P_min = 0.9 * baseline["productivity"]

    results = []

    # Evaluate 3 x 3 operating conditions
    for T_regen in temperatures:
        for t_des in desorption_times:
            metrics = evaluate(T_regen, t_des, dt_s=dt, T_ambC=20.0, RH_in=0.50)

            results.append(
                {
                    "T_regen_C": T_regen,
                    "t_des_s": t_des,
                    "SEC": metrics["SEC"],
                    "Productivity": metrics["productivity"],
                    "Feasible": metrics["productivity"] >= P_min,
                }
            )

    df = pd.DataFrame(results)

    print("\n--- OPERATING SENSITIVITY RESULTS ---")
    print(df.to_string(index=False, float_format="%.4f"))

    df.to_csv("operating_sensitivity_results.csv", index=False)

    # Plot sensitivity heatmaps
    for metric, title in [
        ("SEC", "Thermal SEC [kWh_th/kgCO2]"),
        ("Productivity", "CO2 Productivity [kg/h]"),
    ]:

        matrix = df.pivot(index="T_regen_C", columns="t_des_s", values=metric)

        feasible = df.pivot(index="T_regen_C", columns="t_des_s", values="Feasible")

        fig, ax = plt.subplots(figsize=(8, 5))

        im = ax.imshow(matrix.values, origin="lower", aspect="auto")

        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([f"{x:.0f}" for x in matrix.columns])

        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([f"{y:.0f}" for y in matrix.index])

        # Display results and feasibility
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix.iloc[i, j]
                is_feasible = feasible.iloc[i, j]

                label = f"{value:.3f}"

                if not is_feasible:
                    label += "\n(Infeasible)"

                ax.text(
                    j,
                    i,
                    label,
                    ha="center",
                    va="center",
                    color=(
                        "white"
                        if value < (matrix.values.min() + matrix.values.max()) / 2
                        else "black"
                    ),
                )

        ax.set_xlabel("Desorption Time [s]")
        ax.set_ylabel("Regeneration Temperature [°C]")
        ax.set_title(title)

        fig.colorbar(im, ax=ax)
        plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()
