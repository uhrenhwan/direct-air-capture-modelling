"""Compare thermal performance across a 3 by 3 ambient temperature/humidity grid."""

import matplotlib.pyplot as plt
import pandas as pd

from optimisation import evaluate


def main():
    # Fixed optimised operating conditions
    T_regen = 95.0  # deg C
    t_des = 23077  # s
    dt = 1.0  # s

    # Ambient conditions
    temperatures = [15.0, 20.0, 25.0]
    humidities = [0.30, 0.50, 0.70]

    # Reference: moderate ambient baseline
    baseline = evaluate(100.0, 20000, dt_s=dt)
    P_min = 0.90 * baseline["productivity"]

    results = []

    for T_amb in temperatures:
        for RH in humidities:
            metrics = evaluate(T_regen, t_des, dt_s=dt, T_ambC=T_amb, RH_in=RH)

            results.append(
                {
                    "T_amb_C": T_amb,
                    "RH_pct": RH * 100,
                    "Productivity_kg_h": metrics["productivity"],
                    "SEC_kWh_th_kg": metrics["SEC"],
                    "Heater_kWh": metrics["E_heater"],
                    "Constraint_met": metrics["productivity"] >= P_min,
                }
            )

    df = pd.DataFrame(results)

    print("\n--- AMBIENT SENSITIVITY RESULTS ---")
    print(df.to_string(index=False, float_format="%.4f"))

    # Save results
    df.to_csv("ambient_sensitivity_results.csv", index=False)

    metrics = [
        ("SEC_kWh_th_kg", "Thermal SEC [kWh_th/kgCO2]"),
        ("Productivity_kg_h", "CO2 Productivity [kg/h]"),
    ]

    for column, title in metrics:
        matrix = df.pivot(index="T_amb_C", columns="RH_pct", values=column)

        fig, ax = plt.subplots(figsize=(7, 5))

        im = ax.imshow(matrix.values, origin="lower", aspect="auto")

        # Axis labels
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels([f"{x:.0f}%" for x in matrix.columns])

        ax.set_yticks(range(len(matrix.index)))
        ax.set_yticklabels([f"{y:.0f}" for y in matrix.index])

        # Display numbers inside cells
        midpoint = (matrix.values.min() + matrix.values.max()) / 2

        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix.iloc[i, j]

                ax.text(
                    j,
                    i,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    color="white" if value < midpoint else "black",
                )

        ax.set_xlabel("Relative Humidity")
        ax.set_ylabel("Ambient Temperature [°C]")
        ax.set_title(title)

        fig.colorbar(im, ax=ax)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()
