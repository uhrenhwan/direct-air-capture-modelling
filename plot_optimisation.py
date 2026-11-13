"""Compare heater SEC at the baseline and recorded candidate operating points."""

from pathlib import Path

import matplotlib.pyplot as plt

from optimisation import evaluate

# Evaluate baseline and optimised operating conditions


def main():
    baseline = evaluate(100.0, 20000.0, dt_s=1.0)
    optimised = evaluate(95.0, 23076.65, dt_s=1.0)

    sec_values = [baseline["SEC"], optimised["SEC"]]
    labels = ["Baseline", "Recorded candidate"]

    reduction = (baseline["SEC"] - optimised["SEC"]) / baseline["SEC"] * 100

    fig, ax = plt.subplots(figsize=(7, 5))

    bars = ax.bar(labels, sec_values, width=0.5)

    ax.set_ylabel("Heater SEC (kWh$_{th}$/kg CO$_2$)")
    ax.set_title("TVSA Process Optimisation")
    ax.set_ylim(0, max(sec_values) * 1.2)

    for bar, value in zip(bars, sec_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + 0.03,
            f"{value:.4f}",
            ha="center",
            fontsize=11,
        )

    ax.text(
        0.5,
        max(sec_values) * 1.12,
        f"SEC Reduction: {reduction:.2f}%",
        ha="center",
        fontsize=12,
        fontweight="bold",
    )

    fig.tight_layout()

    output_dir = Path(__file__).resolve().parent / "outputs" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        output_dir / "optimisation_comparison.png", dpi=300, bbox_inches="tight"
    )

    plt.show()


if __name__ == "__main__":
    main()
