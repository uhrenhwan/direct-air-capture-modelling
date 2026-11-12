"""Compare baseline and candidate operating points at 5 s and 1 s time steps."""

from optimisation import evaluate


def main():
    # Baseline and best SciPy solution
    T_opt = 94.84
    t_opt = 23066.14

    for dt in [5.0, 1.0]:
        baseline = evaluate(100, 20000, dt)
        optimal = evaluate(T_opt, t_opt, dt)

        # Recalculate constraint using the same timestep
        P_min = 0.9 * baseline["productivity"]

        reduction = (1 - optimal["SEC"] / baseline["SEC"]) * 100

        print(f"\n--- dt = {dt} s ---")

        print(f"Baseline SEC: {baseline['SEC']:.4f}")
        print(f"Optimised SEC: {optimal['SEC']:.4f}")

        print(f"SEC reduction: {reduction:.2f}%")

        print(f"Baseline productivity: {baseline['productivity']:.4f}")
        print(f"Optimised productivity: {optimal['productivity']:.4f}")

        print(f"Minimum productivity: {P_min:.4f}")
        print(f"Feasible: {optimal['productivity'] >= P_min}")


if __name__ == "__main__":
    main()
