"""Multi-start COBYLA search with a minimum CO2 productivity constraint."""

from functools import lru_cache

from scipy.optimize import minimize

from optimisation import evaluate


def main():
    # Baseline performance
    baseline = evaluate(100, 20000, dt_s=1.0)
    P_min = 0.90 * baseline["productivity"]

    # Convert normalised variables to physical units
    @lru_cache(maxsize=256)
    def evaluate_cached(u_T, u_t):
        T_regen = 80 + 40 * u_T
        t_des = 14000 + 12000 * u_t

        return evaluate(T_regen, t_des, dt_s=1.0)

    def get_metrics(x):
        return evaluate_cached(round(float(x[0]), 8), round(float(x[1]), 8))

    # Objective: minimise Thermal SEC
    def objective(x):
        return get_metrics(x)["SEC"]

    # Constraint: Productivity >= 90% of baseline
    def productivity_constraint(x):
        return get_metrics(x)["productivity"] - P_min

    # Multi-start optimisation
    starting_points = [[0.1, 0.1], [0.1, 0.9], [0.5, 0.5], [0.9, 0.1], [0.9, 0.9]]

    results = []

    for i, x0 in enumerate(starting_points, start=1):
        res = minimize(
            fun=objective,
            x0=x0,
            method="COBYLA",
            bounds=[(0, 1), (0, 1)],
            constraints=[{"type": "ineq", "fun": productivity_constraint}],
            options={"rhobeg": 0.15, "tol": 1e-3, "maxiter": 80, "catol": 1e-4},
        )

        metrics = get_metrics(res.x)

        T = 80 + 40 * res.x[0]
        t = 14000 + 12000 * res.x[1]

        feasible = metrics["productivity"] >= P_min - 1e-4

        print(
            f"Run {i}: "
            f"T = {T:.2f} °C | "
            f"t = {t:.0f} s | "
            f"SEC = {metrics['SEC']:.4f} | "
            f"P = {metrics['productivity']:.3f} | "
            f"Success = {res.success} | "
            f"Feasible = {feasible}"
        )

        if res.success and feasible:
            results.append(res)

    # Select best successful feasible solution
    if not results:
        raise RuntimeError("No successful feasible solution found.")

    result = min(results, key=lambda r: get_metrics(r.x)["SEC"])

    # Extract optimised operating conditions
    T_opt = 80 + 40 * result.x[0]
    t_opt = 14000 + 12000 * result.x[1]

    metrics = get_metrics(result.x)

    print("\n--- SCIPY OPTIMISATION RESULTS ---")

    print(f"Solver success: {result.success}")
    print(f"Solver message: {result.message}")

    print(f"Optimal T_regen: {T_opt:.2f} °C")
    print(f"Optimal t_des: {t_opt:.2f} s")

    print(f"Productivity: {metrics['productivity']:.3f} kg/h")
    print(f"Minimum productivity: {P_min:.3f} kg/h")

    print(f"Thermal SEC: {metrics['SEC']:.3f} kWh_th/kgCO2")

    reduction = (baseline["SEC"] - metrics["SEC"]) / baseline["SEC"] * 100

    print(f"SEC reduction: {reduction:.2f} %")

    feasible = metrics["productivity"] >= P_min - 1e-4

    print(f"Productivity constraint satisfied: {feasible}")


if __name__ == "__main__":
    main()
