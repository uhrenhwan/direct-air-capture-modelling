"""Figure 2 forward-memory TVSA run. Usage: python run_baseline.py"""

from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cycle_driver import extract_adsorption_driver, make_driver
from energy import energy_metrics
from parameters import ambient, baseline_parameters
from simulator import fig2_simulations

OUT = Path(__file__).resolve().parent / "outputs"


def run():
    ctrl, design, phys, opt = baseline_parameters()
    _, cases = fig2_simulations(
        ctrl, design, phys, opt, make_driver, extract_adsorption_driver, ambient
    )
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    (OUT / "data").mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(5, 3, figsize=(13, 9), sharex="col")
    ranges = [(0, 1.6), (0, 14), (0, 105), (0, 1.05), (0, 1.05)]
    rows = []
    for j, (name, amb, adapt, sim) in enumerate(cases):
        drv = sim["drv"]
        p = sim["prof"]
        t = drv["t_s"] / 3600
        E = energy_metrics(sim, amb, ctrl, design, phys, opt)
        traces = [
            p["qCO2_total"],
            p["qH2O_total"],
            p["TC_bed"],
            drv["RH_cmd"],
            drv["P_Pa"] / 1e5,
        ]
        for i, (a, y) in enumerate(zip(axes[:, j], traces)):
            a.plot(t, y, lw=1.2)
            if i == 2:
                a.plot(
                    t, drv["TK_cmd"] - 273.15, "--", c="0.5", lw=1, label="T command"
                )
            a.set_ylim(*ranges[i])
            a.grid(alpha=0.15)
            for edge in drv["edges_s"][1:-1] / 3600:
                a.axvline(edge, lw=0.5, c="0.75")
        axes[0, j].set_title(
            f'Figure 2({chr(97+j)}) — {name}\n{amb["T_ambC"]:.0f} C, {amb["RH_in"]*100:.0f}% RH'
        )
        axes[4, j].set_xlabel("Cycle time [h]")
        df = pd.DataFrame(
            {
                "time_s": drv["t_s"],
                "stage_id": drv["stage_id"],
                "q_co2_total_mol_kg": p["qCO2_total"],
                "q_h2o_total_mol_kg": p["qH2O_total"],
                "q_co2_fast_mol_kg": p["qCO2_fast"],
                "q_co2_slow_mol_kg": p["qCO2_slow"],
                "q_h2o_free_mol_kg": p["qH2O_free"],
                "q_h2o_bound_mol_kg": p["qH2O_bound"],
                "y_co2_g": p["yCO2g"],
                "T_bed_K": p["TK_bed"],
                "T_cmd_K": drv["TK_cmd"],
                "RH_cmd": drv["RH_cmd"],
                "P_Pa": drv["P_Pa"],
                "aw_gas": p["aw_gas"],
                "aw_co2_effective": p["aw_eff"],
                "Qin_m3_s": p["q_in"],
                "Qout_m3_s": p["q_out"],
            }
        )
        df.to_csv(
            OUT / "data" / f'fig2_{j+1}_{name.lower().replace(" ","_")}.csv',
            index=False,
        )
        rows.append(
            dict(
                case=name,
                T_C=amb["T_ambC"],
                RH=amb["RH_in"],
                qco2_start=p["qCO2_total"][0],
                qco2_end=p["qCO2_total"][-1],
                qh2o_start=p["qH2O_total"][0],
                qh2o_end=p["qH2O_total"][-1],
                Prod_tCO2_yr=E["Prod_tCO2_collector_yr"],
                SEC_MWh_tCO2=E["Elec_total_MWh_tCO2"],
            )
        )
    for i, lab in enumerate(
        [
            "CO2 loading [mol/kg]",
            "H2O loading [mol/kg]",
            "Bed temperature [C]",
            "RH command [-]",
            "Bed pressure [bar]",
        ]
    ):
        axes[i, 0].set_ylabel(lab)
    fig.suptitle(
        "Reimplemented single-cycle TVSA forward response after ambient-memory adaptation",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = OUT / "figures" / "figure2_python_reimplementation.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / "data" / "summary_metrics.csv", index=False)
    print(metrics.to_string(index=False))
    print("Wrote", path)
    return metrics


if __name__ == "__main__":
    run()
