"""Compare Figure 2 Python outputs with the exported MATLAB numerical reference."""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REF = ROOT / "validation_reference"
OUT = ROOT / "outputs" / "data"
CASES = ("1_cold_and_dry", "2_moderate", "3_hot_and_humid")


def validate(tolerance=1e-9):
    records = []
    for case in CASES:
        matlab = pd.read_csv(REF / f"matlab_fig2_{case}.csv")
        python = pd.read_csv(OUT / f"fig2_{case}.csv")
        assert len(matlab) == len(python) == 6185, case
        assert np.array_equal(
            matlab["time_s"].to_numpy(), python["time_s"].to_numpy()
        ), case
        assert np.array_equal(
            matlab["stage_id"].to_numpy(), python["stage_id"].to_numpy()
        ), case
        for var in python.columns:
            assert var in matlab.columns, f"MATLAB reference missing {var}"
            error = python[var].to_numpy(dtype=float) - matlab[var].to_numpy(
                dtype=float
            )
            rec = dict(
                case=case,
                variable=var,
                n=len(error),
                rmse=float(np.sqrt(np.mean(error * error))),
                mae=float(np.mean(np.abs(error))),
                max_abs=float(np.max(np.abs(error))),
            )
            records.append(rec)
            assert rec["max_abs"] <= tolerance, (case, var, rec["max_abs"])
    table = pd.DataFrame(records)
    table.to_csv(OUT / "matlab_python_parity_metrics.csv", index=False)

    ref = pd.read_csv(REF / "matlab_summary_metrics.csv")
    got = pd.read_csv(OUT / "summary_metrics.csv")
    assert len(ref) == len(got) == 3
    for col in [
        "qco2_start",
        "qco2_end",
        "qh2o_start",
        "qh2o_end",
        "Prod_tCO2_yr",
        "SEC_MWh_tCO2",
    ]:
        delta = np.abs(ref[col].to_numpy() - got[col].to_numpy())
        assert np.max(delta) <= tolerance, (col, np.max(delta))
    print("PASSED: 3 x 6185 time points; all 17 common columns within", tolerance)
    print("Maximum absolute difference across all columns:", table["max_abs"].max())
    print("MATLAB vs Python productivity / SEC: all 3 cases within", tolerance)
    print(
        table[
            table["variable"].isin(
                [
                    "aw_co2_effective",
                    "q_co2_total_mol_kg",
                    "q_h2o_total_mol_kg",
                    "T_bed_K",
                    "y_co2_g",
                ]
            )
        ][["case", "variable", "rmse"]].to_string(index=False)
    )
    return table


if __name__ == "__main__":
    validate()
