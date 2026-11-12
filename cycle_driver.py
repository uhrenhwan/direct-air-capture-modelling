"""Prescribed pressure, humidity and temperature commands for the six TVSA stages."""

from math import exp

import numpy as np


def clamp(x, lo, hi):
    return min(max(x, lo), hi)


def smoothstep(u):
    u = clamp(u, 0.0, 1.0)
    return u * u * (3 - 2 * u)


def ambient_pressure_bar(alt):
    return 1.01325 * exp(-alt / 8400.0)


STAGE_NAMES = [
    "Vacuum",
    "Heating",
    "Desorption",
    "Cooling",
    "Pressurization",
    "Adsorption",
]


def make_driver(amb, ctrl, opt):
    """Build the six-stage command arrays; pressure is in Pa and temperature in K."""
    dt_s = opt["dt_s"]
    edges_s = np.cumsum(
        [
            0.0,
            ctrl["time_vac_s"],
            ctrl["time_heat_s"],
            ctrl["time_des_s"],
            ctrl["time_cool_s"],
            ctrl["time_press_s"],
            ctrl["time_ads_s"],
        ]
    )
    time_s = np.arange(0.0, edges_s[-1] + dt_s / 2, dt_s)
    stage_ids = np.searchsorted(edges_s[1:6], time_s, side="right").astype(np.uint8) + 1
    temperature_c = np.zeros(time_s.size)
    humidity = np.zeros(time_s.size)
    pressure_bar = np.zeros(time_s.size)
    ambient_temperature_c = amb["T_ambC"]
    inlet_rh = amb["RH_in"]
    ambient_pressure = ambient_pressure_bar(amb["altitude_m"])
    vacuum_pressure = ctrl["p_vac_bar"]
    regeneration_temperature_c = ctrl["T_regenC"]
    heat_end = ambient_temperature_c + 0.78 * (
        regeneration_temperature_c - ambient_temperature_c
    )
    for k, (time, stage_id) in enumerate(zip(time_s, stage_ids)):
        if stage_id == 1:
            tau = max(0.20 * (edges_s[1] - edges_s[0]), 1e-6)
            pressure_bar[k] = vacuum_pressure + (
                ambient_pressure - vacuum_pressure
            ) * exp(-time / tau)
        elif stage_id in (2, 3):
            pressure_bar[k] = vacuum_pressure
        elif stage_id == 4:
            pressure_bar[k] = max(vacuum_pressure - 0.015, 0.001)
        elif stage_id == 5:
            tau = max(0.35 * (edges_s[5] - edges_s[4]), 1e-6)
            pressurisation_time = max(time - edges_s[4], 0.0)
            start = max(vacuum_pressure - 0.015, 0.001)
            pressure_bar[k] = ambient_pressure - (ambient_pressure - start) * exp(
                -pressurisation_time / tau
            )
        else:
            pressure_bar[k] = ambient_pressure
        if stage_id == 1:
            temperature_c[k] = ambient_temperature_c
        elif stage_id == 2:
            u = clamp((time - edges_s[1]) / max(edges_s[2] - edges_s[1], 1e-12), 0, 1)
            if u <= 0.05:
                frac = 0.1 * smoothstep(u / 0.05)
            else:
                frac = 0.10 + 0.90 * ((u - 0.05) / 0.95) ** 1.05
            temperature_c[k] = (
                ambient_temperature_c + (heat_end - ambient_temperature_c) * frac
            )
        elif stage_id == 3:
            u = clamp((time - edges_s[2]) / max(edges_s[3] - edges_s[2], 1e-12), 0, 1)
            g = (1 - exp(-u / 0.55)) / (1 - exp(-1 / 0.55))
            temperature_c[k] = heat_end + (regeneration_temperature_c - heat_end) * g
        else:
            cooling_time = max(time - edges_s[3], 0.0)
            tau = max(0.10 * (edges_s[6] - edges_s[3]), 1e-6)
            temperature_c[k] = ambient_temperature_c + (
                regeneration_temperature_c - ambient_temperature_c
            ) * exp(-cooling_time / tau)
        if stage_id == 1:
            humidity[k] = inlet_rh
        elif stage_id == 2:
            u = clamp((time - edges_s[1]) / max(edges_s[2] - edges_s[1], 1e-12), 0, 1)
            humidity[k] = (1 - u) * inlet_rh + u * min(0.98, inlet_rh * 1.10)
        elif stage_id == 3:
            u = clamp((time - edges_s[2]) / max(edges_s[3] - edges_s[2], 1e-12), 0, 1)
            humidity[k] = 0.05 + (min(0.98, inlet_rh * 1.10) - 0.05) * exp(-5 * u)
        elif stage_id in (4, 5):
            humidity[k] = 0.05 * 0.85
        else:
            u = clamp((time - edges_s[5]) / max(edges_s[6] - edges_s[5], 1e-12), 0, 1)
            rec = 1 - exp(-((u / 0.24) ** 1.55))
            humidity[k] = 0.05 * 0.85 + (inlet_rh - 0.05 * 0.85) * rec
        humidity[k] = clamp(humidity[k], 0, 0.999)
    # Smooth pressure jumps at stage boundaries using the reference ramp shape.
    pressure = np.maximum(pressure_bar * 1e5, 200.0)
    pressure_options = opt["ptraj"]
    if pressure_options["driver_edge_smooth_enable"]:
        ramp = pressure_options["driver_edge_ramp_s"]
        n = len(time_s)
        for k in range(1, n):
            if stage_ids[k] != stage_ids[k - 1]:
                i0 = k - 1
                ramp_steps = max(1, round(ramp / max(dt_s, 1e-12)))
                j2 = min(n - 1, k + ramp_steps)
                start_pressure = pressure[i0]
                end_pressure = pressure[j2]
                m = j2 - i0
                if m >= 2:
                    u = np.arange(m + 1) / m
                    s = u * u * (3 - 2 * u)
                    pressure[i0 : j2 + 1] = (
                        start_pressure + (end_pressure - start_pressure) * s
                    )
    return dict(
        t_s=time_s,
        stage_id=stage_ids,
        edges_s=edges_s,
        stage_names=STAGE_NAMES,
        TK_cmd=temperature_c + 273.15,
        RH_cmd=humidity,
        P_Pa=pressure,
        TC_cmd=temperature_c,
        Pbar=pressure / 1e5,
        meta={"source": "forced_only"},
    )


def extract_adsorption_driver(full, opt):
    """Copy the adsorption segment and reset its time origin for preconditioning."""
    idx = np.where(full["stage_id"] == 6)[0]
    driver = {}
    for key, value in full.items():
        if (
            isinstance(value, np.ndarray)
            and value.ndim == 1
            and value.size == full["t_s"].size
        ):
            driver[key] = value[idx].copy()
        else:
            driver[key] = value
    sample_count = idx.size
    driver["t_s"] = np.arange(sample_count) * opt["dt_s"]
    driver["stage_id"] = np.full(sample_count, 6, dtype=np.uint8)
    driver["TK_cmd"] = full["TK_cmd"][idx].copy()
    driver["TC_cmd"] = driver["TK_cmd"] - 273.15
    driver["RH_cmd"] = full["RH_cmd"][idx].copy()
    driver["P_Pa"] = full["P_Pa"][idx].copy()
    driver["Pbar"] = driver["P_Pa"] / 1e5
    driver["edges_s"] = np.array([0.0, sample_count * opt["dt_s"]])
    driver["stage_names"] = ["Adsorption"]
    driver["meta"] = {"source": "adsorption_preconditioning_driver"}
    return driver
