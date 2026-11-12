"""Baseline operating conditions, sorbent properties and thermal equipment settings."""

from copy import deepcopy

ctrl = {
    "u_feed_m_s": 0.028,  # Feed superficial velocity (m/s)
    "p_vac_bar": 0.05,  # Vacuum pressure (bar)
    "T_regenC": 100,  # Regeneration temperature (C)
    "time_vac_s": 60,  # Vacuum duration (s)
    "time_heat_s": 2400,  # Heating duration (s)
    "time_des_s": 20000,  # Desorption duration (s)
    "time_cool_s": 400,  # Cooling duration (s)
    "time_press_s": 60,  # Pressurization duration (s)
    "time_ads_s": 8000,  # Adsorption duration (s)
}

# design parameters - geometry and bed properties (design)

design = {
    "width_cell_m": 1.43,  # Width of a single cell (m)
    "height_cell_m": 0.10,  # Height of a single cell (m)
    "N_cell": 13 * 88,  # Number of cells in the module
    "L_flow_fixed_m": 0.0172,  # Flow length (m)
    "rb_fixed": 1000 / 2.806,  # Bulk sorbent density (kg/m3)
    "eps": 0.56,  # Bed void fraction (-)
    "dp_m": 7.50e-4,  # Particle diameter (m)
    "A_mult": 2.1,  # Area multiplier (-)
    "Ng_ref": 35,  # Effective gas-volume multiplier (-)
    "dP_misc_Pa": 120,  # Miscellaneous pressure drop (Pa)
    "N_series_pass": 1,  # Number of series passes (-)
}


def finalize_design(design):
    # Area of a single cell [m²]
    A_cell = design["width_cell_m"] * design["height_cell_m"]
    # Total flow and contact areas [m²]
    design["A_flow_m2"] = max(design["N_cell"] * A_cell, 1e-12)
    design["A_contact_m2"] = max(design["A_flow_m2"] * design["A_mult"], 1e-12)
    # Bed dimensions and volume
    design["L_flow_m"] = design["L_flow_fixed_m"]
    design["Vbed_m3"] = design["A_flow_m2"] * design["L_flow_m"]
    # Sorbent mass [kg]
    design["msorb_kg"] = design["rb_fixed"] * design["Vbed_m3"]
    # Effective contact length [m]
    design["L_contact_m"] = design["Vbed_m3"] / max(design["A_contact_m2"], 1e-12)
    # Effective gas volume [m³]
    design["Vg_m3"] = design["Ng_ref"] * design["eps"] * design["Vbed_m3"]
    design["Vg_ref_m3"] = design["Vg_m3"]

    return design


# Calculate derived design parameters
design = finalize_design(design)

# Ambient cases

amb_base = {
    "altitude_m": 0,  # Altitude [m]
    "T_ambC": 20.0,  # Ambient temperature [°C]
    "RH_in": 0.50,  # Relative humidity [-]
    "yCO2_air": 400e-6,  # Atmospheric CO2 mole fraction [-]
}

ambient_cases = [
    {"T_ambC": 4.0, "RH_in": 0.20},  # Cold and dry
    {"T_ambC": 20.0, "RH_in": 0.50},  # Moderate
    {"T_ambC": 40.0, "RH_in": 0.90},  # Hot and humid
]

# Physical parameters

phys = {
    "eq": {
        "R": 8.314462618,
        # Molecular weights [kg/mol]
        "MW_CO2": 44e-3,
        "MW_H2O": 18e-3,
        "MW_air": 28.97e-3,
        # Heat capacities
        "Cp_CO2_g": 37.0,  # [J/mol/K]
        "Cp_H2O_g": 34.0,  # [J/mol/K]
        "Cp_air_g": 29.0,  # [J/mol/K]
        "Cp_s": 1880,  # [J/kg/K]
        # Adsorption enthalpy magnitudes [J/mol]
        "dHads_CO2": 70e3,
        "dHads_H2O": 46e3,
    }
}

# Bed Properties

phys["bed"] = {"rb": design["msorb_kg"] / max(design["Vbed_m3"], 1e-12)}

# Adsorption Isotherm Parameters

phys["iso"] = {
    # CO2 mechanistic isotherm
    "mech": {
        "T0": 298.15,
        "qinf0": 4.86,
        "chi": 0.0,
        "b0": 2.85e-21,
        "DH_dry": 117798,
        "DH_wet": 130155,
        "tau0": 0.209,
        "alpha": 0.523,
        "phi_max": 1.000,
        "fblk_max": 0.433,
        "k_blk": 0.795,
        "phi_dry": 1.000,
        "Aexp": 1.535,
        "n_blk": 1.425,
        "S_cap": 1.0,
        "dry_scale": 1.0,
    },
    # H2O GAB isotherm
    "gab": {"qm": 3.63, "C": 47110, "D": 0.023744, "F": 57706, "G": -47.814},
}

# Adsorption Kinetic Parameters

phys["kin"] = {
    # Reference kinetic constants
    "kCO2_ref": 3.0e-3,
    "kH2O_ref": 8.6e-3,
    # Activation energies [J/mol]
    "Ea_CO2": 35e3,
    "Ea_H2O": 20e3,
    # Reference conditions
    "Tref": 293.15,
    "p_ref_Pa": 40.0,
    # Pressure-dependence parameters
    "beta": 1.0,
    "n_exp": 1.0,
    # Stage multipliers
    # Order: vacuum, heating, desorption,
    #        cooling, pressurization, adsorption
    "mCO2": [0.2, 0.002, 0.009, 0.001, 0.001, 1.0],
    "mH2O": [0.20, 0.18, 0.10, 0.20, 0.25, 1.00],
}

# Figure 2 Parameter Adjustments

phys["iso"]["mech"]["S_cap"] = 1.1

phys["kin"]["kCO2_ref"] *= design["A_mult"] ** 0.5

phys["kin"]["kCO2_ads_mult"] = 1.00

# Additional energy-related parameters
phys["eq"].update(
    {
        "r_steam_CO2": 0.9,
        "TsatC": 100.0,
        "Cp_liq": 4180,
        "hvap": 2257000,
        "eta2nd": 0.50,
        "dTmin_K": 10.0,
        "eta_fan": 0.60,
        "eta_comp": 0.60,
        "P_CO2_out_bar": 150,
        "T_comp_K": 298.15,
    }
)


# Simulation Options
# Reference: TVSA_Fig2.m - Figure 2 baseline

opt = {
    # Numerical settings
    "dt_s": 5.0,
    "max_iters": 1,
    "tol_rel": 5e-4,
    "omega": 1.00,
    # Model configuration
    "mode": "S2",
    "driver_baseline": {"mode": "forced_only"},
    "flow_mode": "Acontact",
    "legacy": {"keepQ_with_Acontact": 1},
    # Water activity
    "aw": {
        "max": 0.999,
        "chain_mode": "bed_consistent",
        "cmd_based_cap_policy": "cap_to_1",
        "mode_for_co2_effects": "invert_qH2O",
    },
    # CO2 partial pressure
    "pco2": {
        "ads_mode": "inlet",
        "regen_mode": "inventory",
        "blend_alpha_inventory": 0.20,
    },
    # Gas-phase CO2 mole fraction
    "yco2": {"y_floor": 1e-6, "y_max": 0.60, "dy_cap": {"enable": 1, "dy_max": 0.010}},
    # Pressure trajectory and flow smoothing
    "ptraj": {
        "enable_smoothing": 1,
        "tau_flow_ramp_s": 25,
        "Qmax_mult_Vg_dt": 20,
        "log_ndry_residual": 1,
        "driver_edge_smooth_enable": 1,
        "driver_edge_ramp_s": 60,
        "enable_ndry_filter": 1,
        "tau_ndry_s": {
            "vac": 120,
            "heat": 180,
            "des": 180,
            "cool": 120,
            "press": 20,
            "ads": 0,
        },
        "ndry_reset_on_stage_change": 1,
        "reset_flow_filters_on_stage_change": 1,
        "enable_ndot_rate_limit": 1,
    },
    # Numerical stabilisation
    "stab": {
        "enable": 1,
        "pdrv_lpf_enable": 1,
        "pdrv_tau_s": 20,
        "pdrv_stagejump_reset": 1,
        "aw_lpf_enable": 1,
        "aw_tau_s": 20,
        "aw_stagejump_reset": 1,
    },
    # Thermal response time constants [s]
    "tauT_s": {"vac": 0, "heat": 900, "des": 2000, "cool": 600, "press": 300, "ads": 0},
    # CO2 two-pool dynamics
    "co2_2pool": {
        "enable": 1,
        "fslow": {
            "f0_ads": 0.18,
            "f0_regen": 0.18,
            "k_aw_ads": 0.10,
            "k_aw_reg": 0.00,
            "f_min": 0.05,
            "f_max": 0.45,
        },
        "kslow": {"ratio": 0.02, "kmin": 2e-6, "kmax": 1e-3},
        "resid": {
            "enable": 1,
            "mode": "ads_scaled",
            "apply_in_stage_ids": [2, 3],
            "f_res": 0.08,
            "q_abs_min": 0.05,
            "k_dry": 0.0,
            "k_T": 0.05,
            "k_wet": 0.35,
            "n_wet": 1.0,
            "Tref_K": 373.15,
            "f_min": 0.02,
            "f_max": 0.25,
            "tol_frac": 1e-4,
        },
    },
    # H2O two-pool dynamics
    "h2o_2pool": {
        "enable": 1,
        "fbound": {
            "f0_ads": 0.15,
            "f0_regen": 0.22,
            "k_aw_ads": 0.10,
            "k_aw_reg": 0.02,
            "f_min": 0.05,
            "f_max": 0.60,
        },
        "kslow": {"ratio": 0.03, "kmin": 5e-7, "kmax": 5e-3},
        "resid": {
            "enable": 1,
            "apply_in_stage_ids": [1, 2, 3, 4, 5],
            "f0": 0.10,
            "f_max": 0.55,
            "k_T": 0.10,
            "k_time": 0.10,
            "Tref_K": 373.15,
            "tref_s": 20000,
            "q_abs_min": 0.05,
            "tol_frac": 1e-4,
        },
    },
    # Optional kinetic corrections
    "ads_kinshape": {"enable": 0, "kboost0": 1.0, "tau_k_s": 400, "fac_max": 1.0},
    "co2_awfloor": {"enable": 1, "apply_below": 0.10, "value": 0.05},
    "co2_dry_scale": {
        "enable": 0,
        "factor": 1.00,
        "aw0": 0.15,
        "p": 1.0,
        "apply_stages": [6],
    },
    "h2o_T_ads_scale": {"enable": 0, "TrefC": 20, "dT_spanC": 30, "k_gain": 0.45},
    # Other simulation settings
    "early_stop": {"enable": 0},
    "pub_trace": {"enable": 0},
}


# Exact Figure 2 options absent from the student's initial data-entry file
# are supplied by MATLAB ensure_opt_defaults() / getfield_safe().
def baseline_parameters():
    """Isolated copies for reproducible runs; edit controls without global side effects."""
    c, d, p, o = deepcopy(ctrl), deepcopy(design), deepcopy(phys), deepcopy(opt)
    # Active behaviour in source. Default values referenced at MATLAB
    # run_one_cycle_forward_base_v11 lines 1415-1426.
    o.setdefault("co2_sat", {"alpha": 0.80, "floor": 0.05})
    o.setdefault("ads_supply", {"enable": False, "keep_y_floor_inventory": True})
    o.setdefault("ads_shape", {"enable": False})
    o.setdefault("energy", {"vac_aw_ref": 0.25, "vac_aw_gain": 0.60})
    return c, d, p, o


AMBIENT_CASES = [
    ("Cold and dry", 4.0, 0.20),
    ("Moderate", 20.0, 0.50),
    ("Hot and humid", 40.0, 0.90),
]


def ambient(temp_c=20.0, rh=0.50):
    return dict(altitude_m=0.0, T_ambC=float(temp_c), RH_in=float(rh), yCO2_air=400e-6)


if __name__ == "__main__":
    print("Bed density:", phys["bed"]["rb"])
    print("CO2 capacity parameter:", phys["iso"]["mech"]["qinf0"])
    print("CO2 reference kinetic constant:", phys["kin"]["kCO2_ref"])
    print("Time step:", opt["dt_s"])
    print("Simulation mode:", opt["mode"])
    print("Heating thermal time constant:", opt["tauT_s"]["heat"])


equipment = {
    "T_h_C": 140.0,
    "T_c_C": 15.0,
    "UA_h_W_K": 862.97,
    "UA_c_W_K": 10086.21,
    "Q_h_max_W": 115e3,
    "Q_c_max_W": 950e3,
    "hysteresis_h_K": 1.0,
    "hysteresis_c_K": 1.0,
}
