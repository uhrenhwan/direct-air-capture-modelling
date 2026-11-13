# Direct Air Capture | Reduced-Order TVSA Modelling and Thermal Optimisation

This project investigates the performance and energy requirements of a temperature-vacuum swing adsorption (TVSA) process for direct air capture (DAC) using Python-based process simulation.

The work began with reproducing a published reduced-order TVSA model in Python and comparing its outputs against MATLAB reference data. The model was subsequently extended to incorporate transient bed energy balances, adsorption heat effects, and heating/cooling equipment constraints. Constrained numerical optimisation was then performed to investigate the trade-off between regeneration energy consumption and CO₂ capture productivity.

**Key result:** A SciPy COBYLA optimisation identified operating conditions that reduced simulated heater-specific energy consumption by **2.31% (1.9951 → 1.9490 kWh_th/kg CO₂)** while maintaining at least 90% of baseline CO₂ productivity at a 1-second simulation timestep. This improvement comes with approximately a 10% reduction in productivity and is subject to the model's thermal assumptions.

## Engineering work

- Reproduced the six-stage adsorption cycle and compared the Python baseline with independent MATLAB exports.
- Replaced the prescribed temperature response with a lumped bed energy balance, including CO2/water sorption heat and gas-flow sensible heat.
- Added capacity-limited heating and cooling with on/off hysteresis control.
- Carried the bed state between cycles and evaluated fifth-cycle performance.
- Used SciPy COBYLA to minimise heater specific energy consumption (SEC), subject to at least 90% of baseline CO2 productivity.
- Investigated ambient conditions, operating conditions and time-step sensitivity.


## Optimisation Results

Constrained optimisation was performed using SciPy COBYLA to minimise heater-specific energy consumption (SEC), subject to maintaining at least 90% of baseline CO₂ productivity.

| Parameter | Baseline | Optimised |
|---|---:|---:|
| Regeneration temperature | 100°C | 95°C |
| Desorption duration | 20,000 s | 23,076.65 s |
| Heater SEC | 1.9951 kWh_th/kg CO₂ | 1.9490 kWh_th/kg CO₂ |
| CO₂ productivity | 6.3535 kg/h | 5.7183 kg/h |


### Energy–Productivity Trade-off

![TVSA energy-productivity trade-off](outputs/figures/sec_productivity_tradeoff.png)

The figure illustrates how regeneration temperature and desorption duration influence simulated heater-specific energy consumption and CO₂ productivity.

The selected operating point achieves a 2.31% reduction in heater SEC while maintaining approximately 90% of baseline productivity. The plotted operating conditions represent a discrete sensitivity study, not a globally optimised Pareto frontier.


### Dynamic Thermal Model

![TVSA dynamic bed temperature and equipment duties](outputs/figures/tvsa_cycle_temperature.png)

The dynamic thermal model incorporates a lumped bed energy balance, accounting for CO₂ and H₂O sorption heat, gas-flow sensible heat, and heating/cooling equipment duties.

Unlike the reference model, which uses a prescribed temperature response, the thermal extension calculates the bed temperature dynamically throughout the six-stage TVSA cycle.

The results demonstrate:

- **Heating:** Transient bed temperature response during regeneration, including deviations from the prescribed temperature trajectory.
- **Desorption:** Temperature regulation near 100°C through hysteresis-based heater control.
- **Cooling:** Rapid temperature reduction following desorption, accompanied by a high initial cooling duty.
- **Adsorption:** Gradual cooling towards ambient conditions.

The model predicts the thermal response under simplified equipment and heat-transfer assumptions. The calculated peak duties have not been independently validated against experimental measurements.


### MATLAB–Python Reference Model Comparison

![MATLAB–Python TVSA model comparison](outputs/figures/matlab_python_validation.png)

The Python implementation of the published reduced-order TVSA model was compared against reference MATLAB outputs under three ambient conditions:

- **Cold and dry:** 4°C, 20% RH
- **Moderate:** 20°C, 50% RH
- **Hot and humid:** 40°C, 90% RH

The comparison evaluates transient CO₂ loading and bed temperature throughout the six-stage TVSA cycle.

Across the three conditions, the comparison produced:

| Variable | RMSE range |
|---|---|
| CO₂ loading | 1.73–3.89 × 10⁻¹⁵ mol/kg |
| Bed temperature | 2.89–2.93 × 10⁻¹³ °C |

The near-identical trajectories demonstrate numerical agreement between the Python implementation and the stored MATLAB reference outputs.

**Validation scope:** This comparison verifies consistency with the reference data. It does not independently validate the extended thermal model or establish agreement with experimental measurements.


**Key finding:** The optimised operating conditions achieved a **2.31% reduction in simulated heater SEC**, while satisfying the minimum productivity constraint at the reported precision.

The results were reproduced using `run_scipy.py` and `check_timestep.py` with a 1-second simulation timestep.

**Limitations:** The solution is a local optimum under simplified thermal assumptions. Productivity is very close to the constraint boundary, numerical timestep convergence has not been demonstrated, and the reported SEC excludes vacuum, fan and cooling utility energy.

Further discussion is provided in [Results and Interpretation](project_result.md).



## Installation and Reproducibility

### Requirements

- Python 3.10 or later
- NumPy
- SciPy
- Pandas
- Matplotlib

The required Python packages are listed in `requirements.txt`.

### 1. Clone the Repository

```bash
git clone https://github.com/uhrenhwan/direct-air-capture-modelling.git
cd direct-air-capture-modelling
```

### 2. Set Up the Python Environment

On Windows (PowerShell):

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

On macOS/Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

The commands below use the Windows virtual-environment Python executable. On macOS/Linux, replace it with `.venv/bin/python`.

### 3. Run the Simulations

**Reduced-order reference TVSA model:**

```powershell
.\.venv\Scripts\python.exe run_baseline.py
```

**Dynamic thermal model:**

```powershell
.\.venv\Scripts\python.exe run_cyclic_thermal.py
```

**Constrained process optimisation:**

```powershell
.\.venv\Scripts\python.exe run_scipy.py
```

The optimisation uses SciPy COBYLA to minimise simulated heater-specific energy consumption while maintaining at least 90% of baseline CO₂ productivity.

The reported optimisation was performed using a 1-second simulation timestep.

### 4. Reproduce the Engineering Figures

Run the following scripts individually:

| Script | Output |
|---|---|
| `plot_optimisation.py` | Baseline vs optimised heater SEC |
| `run_baseline.py` | Reduced-order TVSA cycle dynamics |
| `run_thermal.py` | Adsorption thermal energy balance |
| `plot_cycle_temperature.py` | Dynamic bed temperature and equipment duties |
| `plot_tradeoff.py` | Energy–productivity trade-off |
| `plot_sensitivities.py` | Operating and ambient sensitivity analysis |
| `plot_model_validation.py` | MATLAB–Python reference model comparison |

For example:

```powershell
.\.venv\Scripts\python.exe plot_tradeoff.py
```

Generated figures are saved in `outputs/figures/`. Supporting numerical results are saved in `outputs/data/`.

Some plotting scripts perform multiple five-cycle simulations and may require substantial computation time.

### 5. Verify the Results

To verify the recorded optimisation candidate:

```powershell
.\.venv\Scripts\python.exe check_timestep.py
```

To regenerate the reference-model outputs and compare them with the stored MATLAB reference data:

```powershell
.\.venv\Scripts\python.exe run_baseline.py
.\.venv\Scripts\python.exe validate_matlab_parity.py
```

These checks evaluate numerical consistency, not experimental accuracy.

### Reproducibility Notes

- The reported optimisation results use a 1-second timestep.
- The model evaluates performance after five consecutive TVSA cycles.
- Numerical timestep convergence and cyclic steady-state convergence have not been fully established.
- The optimisation is a local constrained search; global optimality is not guaranteed.
- Dependency minimum versions are specified in `requirements.txt`, but an exact environment lockfile has not yet been established.

```



## Modelling Methodology

### 1. Reduced-Order TVSA Process Model

The process is represented by a six-stage temperature-vacuum swing adsorption (TVSA) cycle:

**Vacuum → Heating → Desorption → Cooling → Pressurisation → Adsorption**

The Python implementation reproduces the published reduced-order adsorption framework, including:

- CO₂ and H₂O adsorption equilibrium relationships.
- Linear driving force (LDF) adsorption/desorption kinetics.
- Fast and slow CO₂ adsorption states and free/bound H₂O loading states.
- Gas-phase CO₂ inventory and flow balances.
- Prescribed pressure, humidity and temperature trajectories.

The model tracks six state variables: four sorbent loading states, gas-phase CO₂ mole fraction and bed temperature.

The reference implementation uses a prescribed temperature response, providing a basis for comparison with the subsequent thermal extension.

### 2. Dynamic Thermal Model

The model was extended by replacing the prescribed bed temperature response with a lumped dynamic energy balance:

$$
C_{\mathrm{eff}}\frac{dT_{\mathrm{bed}}}{dt}
=
\dot Q_{\mathrm{heater}}
+
\dot Q_{\mathrm{sorption}}
+
\dot Q_{\mathrm{gas}}
-
\dot Q_{\mathrm{cooler}}
-
\dot Q_{\mathrm{loss}}
$$

where:

- $C_{\mathrm{eff}}$: effective bed heat capacity [J/K]
- $\dot Q_{\mathrm{heater}}$: heater input [W]
- $\dot Q_{\mathrm{sorption}}$: heat released or absorbed during CO₂/H₂O sorption [W]
- $\dot Q_{\mathrm{gas}}$: sensible heat transfer associated with inlet gas [W]
- $\dot Q_{\mathrm{cooler}}$: cooling duty [W]
- $\dot Q_{\mathrm{loss}}$: heat loss to surroundings [W]

Sorption heat is calculated from the simulated loading rates:

$$
\dot Q_{\mathrm{sorption}}
=
m_s\left(
\Delta H_{\mathrm{CO_2}}\frac{dq_{\mathrm{CO_2}}}{dt}
+
\Delta H_{\mathrm{H_2O}}\frac{dq_{\mathrm{H_2O}}}{dt}
\right)
$$

Heating and cooling duties are constrained by equipment heat-transfer conductances and maximum thermal capacities. On/off hysteresis controllers regulate the bed temperature during the relevant operating stages.

The thermal balance is integrated through explicit timestep updates. The current implementation assumes zero hardware heat capacity and zero ambient heat-loss conductance.

### 3. Cycle Simulation and Performance Evaluation

Five consecutive TVSA cycles are simulated, transferring the final bed state of each cycle to the next.

Performance is evaluated using the fifth cycle. This is a fixed evaluation procedure rather than a demonstrated cyclic steady state.

CO₂ productivity is calculated as:

$$
P_{\mathrm{CO_2}}=
\frac{m_{\mathrm{CO_2,captured}}}{t_{\mathrm{cycle}}}
$$

Heater-specific energy consumption is calculated as:

$$
SEC_{\mathrm{heater}}=
\frac{E_{\mathrm{heater}}}{m_{\mathrm{CO_2,captured}}}
$$

The reported SEC represents heater thermal energy only. It excludes vacuum/fan electricity and cooling utility energy.

### 4. Constrained Process Optimisation

SciPy COBYLA was used to minimise simulated heater SEC by varying:

- Regeneration temperature
- Desorption duration

The optimisation was subject to:

$$
P_{\mathrm{CO_2}}\geq0.90P_{\mathrm{baseline}}
$$

Five initial guesses were used to explore candidate operating conditions. The reported solution is a local optimum within the model's specified assumptions and operating bounds.

### 5. Code Organisation

| Module | Responsibility |
|---|---|
| `parameters.py` | Process, sorbent and equipment parameters |
| `cycle_driver.py` | Six-stage operating schedule |
| `isotherms.py`, `kinetics.py` | Adsorption equilibrium and kinetics |
| `gas_balance.py` | Gas-phase inventory and flow calculations |
| `simulator.py` | Python reproduction of the reference model |
| `simulator_thermal.py` | Extended dynamic thermal simulation |
| `thermal.py` | Energy balance and equipment controllers |
| `optimisation.py` | Performance evaluation |
| `run_scipy.py` | Constrained numerical optimisation |

For a detailed explanation of the implementation, see the corresponding source files.


## Baseline example

![Baseline cycle profiles for cold/dry, moderate and hot/humid air](outputs/figures/figure2_python_reimplementation.png)

The reference run shows sorbent loadings, bed temperature, humidity and pressure through each stage. This figure uses the prescribed-temperature baseline.


## Model Assumptions and Limitations

The model is intended for preliminary process analysis and operating-condition optimisation rather than detailed equipment design.

### Thermal Modelling Assumptions

- **Lumped thermal model:** The bed is represented by a single temperature state, without resolving spatial temperature gradients.
- **Adiabatic assumption:** Heat loss to the surroundings is neglected (UA = 0).
- **Equipment thermal inertia:** Hardware heat capacity is assumed to be zero.
- **Simplified equipment representation:** Heater and cooler performance is represented using heat-transfer conductances, maximum duties and hysteresis-based control.

These assumptions reduce computational complexity but may influence predicted temperature transients and equipment duties.

### Numerical Limitations

- The optimisation results were evaluated using a 1-second simulation timestep.
- Timestep sensitivity was investigated, but numerical convergence has not been fully established.
- Five consecutive cycles are simulated, with performance evaluated using the fifth cycle. Periodic cyclic steady state is not explicitly verified.
- The predicted initial adsorption heat spike and peak cooling duty have not been independently validated.

### Optimisation Limitations

- SciPy COBYLA performs a local constrained optimisation; global optimality is not established.
- The reported 2.31% reduction in heater SEC is accompanied by an approximately 10% reduction in CO₂ productivity.
- The optimised operating point lies very close to the minimum productivity constraint.
- The objective considers heater thermal energy only. Vacuum electricity, fan electricity and cooling utility energy are excluded.

Consequently, the reported SEC reduction should not be interpreted as an equivalent reduction in total plant energy consumption.

### Validation Scope

The Python implementation of the reference TVSA model was compared against stored MATLAB reference outputs across three ambient conditions.

This demonstrates numerical agreement with the reference dataset but does not independently establish experimental accuracy.

The extended dynamic thermal model has not been validated against experimental DAC operating data.

### Future Improvements

Potential improvements include:

- Incorporating equipment thermal inertia and ambient heat losses.
- Performing numerical timestep and cyclic steady-state convergence studies.
- Extending energy accounting to vacuum, fan and cooling utilities.
- Validating predicted temperature profiles and equipment duties against experimental or independently established reference data.


## Data and source material

- `validation_reference/`: MATLAB CSV exports used for regression tests; retain these in GitHub.
- `outputs/figures/`: example baseline comparison figure.
- `outputs/data/summary_metrics.csv`: compact baseline results. Large regenerated time-series CSVs are ignored by Git.
- Root sensitivity CSVs: recorded study outputs; rerunning a sweep overwrites the corresponding file.
- `reference_matlab/`: upstream MATLAB source.

## Development checks

```powershell
.\.venv\Scripts\python.exe -m pip install -r requirements-dev.txt
.\.venv\Scripts\python.exe -m ruff check .
.\.venv\Scripts\python.exe -m black --check .
.\.venv\Scripts\python.exe -m pytest -q
```

GitHub Actions runs these checks on pushes and pull requests.

## Attribution

The reference adsorption framework, parameters and Figure 2 protocol come from:

An, K., Wald, D., Chen, K. & King, R. (2026). *Cycle-resolved dynamic modeling of atmospheric CO2 separation systems under time-varying ambient and electricity conditions*. Frontiers in Environmental Engineering 5:1881554. https://doi.org/10.3389/fenve.2026.1881554

The Python translation and subsequent thermal/optimisation work are described separately above. The upstream MATLAB files retain their original content and attribution. Their redistribution remains subject to the upstream licence or permission.
