# Thermal model and optimisation results

The study minimises heater specific energy consumption (SEC) while maintaining at least 90% of the reference CO2 productivity. SciPy COBYLA searches over regeneration temperature (80-120 C) and desorption duration (14,000-26,000 s), using five starting points.

Each evaluation runs five consecutive thermal cycles at 20 C and 50% relative humidity, passing the final state into the next cycle. Metrics use the fifth cycle. The optimisation and sensitivity scripts use a 1 s time step.

## Recorded optimisation result

| Quantity | Baseline | Candidate |
| --- | ---: | ---: |
| Regeneration temperature [C] | 100 | approximately 95 |
| Desorption duration [s] | 20,000 | approximately 23,077 |
| CO2 productivity [kg/h] | 6.354 | 5.718 |
| Heater SEC [kWh thermal/kg CO2] | 1.995 | 1.949 |

The recorded SEC reduction is 2.31%. Lower regeneration temperature and longer desorption reduced heater energy per captured mass, with a productivity penalty close to the constraint boundary. These values come from the earlier project study; the full solver log and exact dependency versions were not archived. Re-run `run_scipy.py` to establish results for the current code and environment.

## Thermal extension

The bed temperature follows a lumped energy balance including sorbent/gas heat capacity, CO2 and water sorption heat, capacity-limited heating and cooling, hysteresis control, and sensible heat transfer during adsorption. `equipment_sizing.py` provides preliminary sensible-only sizing estimates, rather than a detailed exchanger design.

## Sensitivity studies

`ambient_sensitivity.py` evaluates 15, 20 and 25 C at 30%, 50% and 70% RH using the candidate operating point. `operating_sensitivity.py` evaluates nearby regeneration temperatures and desorption durations. The checked-in CSVs contain the recorded outputs.

The ambient study suggests a strong humidity effect on heater SEC and a temperature effect on productivity within the sampled range. The operating study illustrates that the lowest SEC and largest productivity margin need not occur at the same point. These observations describe the sampled model outputs, rather than general DAC performance laws.

## Interpretation limits

Heater SEC excludes auxiliary electricity and cooling utility energy. The model has simplified thermal equipment and currently assumes no hardware heat capacity or ambient heat loss. Five cycles do not establish periodic convergence. The search identifies a candidate local solution; it does not establish a global optimum. MATLAB numerical agreement was checked for the separate prescribed-temperature baseline, not this thermal extension.
