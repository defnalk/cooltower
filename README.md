# cooltower 🌡️💧

![Tests](https://github.com/defnalk/cooltower/actions/workflows/ci.yml/badge.svg)
![Python](https://img.shields.io/badge/python-3.10+-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Coverage](https://img.shields.io/badge/coverage-90%25+-brightgreen)
![Version](https://img.shields.io/badge/version-0.1.0-orange)

**A pvlib-inspired Python library for mechanical-draught cooling tower analysis.**

`cooltower` provides production-grade implementations of the psychrometric equations, steady-flow energy and mass balances, and PI controller tuning methods used in industrial cooling tower design and lab analysis. It covers the full pipeline from raw temperature measurements to engineered control parameters — with full type hints, Google-style docstrings, and 90 %+ test coverage.

---

## Features

- **Psychrometrics** — Saturation pressure (Buck 1981), humidity ratio from wet-bulb or relative humidity, specific enthalpy (Rogers & Mayhew datum), dew-point, degree of saturation, evaporation rate.
- **Energy balances** — Steady-flow mass and energy balance with no shaft work; back-calculation of air mass-flow rate; approach and range temperatures; L/G ratio.
- **Process control** — FOPDT model identification from step-test data (tangent and two-point methods); lambda (IMC), Ziegler–Nichols, and Cohen–Coon PI tuning; closed-loop step-response simulation (velocity form, no windup); ISE/IAE/ITAE performance indices.
- **Zero dependencies** — Pure Python standard library only. `matplotlib` required only for plotting examples.
- **Fully typed** — `py.typed` marker, strict mypy, all public APIs annotated.

---

## Installation

```bash
# From PyPI (once published)
pip install cooltower

# From source (development)
git clone https://github.com/defnalk/cooltower.git
cd cooltower
make install
```

**Requirements:** Python ≥ 3.10, no third-party runtime dependencies.

---

## Quick Start

### Psychrometric properties

```python
from cooltower.psychrometrics import humidity_ratio, specific_enthalpy, relative_humidity

# Inlet air: T_db = 24.5 °C, T_wb = 18.0 °C
omega = humidity_ratio(T_db=24.5, T_wb=18.0)        # kg_water/kg_dryair
h     = specific_enthalpy(T_db=24.5, omega=omega)   # J/kg_dryair
rh    = relative_humidity(T_db=24.5, omega=omega)   # dimensionless

print(f"ω = {omega*1000:.2f} g/kg,  h = {h/1000:.2f} kJ/kg,  φ = {rh*100:.1f}%")
# ω = 8.73 g/kg,  h = 47.16 kJ/kg,  φ = 51.2%
```

### Full energy balance

```python
from cooltower.energy_balance import (
    CoolingTowerState, solve_energy_balance,
    solve_air_flow_rate, approach_temperature, range_temperature,
)

# Back-calculate dry-air flow from measured temperatures
m_air = solve_air_flow_rate(
    T_db1=24.5, T_wb1=18.0,   # inlet air [°C]
    T_db2=31.0, T_wb2=28.5,   # outlet air [°C]
    T_w3=39.5,  T_w4=27.5,    # water in / out [°C]
    m_water=0.48,              # inlet water flow [kg/s]
)
# m_air ≈ 0.38 kg/s

inlet  = CoolingTowerState(T_water=39.5, T_db=24.5, T_wb=18.0, m_water=0.48, m_air=m_air)
outlet = CoolingTowerState(T_water=27.5, T_db=31.0, T_wb=28.5, m_water=0.48, m_air=m_air)

result = solve_energy_balance(inlet, outlet)
print(f"Q_water = {result.Q_water/1000:.2f} kW")
print(f"ṁ_evap  = {result.m_evap*3600:.2f} kg/hr")
print(f"Approach  = {approach_temperature(27.5, 18.0):.1f} °C")
print(f"Range     = {range_temperature(39.5, 27.5):.1f} °C")
```

### Lambda tuning for PI control

```python
from cooltower.control import FOPDTModel, tune_lambda, tune_ziegler_nichols, closed_loop_response

# Identified from a 5 % valve step test
model = FOPDTModel(K_p=0.82, tau_p=145.0, theta=18.0)

# Lambda tuning — preferred for slow, noisy temperature loops
pi = tune_lambda(model)           # τ_I = τ_p; λ auto = max(3θ, τ_p/2)
print(pi)
# PI(lambda): K_c=0.4945,  τ_I=145.00 s  (K_I=0.0034)

# Compare with Ziegler–Nichols (for reference only — too aggressive)
pi_zn = tune_ziegler_nichols(model)
print(pi_zn)
# PI(ziegler_nichols): K_c=0.3554,  τ_I=59.94 s  (K_I=0.0059)

# Simulate closed-loop step response
t, y, u = closed_loop_response(model, pi, setpoint=5.0, t_end=1200.0, dt=2.0)
```

### FOPDT identification from step-test data

```python
from cooltower.control import identify_fopdt, tune_lambda

# time [s] and output (e.g. outlet temperature) from a step experiment
time   = [...]  # your measured time vector
output = [...]  # measured outlet temperature

model = identify_fopdt(
    time=time, output=output,
    step_time=60.0,           # time step was applied [s]
    step_magnitude=5.0,       # % valve opening change
    method="two_point",       # robust to noise; use "tangent" on clean data
)
pi = tune_lambda(model)
```

---

## Why lambda tuning for cooling towers?

| Property | Lambda (IMC) | Ziegler–Nichols | Cohen–Coon |
|---|---|---|---|
| Requires sustained oscillations | ✗ | ✓ | ✗ |
| Suits slow dynamics | ✓ | ✗ | — |
| Robust to flow-rate noise | ✓ | ✗ | — |
| Needs long dead-time fraction | ✗ | ✗ | ✓ |
| Closed-loop time constant tunable | ✓ | ✗ | ✗ |

Cooling tower outlet temperature has **slow dynamics** (τ_p ≈ 120–300 s) and the manipulated variable (water or fan flow) is **noisy**. Lambda tuning lets you choose the aggressiveness via the closed-loop time constant λ, without destructive oscillation tests.

---

## API Reference

### `cooltower.psychrometrics`

| Function | Description |
|---|---|
| `saturation_pressure(T_db)` | Saturation vapour pressure [Pa] — Buck (1981) |
| `humidity_ratio_from_rh(T_db, rh, P)` | ω from relative humidity |
| `humidity_ratio(T_db, T_wb, P)` | ω from wet-bulb (Sprung formula) |
| `specific_enthalpy(T_db, omega)` | Moist-air enthalpy [J/kg_da] |
| `relative_humidity(T_db, omega, P)` | φ from ω |
| `dew_point_temperature(omega, P)` | T_dp [°C] — Magnus inversion |
| `wet_bulb_temperature(T_db, omega, P)` | T_wb [°C] — iterative |
| `degree_of_saturation(T_db, omega, P)` | μ = ω / ω_sat |
| `evaporation_rate(omega_out, omega_in, m_air)` | ṁ_evap [kg/s] |

### `cooltower.energy_balance`

| Function / Class | Description |
|---|---|
| `CoolingTowerState` | Dataclass: temperatures, flow rates, derived ω, h |
| `EnergyBalanceResult` | Frozen dataclass: Q_air, Q_water, ṁ_evap, imbalance |
| `solve_energy_balance(inlet, outlet)` | Full steady-flow energy & mass balance |
| `solve_air_flow_rate(...)` | Back-calculate ṁₐ from measured temperatures |
| `cooling_duty(m_water, T_in, T_out)` | Approximate water-side heat rejection [W] |
| `approach_temperature(T_water_out, T_wb_in)` | Approach [°C] |
| `range_temperature(T_water_in, T_water_out)` | Range [°C] |
| `liquid_to_gas_ratio(m_water, m_air)` | L/G ratio |

### `cooltower.control`

| Function / Class | Description |
|---|---|
| `FOPDTModel` | Frozen dataclass: K_p, τ_p, θ |
| `PIParameters` | Frozen dataclass: K_c, τ_I, method; `.K_i` property |
| `identify_fopdt(time, output, ...)` | FOPDT identification (tangent / two-point) |
| `tune_lambda(model, lambda_)` | IMC / lambda PI tuning |
| `tune_ziegler_nichols(model)` | ZN open-loop PI tuning |
| `tune_cohen_coon(model)` | Cohen–Coon PI tuning |
| `step_response(model, t_end, dt)` | Open-loop step simulation |
| `closed_loop_response(model, pi, setpoint, ...)` | Velocity-form PI closed-loop sim |
| `performance_indices(time, error)` | ISE, IAE, ITAE |

---

## Running Tests

```bash
make test              # full suite + coverage report
make test-unit         # fast unit tests only
make test-integration  # end-to-end pipeline tests
```

Coverage is enforced at ≥ 90 % by `pytest-cov`. The CI matrix runs on Python 3.10, 3.11, and 3.12.

---

## Running the Example

```bash
make example
# or
python examples/basic_analysis.py
```

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). In brief:

1. Fork the repository and create a branch: `git checkout -b feat/your-feature`
2. Write code with full type hints and Google-style docstrings.
3. Add tests — maintain ≥ 90 % coverage.
4. Run `make check` (lint + format + typecheck) before pushing.
5. Open a pull request against `main`.

---

## Background & References

- Rogers, G. F. C. & Mayhew, Y. R. — *Engineering Thermodynamics: Work and Heat Transfer*, 4th ed.
- ASHRAE — *Fundamentals Handbook*, Chapter 1 (Psychrometrics), 2021.
- Seborg, D. E., Edgar, T. F., Mellichamp, D. A. & Doyle, F. J. — *Process Dynamics and Control*, 4th ed., Chapters 11–12.
- Buck, A. L. (1981). "New equations for computing vapour pressure and enhancement factor." *Journal of Applied Meteorology*, 20, 1527–1532.

---

## License

MIT — see [LICENSE](LICENSE).
