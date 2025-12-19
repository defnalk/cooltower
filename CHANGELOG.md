# Changelog

All notable changes to `cooltower` will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Planned
- Merkel number (NTU) calculation and tower sizing.
- Cooling tower characteristic curve fitting (Chebyshev polynomial).
- Drift loss estimation.
- `matplotlib` plotting helpers (`cooltower.plotting`).

---

## [0.1.0] — 2024-06-01

### Added
- `cooltower.psychrometrics` module:
  - `saturation_pressure` — Buck (1981) formulation with ice-surface branch below 0 °C.
  - `humidity_ratio_from_rh` — from relative humidity and barometric pressure.
  - `humidity_ratio` — from dry- and wet-bulb temperatures (Sprung psychrometric formula).
  - `specific_enthalpy` — moist-air enthalpy relative to 0 °C liquid datum.
  - `relative_humidity`, `dew_point_temperature`, `wet_bulb_temperature`.
  - `degree_of_saturation`, `evaporation_rate`.
- `cooltower.energy_balance` module:
  - `CoolingTowerState` dataclass with derived psychrometric properties.
  - `EnergyBalanceResult` frozen dataclass.
  - `solve_energy_balance` — steady-flow energy and mass balance with closure check.
  - `solve_air_flow_rate` — back-calculation of ṁₐ from measured temperatures.
  - `cooling_duty`, `approach_temperature`, `range_temperature`, `liquid_to_gas_ratio`.
- `cooltower.control` module:
  - `FOPDTModel` and `PIParameters` frozen dataclasses.
  - `identify_fopdt` — tangent and two-point step-test identification.
  - `tune_lambda` (IMC / lambda), `tune_ziegler_nichols`, `tune_cohen_coon`.
  - `step_response` — open-loop FOPDT simulation.
  - `closed_loop_response` — velocity-form PI simulation with optional load disturbance.
  - `performance_indices` — ISE, IAE, ITAE.
- `cooltower.constants` — all physical constants (SI), no hardcoded values in modules.
- Full pytest suite: unit tests (psychrometrics, energy balance, control) + integration pipeline tests.
- GitHub Actions CI: lint (ruff + mypy) + test matrix (Python 3.10, 3.11, 3.12).
- `Makefile` with `install`, `test`, `lint`, `format`, `typecheck`, `example`, `clean` targets.
- `examples/basic_analysis.py` — end-to-end lab analysis script.

[Unreleased]: https://github.com/defnalk/cooltower/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/defnalk/cooltower/releases/tag/v0.1.0
