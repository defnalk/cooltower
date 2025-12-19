"""Basic cooling tower analysis example.

Demonstrates the full cooltower workflow:
  1. Compute psychrometric properties from measured temperatures.
  2. Back-calculate the air-side mass-flow rate.
  3. Close the energy balance and extract performance metrics.
  4. Identify a process model from a step test and design a PI controller.

Run:
    python examples/basic_analysis.py
"""

import math

from cooltower.control import (
    FOPDTModel,
    closed_loop_response,
    performance_indices,
    tune_lambda,
    tune_ziegler_nichols,
)
from cooltower.energy_balance import (
    CoolingTowerState,
    approach_temperature,
    liquid_to_gas_ratio,
    range_temperature,
    solve_air_flow_rate,
    solve_energy_balance,
)
from cooltower.psychrometrics import (
    evaporation_rate,
    humidity_ratio,
    relative_humidity,
    specific_enthalpy,
)

# ── 1. Measured operating point ──────────────────────────────────────────────

print("=" * 60)
print("  COOLING TOWER LAB — ENERGY BALANCE ANALYSIS")
print("=" * 60)

# Air-side measurements
T_db1, T_wb1 = 24.5, 18.0   # Inlet air dry-bulb / wet-bulb [°C]
T_db2, T_wb2 = 31.0, 28.5   # Outlet air [°C]

# Water-side measurements
T_w3, T_w4 = 39.5, 27.5     # Hot inlet / cold outlet [°C]
m_w3 = 0.48                  # Inlet water mass-flow rate [kg/s]

# ── 2. Psychrometric properties ──────────────────────────────────────────────

omega1 = humidity_ratio(T_db1, T_wb1)
omega2 = humidity_ratio(T_db2, T_wb2)
h1 = specific_enthalpy(T_db1, omega1)
h2 = specific_enthalpy(T_db2, omega2)
rh1 = relative_humidity(T_db1, omega1)

print(f"\n  Inlet air:   T_db={T_db1}°C, T_wb={T_wb1}°C")
print(f"               ω₁ = {omega1*1000:.2f} g/kg_da,  φ₁ = {rh1*100:.1f}%")
print(f"               h₁ = {h1/1000:.2f} kJ/kg_da")
print(f"\n  Outlet air:  T_db={T_db2}°C, T_wb={T_wb2}°C")
print(f"               ω₂ = {omega2*1000:.2f} g/kg_da")
print(f"               h₂ = {h2/1000:.2f} kJ/kg_da")

# ── 3. Back-calculate air mass-flow rate ─────────────────────────────────────

m_air = solve_air_flow_rate(T_db1, T_wb1, T_db2, T_wb2, T_w3, T_w4, m_w3)
print(f"\n  Back-calculated ṁₐ = {m_air:.4f} kg/s")

# ── 4. Energy balance ────────────────────────────────────────────────────────

inlet = CoolingTowerState(
    T_water=T_w3, T_db=T_db1, T_wb=T_wb1, m_water=m_w3, m_air=m_air
)
outlet = CoolingTowerState(
    T_water=T_w4, T_db=T_db2, T_wb=T_wb2, m_water=m_w3, m_air=m_air
)

result = solve_energy_balance(inlet, outlet)
print(f"\n  Energy balance:")
print(f"    Q_air   = {result.Q_air/1000:+.3f} kW  (air gains heat)")
print(f"    Q_water = {result.Q_water/1000:+.3f} kW  (water rejects heat)")
print(f"    ṁ_evap  = {result.m_evap*3600:.2f} kg/hr")
print(f"    Closure imbalance = {result.imbalance:.2f} W")

# ── 5. Performance metrics ───────────────────────────────────────────────────

approach = approach_temperature(T_w4, T_wb1)
rng = range_temperature(T_w3, T_w4)
lg = liquid_to_gas_ratio(m_w3, m_air)

print(f"\n  Performance metrics:")
print(f"    Range       = {rng:.1f} °C")
print(f"    Approach    = {approach:.1f} °C")
print(f"    L/G ratio   = {lg:.3f}")

# ── 6. Process control design ────────────────────────────────────────────────

print("\n" + "=" * 60)
print("  PROCESS CONTROL DESIGN — LAMBDA TUNING")
print("=" * 60)

# Identified FOPDT model from open-loop step test on water flow
model = FOPDTModel(K_p=0.82, tau_p=145.0, theta=18.0)
print(f"\n  FOPDT model: K_p={model.K_p}, τ_p={model.tau_p} s, θ={model.theta} s")
print(f"  Dead-time ratio θ/τ_p = {model.theta/model.tau_p:.3f}")

# Lambda tuning (preferred for slow, noisy temperature loops)
pi_lam = tune_lambda(model)
pi_zn = tune_ziegler_nichols(model)

print(f"\n  {pi_lam}")
print(f"  {pi_zn}")
print(
    "\n  → Lambda chosen: conservative gain avoids windup against valve envelope.\n"
    "    Ziegler–Nichols would require observing sustained oscillations\n"
    "    and produces ~25% overshoot — unacceptable for thermal control."
)

# Closed-loop simulation
sp = 5.0
t_lam, y_lam, u_lam = closed_loop_response(model, pi_lam, setpoint=sp, t_end=1200.0, dt=2.0)
e_lam = [sp - y for y in y_lam]
idx_lam = performance_indices(t_lam, e_lam)

print(f"\n  Lambda CL performance (setpoint step = {sp} °C):")
print(f"    ISE  = {idx_lam['ISE']:.1f}")
print(f"    IAE  = {idx_lam['IAE']:.1f}")
print(f"    ITAE = {idx_lam['ITAE']:.1f}")
print(f"    Steady-state output = {y_lam[-1]:.4f} (target {sp})")
print()
