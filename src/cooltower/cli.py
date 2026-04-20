"""
cooltower.cli — config-driven entry point for reproducible cooling-tower analysis.

Usage:
    python -m cooltower.cli --config config/default.yaml
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import yaml

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
from cooltower.psychrometrics import humidity_ratio, relative_humidity, specific_enthalpy


def load_config(path: Path) -> dict:
    with path.open() as f:
        return yaml.safe_load(f)


def run(config_path: Path) -> dict:
    cfg = load_config(config_path)
    random.seed(cfg["output"].get("random_seed", 0))

    a = cfg["measurements"]["air"]
    w = cfg["measurements"]["water"]

    T_db1, T_wb1 = a["T_db_in_C"], a["T_wb_in_C"]
    T_db2, T_wb2 = a["T_db_out_C"], a["T_wb_out_C"]
    T_w3, T_w4 = w["T_hot_in_C"], w["T_cold_out_C"]
    m_w3 = w["m_in_kg_s"]

    omega1 = humidity_ratio(T_db1, T_wb1)
    omega2 = humidity_ratio(T_db2, T_wb2)
    h1 = specific_enthalpy(T_db1, omega1)
    h2 = specific_enthalpy(T_db2, omega2)
    rh1 = relative_humidity(T_db1, omega1)

    m_air = solve_air_flow_rate(T_db1, T_wb1, T_db2, T_wb2, T_w3, T_w4, m_w3)

    inlet = CoolingTowerState(T_water=T_w3, T_db=T_db1, T_wb=T_wb1, m_water=m_w3, m_air=m_air)
    outlet = CoolingTowerState(T_water=T_w4, T_db=T_db2, T_wb=T_wb2, m_water=m_w3, m_air=m_air)
    eb = solve_energy_balance(inlet, outlet)

    rng = range_temperature(T_w3, T_w4)
    appr = approach_temperature(T_w4, T_wb1)
    lg = liquid_to_gas_ratio(m_w3, m_air)

    cc = cfg["control"]["fopdt"]
    model = FOPDTModel(K_p=cc["K_p"], tau_p=cc["tau_p_s"], theta=cc["theta_s"])
    pi_lam = tune_lambda(model)
    pi_zn = tune_ziegler_nichols(model)

    sp = cfg["control"]["setpoint_step_C"]
    t_end = cfg["control"]["t_end_s"]
    dt = cfg["control"]["dt_s"]
    t_lam, y_lam, _ = closed_loop_response(model, pi_lam, setpoint=sp, t_end=t_end, dt=dt)
    e_lam = [sp - y for y in y_lam]
    idx = performance_indices(t_lam, e_lam)

    metrics = {
        "config": str(config_path),
        "psychrometrics": {
            "omega_in_g_kg": omega1 * 1000,
            "omega_out_g_kg": omega2 * 1000,
            "h_in_kJ_kg": h1 / 1000,
            "h_out_kJ_kg": h2 / 1000,
            "rh_in": rh1,
        },
        "energy_balance": {
            "m_air_kg_s": m_air,
            "Q_air_kW": eb.Q_air / 1000,
            "Q_water_kW": eb.Q_water / 1000,
            "m_evap_kg_hr": eb.m_evap * 3600,
            "imbalance_W": eb.imbalance,
        },
        "performance": {
            "range_C": rng,
            "approach_C": appr,
            "L_over_G": lg,
        },
        "control": {
            "lambda": str(pi_lam),
            "ziegler_nichols": str(pi_zn),
            "ISE": idx["ISE"],
            "IAE": idx["IAE"],
            "ITAE": idx["ITAE"],
            "y_final": y_lam[-1],
        },
    }

    out_dir = Path(cfg["output"]["results_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = out_dir / cfg["output"]["metrics_name"]
    with metrics_path.open("w") as f:
        json.dump(metrics, f, indent=2, default=str)

    _render_figure(t_lam, y_lam, sp, out_dir / cfg["output"]["figure_name"])
    print(f"✓ metrics → {metrics_path}")
    print(f"✓ figure  → {out_dir / cfg['output']['figure_name']}")
    return metrics


def _render_figure(t, y, sp, out_path: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t, y, lw=2, label="Closed-loop response")
    ax.axhline(sp, ls="--", color="grey", label=f"setpoint = {sp}")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Output (°C)")
    ax.set_title("Cooltower — λ-tuned PI closed-loop response")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(prog="python -m cooltower.cli")
    parser.add_argument("--config", type=Path, default=Path("config/default.yaml"))
    args = parser.parse_args()
    run(args.config)


if __name__ == "__main__":
    main()
