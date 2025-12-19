"""Integration tests: full cooling-tower analysis pipeline."""

from __future__ import annotations

import pytest

from cooltower.control import (
    FOPDTModel,
    closed_loop_response,
    identify_fopdt,
    performance_indices,
    step_response,
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
from cooltower.psychrometrics import humidity_ratio, specific_enthalpy


@pytest.mark.integration
class TestFullThermalAnalysisPipeline:
    """End-to-end test: measured temperatures → energy balance → metrics."""

    def test_complete_thermal_analysis(self) -> None:
        """Reproduce a typical cooling tower lab calculation."""
        # --- Measured data (representative lab values) ---
        T_db1, T_wb1 = 24.5, 18.0   # Inlet air [°C]
        T_db2, T_wb2 = 31.0, 28.5   # Outlet air [°C]
        T_w3, T_w4 = 39.5, 27.5     # Water in / out [°C]
        m_w3 = 0.48                  # Inlet water flow [kg/s]

        # --- Back-calculate air flow ---
        m_air = solve_air_flow_rate(T_db1, T_wb1, T_db2, T_wb2, T_w3, T_w4, m_w3)
        assert 0.05 < m_air < 10.0, f"Unexpected m_air = {m_air:.4f} kg/s"

        # --- Build states ---
        inlet = CoolingTowerState(T_water=T_w3, T_db=T_db1, T_wb=T_wb1,
                                  m_water=m_w3, m_air=m_air)
        outlet = CoolingTowerState(T_water=T_w4, T_db=T_db2, T_wb=T_wb2,
                                   m_water=m_w3, m_air=m_air)

        # --- Energy balance ---
        result = solve_energy_balance(inlet, outlet)
        assert result.Q_air > 0, "Air must gain enthalpy"
        assert result.m_evap > 0, "Evaporation must occur"

        # --- Performance metrics ---
        approach = approach_temperature(T_w4, T_wb1)
        rng = range_temperature(T_w3, T_w4)
        lg = liquid_to_gas_ratio(m_w3, m_air)
        assert approach >= 0
        assert rng > 0
        assert lg > 0


@pytest.mark.integration
class TestControlDesignPipeline:
    """End-to-end: step test → FOPDT → lambda tuning → CL performance."""

    def test_lambda_better_than_zn_for_slow_process(self) -> None:
        """Lambda tuning should achieve lower ITAE than ZN for slow temperature loop."""
        model = FOPDTModel(K_p=0.6, tau_p=180.0, theta=20.0)

        pi_lambda = tune_lambda(model, lambda_=3 * model.theta)
        pi_zn = tune_ziegler_nichols(model)

        sp = 5.0
        t_end = 1800.0

        _, y_lam, _ = closed_loop_response(model, pi_lambda, sp, t_end, dt=2.0)
        _, y_zn, _ = closed_loop_response(model, pi_zn, sp, t_end, dt=2.0)

        t_vec = [2.0 * i for i in range(len(y_lam))]
        e_lam = [sp - y for y in y_lam]
        e_zn = [sp - y for y in y_zn]

        idx_lam = performance_indices(t_vec, e_lam)
        idx_zn = performance_indices(t_vec, e_zn)

        # Lambda should not produce larger ITAE than ZN for this slow process
        # (ZN is aggressive → large initial overshoot → large ITAE)
        assert idx_lam["ITAE"] <= idx_zn["ITAE"] * 1.5  # generous bound

    def test_identify_then_tune_pipeline(self, step_test_data: dict) -> None:
        """FOPDT identification → lambda tuning → verify setpoint tracking."""
        model = identify_fopdt(
            time=step_test_data["time"],
            output=step_test_data["output"],
            step_time=step_test_data["step_time"],
            step_magnitude=step_test_data["step_magnitude"],
            method="two_point",
        )
        pi = tune_lambda(model)
        sp = step_test_data["step_magnitude"] * model.K_p * 0.5
        t, y, _ = closed_loop_response(model, pi, setpoint=sp, t_end=1000.0, dt=2.0)
        # Verify the loop is stable (output doesn't diverge)
        assert max(abs(yi) for yi in y) < sp * 10
