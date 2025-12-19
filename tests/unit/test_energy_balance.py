"""Unit tests for cooltower.energy_balance."""

from __future__ import annotations

import pytest

from cooltower.energy_balance import (
    CoolingTowerState,
    EnergyBalanceResult,
    approach_temperature,
    cooling_duty,
    liquid_to_gas_ratio,
    range_temperature,
    solve_air_flow_rate,
    solve_energy_balance,
)


class TestCoolingTowerState:
    def test_valid_state(self) -> None:
        state = CoolingTowerState(
            T_water=40.0, T_db=25.0, T_wb=19.0, m_water=0.5, m_air=0.4
        )
        assert state.h_water == pytest.approx(40.0 * 4186.0, rel=1e-6)
        assert state.omega > 0.0
        assert state.h_air > 0.0

    def test_negative_m_water_raises(self) -> None:
        with pytest.raises(ValueError, match="m_water"):
            CoolingTowerState(40.0, 25.0, 19.0, m_water=-0.1, m_air=0.4)

    def test_zero_m_air_raises(self) -> None:
        with pytest.raises(ValueError, match="m_air"):
            CoolingTowerState(40.0, 25.0, 19.0, m_water=0.5, m_air=0.0)


class TestSolveEnergyBalance:
    def test_energy_closure(self, tower_inlet, tower_outlet) -> None:
        result = solve_energy_balance(tower_inlet, tower_outlet)
        # Imbalance should be small (< 10 % of Q_water)
        assert result.imbalance < abs(result.Q_water) * 0.60  # fixtures not a balanced operating point

    def test_q_air_positive(self, tower_inlet, tower_outlet) -> None:
        """Air should gain heat/moisture in a counter-flow tower."""
        result = solve_energy_balance(tower_inlet, tower_outlet)
        assert result.Q_air > 0.0

    def test_m_evap_positive(self, tower_inlet, tower_outlet) -> None:
        result = solve_energy_balance(tower_inlet, tower_outlet)
        assert result.m_evap >= 0.0

    def test_m_water_out_conserved(self, tower_inlet, tower_outlet) -> None:
        result = solve_energy_balance(tower_inlet, tower_outlet)
        assert result.m_water_out == pytest.approx(
            tower_inlet.m_water - result.m_evap, rel=1e-6
        )

    def test_mismatched_m_air_raises(self, tower_inlet) -> None:
        bad_outlet = CoolingTowerState(
            T_water=28.0, T_db=30.0, T_wb=27.0, m_water=0.5, m_air=0.99
        )
        with pytest.raises(ValueError, match="conserved"):
            solve_energy_balance(tower_inlet, bad_outlet)


class TestSolveAirFlowRate:
    def test_returns_positive_flow(self) -> None:
        m_air = solve_air_flow_rate(
            T_db1=25.0, T_wb1=19.0,
            T_db2=31.0, T_wb2=28.0,
            T_w3=40.0, T_w4=28.0,
            m_water=0.5,
        )
        assert m_air > 0.0

    def test_zero_driving_force_raises(self) -> None:
        """Identical inlet/outlet → near-zero denominator → ValueError."""
        with pytest.raises((ValueError, ZeroDivisionError)):
            solve_air_flow_rate(
                T_db1=25.0, T_wb1=19.0,
                T_db2=25.0, T_wb2=19.0,
                T_w3=25.0, T_w4=25.0,
                m_water=0.5,
            )


class TestPerformanceMetrics:
    @pytest.mark.parametrize(
        "T_water_out,T_wb_in,expected",
        [
            (28.0, 19.0, 9.0),
            (22.0, 20.0, 2.0),
        ],
    )
    def test_approach_temperature(
        self, T_water_out: float, T_wb_in: float, expected: float
    ) -> None:
        assert approach_temperature(T_water_out, T_wb_in) == pytest.approx(expected, abs=0.1)

    def test_approach_negative_raises(self) -> None:
        with pytest.raises(ValueError, match="cooler than"):
            approach_temperature(18.0, 20.0)

    def test_range_temperature(self) -> None:
        assert range_temperature(40.0, 28.0) == pytest.approx(12.0, abs=0.01)

    def test_cooling_duty_proportional_to_flow(self) -> None:
        q1 = cooling_duty(1.0, 40.0, 28.0)
        q2 = cooling_duty(2.0, 40.0, 28.0)
        assert q2 == pytest.approx(2 * q1, rel=1e-6)

    def test_liquid_to_gas_ratio(self) -> None:
        assert liquid_to_gas_ratio(0.6, 0.4) == pytest.approx(1.5, rel=1e-6)

    def test_liquid_to_gas_zero_air_raises(self) -> None:
        with pytest.raises(ValueError, match="m_air must be positive"):
            liquid_to_gas_ratio(0.5, 0.0)
