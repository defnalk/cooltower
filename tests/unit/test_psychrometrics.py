"""Unit tests for cooltower.psychrometrics."""

from __future__ import annotations

import math

import pytest

from cooltower.constants import P_STD
from cooltower.psychrometrics import (
    degree_of_saturation,
    dew_point_temperature,
    evaporation_rate,
    humidity_ratio,
    humidity_ratio_from_rh,
    relative_humidity,
    saturation_pressure,
    specific_enthalpy,
    wet_bulb_temperature,
)


# ---------------------------------------------------------------------------
# saturation_pressure
# ---------------------------------------------------------------------------


class TestSaturationPressure:
    @pytest.mark.parametrize(
        "T_db, expected, tol_pct",
        [
            (0.0, 611.2, 0.5),    # Triple point ≈ 611.7 Pa
            (25.0, 3169.9, 0.3),  # Standard reference
            (100.0, 101_325.0, 2.0),  # Boiling point at 1 atm
            (-10.0, 259.9, 1.0),  # Sub-zero (ice surface)
        ],
    )
    def test_known_values(self, T_db: float, expected: float, tol_pct: float) -> None:
        result = saturation_pressure(T_db)
        assert abs(result - expected) / expected * 100 < tol_pct, (
            f"saturation_pressure({T_db}) = {result:.2f}, expected ~{expected:.2f}"
        )

    def test_monotonically_increasing(self) -> None:
        temps = list(range(-20, 60, 5))
        pressures = [saturation_pressure(float(T)) for T in temps]
        for i in range(len(pressures) - 1):
            assert pressures[i] < pressures[i + 1]

    def test_invalid_temperature_too_high(self) -> None:
        with pytest.raises(ValueError, match="outside the physically valid range"):
            saturation_pressure(250.0)

    def test_invalid_temperature_too_low(self) -> None:
        with pytest.raises(ValueError, match="outside the physically valid range"):
            saturation_pressure(-150.0)


# ---------------------------------------------------------------------------
# humidity_ratio_from_rh
# ---------------------------------------------------------------------------


class TestHumidityRatioFromRH:
    def test_saturated_air(self) -> None:
        """At RH=1 the result should match saturation."""
        omega_sat = humidity_ratio_from_rh(25.0, 1.0)
        assert omega_sat > 0.0
        # Approximately 0.020 kg/kg at 25 °C, 101 325 Pa
        assert 0.018 < omega_sat < 0.022

    def test_dry_air(self) -> None:
        assert humidity_ratio_from_rh(25.0, 0.0) == pytest.approx(0.0, abs=1e-9)

    def test_invalid_rh_above_one(self) -> None:
        with pytest.raises(ValueError, match="Relative humidity"):
            humidity_ratio_from_rh(25.0, 1.1)

    def test_invalid_rh_negative(self) -> None:
        with pytest.raises(ValueError, match="Relative humidity"):
            humidity_ratio_from_rh(25.0, -0.1)

    def test_invalid_pressure(self) -> None:
        with pytest.raises(ValueError, match="Pressure must be positive"):
            humidity_ratio_from_rh(25.0, 0.5, P=0.0)


# ---------------------------------------------------------------------------
# humidity_ratio (wet-bulb method)
# ---------------------------------------------------------------------------


class TestHumidityRatio:
    def test_saturated_at_equality(self) -> None:
        """When T_db == T_wb the air is saturated."""
        omega = humidity_ratio(25.0, 25.0)
        omega_sat = humidity_ratio_from_rh(25.0, 1.0)
        assert omega == pytest.approx(omega_sat, rel=1e-3)

    def test_wb_exceeds_db_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot exceed"):
            humidity_ratio(20.0, 25.0)

    @pytest.mark.parametrize("T_db,T_wb", [(30.0, 22.0), (15.0, 12.0), (40.0, 28.0)])
    def test_non_negative(self, T_db: float, T_wb: float) -> None:
        assert humidity_ratio(T_db, T_wb) >= 0.0

    def test_decreases_with_temperature_depression(self) -> None:
        """More depression → lower humidity ratio."""
        omega1 = humidity_ratio(30.0, 28.0)
        omega2 = humidity_ratio(30.0, 22.0)
        assert omega1 > omega2


# ---------------------------------------------------------------------------
# specific_enthalpy
# ---------------------------------------------------------------------------


class TestSpecificEnthalpy:
    def test_zero_humidity_dry_air(self) -> None:
        """At ω=0 enthalpy should equal cp_a * T."""
        from cooltower.constants import CP_AIR
        h = specific_enthalpy(20.0, 0.0)
        assert h == pytest.approx(CP_AIR * 20.0, rel=1e-6)

    def test_increases_with_humidity(self) -> None:
        h_dry = specific_enthalpy(25.0, 0.005)
        h_humid = specific_enthalpy(25.0, 0.020)
        assert h_humid > h_dry

    def test_increases_with_temperature(self) -> None:
        omega = 0.01
        h_low = specific_enthalpy(10.0, omega)
        h_high = specific_enthalpy(40.0, omega)
        assert h_high > h_low

    def test_negative_omega_raises(self) -> None:
        with pytest.raises(ValueError, match="cannot be negative"):
            specific_enthalpy(25.0, -0.001)

    @pytest.mark.parametrize(
        "T_db,omega,expected_kJ",
        [
            (0.0, 0.0, 0.0),           # Datum: both terms zero
            (25.0, 0.010, 50.72),      # Rogers & Mayhew approximate check
        ],
    )
    def test_approximate_values(self, T_db: float, omega: float, expected_kJ: float) -> None:
        h = specific_enthalpy(T_db, omega) / 1000.0
        assert abs(h - expected_kJ) < 1.5, f"h({T_db}°C, ω={omega}) = {h:.2f} kJ/kg"


# ---------------------------------------------------------------------------
# Round-trip: humidity_ratio → relative_humidity → humidity_ratio
# ---------------------------------------------------------------------------


class TestRoundTrip:
    @pytest.mark.parametrize("T_db,T_wb", [(20.0, 15.0), (35.0, 28.0), (10.0, 8.0)])
    def test_humidity_ratio_round_trip(self, T_db: float, T_wb: float) -> None:
        omega = humidity_ratio(T_db, T_wb)
        rh = relative_humidity(T_db, omega)
        omega_recovered = humidity_ratio_from_rh(T_db, rh)
        assert omega_recovered == pytest.approx(omega, rel=1e-3)

    @pytest.mark.parametrize("T_db,rh", [(20.0, 0.6), (30.0, 0.8), (15.0, 0.4)])
    def test_dew_point_below_dry_bulb(self, T_db: float, rh: float) -> None:
        omega = humidity_ratio_from_rh(T_db, rh)
        T_dp = dew_point_temperature(omega)
        assert T_dp < T_db


# ---------------------------------------------------------------------------
# evaporation_rate
# ---------------------------------------------------------------------------


class TestEvaporationRate:
    def test_basic(self) -> None:
        m_evap = evaporation_rate(omega_out=0.020, omega_in=0.010, m_air=0.5)
        assert m_evap == pytest.approx(0.005, rel=1e-6)

    def test_no_evaporation_when_no_delta_omega(self) -> None:
        assert evaporation_rate(0.010, 0.010, 1.0) == pytest.approx(0.0)

    def test_condensation_clipped_to_zero(self) -> None:
        """When omega_out < omega_in the result should be zero, not negative."""
        assert evaporation_rate(0.005, 0.010, 1.0) == pytest.approx(0.0)

    def test_negative_omega_raises(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            evaporation_rate(-0.01, 0.01, 1.0)

    def test_negative_m_air_raises(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            evaporation_rate(0.02, 0.01, -1.0)


# ---------------------------------------------------------------------------
# wet_bulb_temperature (iterative inversion)
# ---------------------------------------------------------------------------


class TestWetBulbTemperature:
    def test_round_trip_with_humidity_ratio(self) -> None:
        T_db, T_wb_true = 30.0, 22.0
        omega = humidity_ratio(T_db, T_wb_true)
        T_wb_recovered = wet_bulb_temperature(T_db, omega)
        assert T_wb_recovered == pytest.approx(T_wb_true, abs=0.5)

    def test_saturated_returns_dry_bulb(self) -> None:
        from cooltower.psychrometrics import humidity_ratio_from_rh
        T_db = 25.0
        omega_sat = humidity_ratio_from_rh(T_db, 1.0)
        T_wb = wet_bulb_temperature(T_db, omega_sat)
        assert T_wb == pytest.approx(T_db, abs=1.0)


# ---------------------------------------------------------------------------
# degree_of_saturation
# ---------------------------------------------------------------------------


class TestDegreeOfSaturation:
    def test_saturated_is_one(self) -> None:
        from cooltower.psychrometrics import humidity_ratio_from_rh
        omega_sat = humidity_ratio_from_rh(25.0, 1.0)
        mu = degree_of_saturation(25.0, omega_sat)
        assert mu == pytest.approx(1.0, rel=1e-3)

    def test_dry_air_is_zero(self) -> None:
        assert degree_of_saturation(25.0, 0.0) == pytest.approx(0.0, abs=1e-9)

    def test_between_zero_and_one(self) -> None:
        omega = humidity_ratio(30.0, 22.0)
        mu = degree_of_saturation(30.0, omega)
        assert 0.0 < mu < 1.0
