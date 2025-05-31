"""Shared pytest fixtures for cooltower test suite."""

from __future__ import annotations

import math

import pytest

from cooltower.control import FOPDTModel, PIParameters
from cooltower.energy_balance import CoolingTowerState
from cooltower.psychrometrics import humidity_ratio, specific_enthalpy

# ---------------------------------------------------------------------------
# Psychrometric fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def standard_conditions() -> dict:
    """Standard atmospheric conditions for a temperate summer day."""
    return {"T_db": 25.0, "T_wb": 19.0, "P": 101_325.0}


@pytest.fixture
def humid_conditions() -> dict:
    """High-humidity tropical inlet conditions."""
    return {"T_db": 32.0, "T_wb": 28.0, "P": 101_325.0}


@pytest.fixture
def omega_standard(standard_conditions: dict) -> float:
    """Humidity ratio at standard conditions."""
    return humidity_ratio(
        T_db=standard_conditions["T_db"],
        T_wb=standard_conditions["T_wb"],
        P=standard_conditions["P"],
    )


@pytest.fixture
def enthalpy_standard(standard_conditions: dict, omega_standard: float) -> float:
    """Specific enthalpy at standard conditions."""
    return specific_enthalpy(T_db=standard_conditions["T_db"], omega=omega_standard)


# ---------------------------------------------------------------------------
# Energy balance fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tower_inlet() -> CoolingTowerState:
    """Cooling tower inlet state (air enters cool, water enters hot)."""
    return CoolingTowerState(
        T_water=40.0,
        T_db=25.0,
        T_wb=19.0,
        m_water=0.5,
        m_air=0.4,
    )


@pytest.fixture
def tower_outlet() -> CoolingTowerState:
    """Cooling tower outlet state (air exits humid, water exits cooler)."""
    return CoolingTowerState(
        T_water=28.0,
        T_db=30.0,
        T_wb=27.0,
        m_water=0.5,
        m_air=0.4,
    )


# ---------------------------------------------------------------------------
# Control fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def fopdt_typical() -> FOPDTModel:
    """Typical FOPDT model for a cooling tower temperature loop."""
    return FOPDTModel(K_p=0.75, tau_p=120.0, theta=15.0)


@pytest.fixture
def fopdt_fast() -> FOPDTModel:
    """Fast process with large dead time (challenging for Ziegler–Nichols)."""
    return FOPDTModel(K_p=1.2, tau_p=30.0, theta=25.0)


@pytest.fixture
def fopdt_slow() -> FOPDTModel:
    """Slow, noisy process typical of outlet temperature control."""
    return FOPDTModel(K_p=0.5, tau_p=300.0, theta=8.0)


@pytest.fixture
def step_test_data(fopdt_typical: FOPDTModel) -> dict:
    """Synthetic step-test data generated from the typical FOPDT model."""
    from cooltower.control import step_response

    time, output = step_response(fopdt_typical, t_end=600.0, dt=2.0, step_magnitude=5.0)
    # Add small reproducible noise (deterministic, no random seed needed)
    noisy_output = [y + 0.05 * math.sin(i * 0.3) for i, y in enumerate(output)]
    return {
        "time": time,
        "output": noisy_output,
        "step_time": 0.0,
        "step_magnitude": 5.0,
        "model": fopdt_typical,
    }
