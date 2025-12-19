"""Steady-flow energy and mass balances for cooling tower systems.

The sign convention throughout is:
  - Subscripts 1/2 refer to the air stream (inlet / outlet).
  - Subscripts 3/4 refer to the water stream (inlet / outlet).
  - All enthalpies are in J per unit mass of **dry air** for the air side
    and J/kg for the water side.

The core identity is the steady-flow energy balance with no shaft work:

.. math::

    \\dot{m}_a h_1 + \\dot{m}_3 h_3 = \\dot{m}_a h_2 + \\dot{m}_4 h_4

Combined with the water-side continuity:

.. math::

    \\dot{m}_4 = \\dot{m}_3 - (\\omega_2 - \\omega_1) \\dot{m}_a

Reference: Rogers & Mayhew, *Engineering Thermodynamics: Work and Heat
Transfer*, 4th ed., §13.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

from cooltower.constants import CP_WATER, H_FG_0
from cooltower.psychrometrics import (
    evaporation_rate,
    humidity_ratio_from_rh,
    specific_enthalpy,
)

__all__ = [
    "CoolingTowerState",
    "EnergyBalanceResult",
    "solve_air_flow_rate",
    "solve_energy_balance",
    "cooling_duty",
    "approach_temperature",
    "range_temperature",
    "liquid_to_gas_ratio",
]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass
class CoolingTowerState:
    """Thermodynamic state at a single cross-section of a cooling tower.

    Args:
        T_water: Water temperature  [°C].
        T_db:    Air dry-bulb temperature  [°C].
        T_wb:    Air wet-bulb temperature  [°C].
        m_water: Water mass-flow rate  [kg/s].
        m_air:   Dry-air mass-flow rate  [kg/s].
        P:       Barometric pressure  [Pa].  Defaults to 101 325 Pa.
        rh:      Relative humidity of air (0–1).  Computed from *T_db* /
                 *T_wb* if not supplied.
    """

    T_water: float
    T_db: float
    T_wb: float
    m_water: float
    m_air: float
    P: float = 101_325.0
    rh: float | None = None

    # Derived fields — populated by __post_init__
    omega: float = field(init=False)
    h_air: float = field(init=False)
    h_water: float = field(init=False)

    def __post_init__(self) -> None:
        from cooltower.psychrometrics import humidity_ratio

        if self.m_water < 0:
            raise ValueError(f"m_water must be non-negative, got {self.m_water}.")
        if self.m_air <= 0:
            raise ValueError(f"m_air must be positive, got {self.m_air}.")

        self.omega = humidity_ratio(self.T_db, self.T_wb, self.P)
        self.h_air = specific_enthalpy(self.T_db, self.omega)
        # Water enthalpy relative to 0 °C liquid datum: h = cp * T
        self.h_water = CP_WATER * self.T_water


@dataclass(frozen=True)
class EnergyBalanceResult:
    """Results of a complete cooling tower energy-balance calculation.

    Attributes:
        Q_air:       Heat gained by the air stream  [W].
        Q_water:     Heat rejected by the water stream  [W].
        m_evap:      Evaporation rate  [kg/s].
        m_water_out: Outlet water mass-flow rate  [kg/s].
        imbalance:   Absolute energy imbalance |Q_air − Q_water|  [W].
                     Should be < 1 W for a well-converged balance.
    """

    Q_air: float
    Q_water: float
    m_evap: float
    m_water_out: float
    imbalance: float


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def solve_energy_balance(
    inlet: CoolingTowerState,
    outlet: CoolingTowerState,
) -> EnergyBalanceResult:
    """Evaluate the steady-flow energy balance across a cooling tower.

    The air and water sides are treated separately and the closure is
    checked.  Shaft work is neglected (fan power << thermal duty).

    Args:
        inlet:  State at the tower inlet (air side 1, water side 3).
        outlet: State at the tower outlet (air side 2, water side 4).

    Returns:
        An :class:`EnergyBalanceResult` containing heat duties, evaporation
        rate, outlet water flow, and closure imbalance.

    Raises:
        ValueError: If air mass-flow rates are inconsistent between inlet
            and outlet (dry-air flow must be conserved).
    """
    if abs(inlet.m_air - outlet.m_air) / (inlet.m_air + 1e-9) > 0.01:
        raise ValueError(
            f"Dry-air mass-flow rate should be conserved across the tower. "
            f"inlet.m_air={inlet.m_air:.4f}, outlet.m_air={outlet.m_air:.4f}."
        )

    m_a = inlet.m_air  # dry-air flow [kg/s]
    m_evap = evaporation_rate(outlet.omega, inlet.omega, m_a)
    m_water_out = inlet.m_water - m_evap

    if m_water_out < 0:
        raise ValueError(
            f"Computed outlet water flow m4 = {m_water_out:.4f} kg/s < 0. "
            "Evaporation rate exceeds inlet water flow — check inputs."
        )

    Q_air = m_a * (outlet.h_air - inlet.h_air)  # heat gained by air [W]
    Q_water = inlet.m_water * inlet.h_water - m_water_out * outlet.h_water  # heat lost by water [W]
    imbalance = abs(Q_air - Q_water)

    if imbalance > 100.0:
        logger.warning(
            "Energy balance imbalance = %.2f W (>100 W). "
            "Verify that inlet/outlet states are measured at the same steady-state.",
            imbalance,
        )

    logger.debug(
        "Energy balance: Q_air=%.1f W, Q_water=%.1f W, imbalance=%.2f W",
        Q_air,
        Q_water,
        imbalance,
    )
    return EnergyBalanceResult(
        Q_air=Q_air,
        Q_water=Q_water,
        m_evap=m_evap,
        m_water_out=m_water_out,
        imbalance=imbalance,
    )


def solve_air_flow_rate(
    T_db1: float,
    T_wb1: float,
    T_db2: float,
    T_wb2: float,
    T_w3: float,
    T_w4: float,
    m_water: float,
    P: float = 101_325.0,
) -> float:
    """Back-calculate dry-air mass-flow rate from measured temperatures.

    Rearranges the combined energy and mass balance to solve for ṁₐ:

    .. math::

        \\dot{m}_a = \\frac{\\dot{m}_3 (h_3 - h_4^{*})}{
            (h_2 - h_1) + (\\omega_2 - \\omega_1) h_4^{*}}

    where h₄* = cₚ,w · T_w4 is the enthalpy of the make-up water.

    Args:
        T_db1:   Inlet air dry-bulb temperature  [°C].
        T_wb1:   Inlet air wet-bulb temperature  [°C].
        T_db2:   Outlet air dry-bulb temperature  [°C].
        T_wb2:   Outlet air wet-bulb temperature  [°C].
        T_w3:    Inlet water temperature  [°C].
        T_w4:    Outlet water temperature  [°C].
        m_water: Inlet water mass-flow rate ṁ₃  [kg/s].
        P:       Barometric pressure  [Pa].

    Returns:
        Dry-air mass-flow rate ṁₐ  [kg/s].

    Raises:
        ValueError: If the denominator is near zero (no driving force).
    """
    from cooltower.psychrometrics import humidity_ratio

    omega1 = humidity_ratio(T_db1, T_wb1, P)
    omega2 = humidity_ratio(T_db2, T_wb2, P)
    h1 = specific_enthalpy(T_db1, omega1)
    h2 = specific_enthalpy(T_db2, omega2)
    h3 = CP_WATER * T_w3
    h4 = CP_WATER * T_w4

    numerator = m_water * (h3 - h4)
    denominator = (h2 - h1) + (omega2 - omega1) * h4

    if abs(denominator) < 1e-3:
        raise ValueError(
            "Denominator in air-flow-rate calculation is near zero. "
            "There is no thermodynamic driving force — check that inlet and "
            "outlet conditions differ meaningfully."
        )

    m_air = numerator / denominator
    if m_air <= 0:
        raise ValueError(
            f"Computed ṁₐ = {m_air:.4f} kg/s ≤ 0.  "
            "Ensure T_w3 > T_w4 (water is being cooled) and that air exits "
            "with higher enthalpy than it enters."
        )

    logger.debug("solve_air_flow_rate → ṁₐ = %.4f kg/s", m_air)
    return m_air


def cooling_duty(
    m_water: float,
    T_in: float,
    T_out: float,
) -> float:
    """Compute the water-side heat rejection (cooling duty).

    .. math::

        Q = \\dot{m}_3 \\, c_{p,w} \\, (T_3 - T_4)

    Note: This is an approximation that neglects the enthalpy carried
    away by the evaporated fraction.  Use :func:`solve_energy_balance`
    for the rigorous calculation.

    Args:
        m_water: Water mass-flow rate  [kg/s].
        T_in:    Inlet water temperature  [°C].
        T_out:   Outlet water temperature  [°C].

    Returns:
        Approximate cooling duty  [W].
    """
    if m_water < 0:
        raise ValueError(f"m_water must be non-negative, got {m_water}.")
    return m_water * CP_WATER * (T_in - T_out)


def approach_temperature(T_water_out: float, T_wb_in: float) -> float:
    """Compute the cooling-tower approach temperature.

    The approach is the difference between the cold-water outlet
    temperature and the inlet wet-bulb temperature.  A smaller approach
    indicates better tower performance but requires more transfer units.

    Args:
        T_water_out: Outlet (cold) water temperature  [°C].
        T_wb_in:     Inlet air wet-bulb temperature  [°C].

    Returns:
        Approach temperature  [K or °C difference].

    Raises:
        ValueError: If approach < 0 (thermodynamically impossible in a
            conventional cooling tower).
    """
    approach = T_water_out - T_wb_in
    if approach < -0.5:
        raise ValueError(
            f"Approach temperature = {approach:.2f} °C < 0.  The cold-water "
            "outlet cannot be cooler than the inlet wet-bulb temperature in a "
            "direct-contact cooling tower."
        )
    logger.debug("approach_temperature = %.2f °C", approach)
    return approach


def range_temperature(T_water_in: float, T_water_out: float) -> float:
    """Compute the cooling range (hot → cold water temperature difference).

    Args:
        T_water_in:  Hot water inlet temperature  [°C].
        T_water_out: Cold water outlet temperature  [°C].

    Returns:
        Cooling range  [K or °C difference].
    """
    rng = T_water_in - T_water_out
    logger.debug("range_temperature = %.2f °C", rng)
    return rng


def liquid_to_gas_ratio(m_water: float, m_air: float) -> float:
    """Compute the liquid-to-gas (L/G) mass-flow ratio.

    The L/G ratio determines the operating line on the enthalpy–
    temperature diagram and is a primary design variable for tower
    sizing.  Typical values are 0.75–1.50.

    Args:
        m_water: Water mass-flow rate  [kg/s].
        m_air:   Dry-air mass-flow rate  [kg/s].

    Returns:
        L/G ratio (dimensionless).

    Raises:
        ValueError: If *m_air* ≤ 0.
    """
    if m_air <= 0:
        raise ValueError(f"m_air must be positive, got {m_air}.")
    if m_water < 0:
        raise ValueError(f"m_water must be non-negative, got {m_water}.")
    return m_water / m_air
