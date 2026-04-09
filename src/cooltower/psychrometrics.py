"""Psychrometric calculations for moist-air / cooling-tower analysis.

All functions operate on SI units.  Temperature inputs are in degrees
Celsius for readability; internal conversions to Kelvin are explicit.

Equations follow Rogers & Mayhew *Engineering Thermodynamics: Work and
Heat Transfer* (4th ed.) and ASHRAE Fundamentals (2021).

Example
-------
>>> from cooltower.psychrometrics import humidity_ratio, specific_enthalpy
>>> omega = humidity_ratio(T_db=30.0, T_wb=25.0)
>>> h = specific_enthalpy(T_db=30.0, omega=omega)
"""

from __future__ import annotations

import logging
import math

from cooltower.constants import (
    CP_AIR,
    CP_VAPOUR,
    CP_WATER,
    EPSILON,
    H_FG_0,
    MR_RATIO,
    P_STD,
    T_REF,
)

__all__ = [
    "saturation_pressure",
    "humidity_ratio_from_rh",
    "humidity_ratio",
    "specific_enthalpy",
    "wet_bulb_temperature",
    "dew_point_temperature",
    "relative_humidity",
    "degree_of_saturation",
    "evaporation_rate",
]

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _celsius_to_kelvin(T_C: float) -> float:
    """Convert Celsius to Kelvin."""
    return T_C + T_REF


def _validate_temperature(T_C: float, name: str = "T") -> None:
    """Raise ValueError if temperature is physically implausible."""
    if T_C < -100.0 or T_C > 200.0:
        raise ValueError(
            f"{name} = {T_C} °C is outside the physically valid range "
            "[-100, 200] °C for moist-air psychrometrics."
        )


def _validate_pressure(P: float) -> None:
    """Raise ValueError if pressure is non-positive or unrealistically large."""
    if P <= 0:
        raise ValueError(f"Pressure must be positive, got P = {P} Pa.")
    if P > 2e6:
        raise ValueError(
            f"Pressure P = {P} Pa is unrealistically large for atmospheric "
            "psychrometrics (max sensible value ≈ 2 MPa)."
        )


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def saturation_pressure(T_db: float) -> float:
    """Compute the saturation vapour pressure using the Antoine / Magnus form.

    Uses the Tetens (1930) approximation, accurate to ±0.1 % over
    −40 °C to 60 °C.  The Buck (1981) coefficients are applied for
    temperatures above 0 °C and the ice-surface form below 0 °C.

    Args:
        T_db: Dry-bulb temperature  [°C].

    Returns:
        Saturation vapour pressure  [Pa].

    Raises:
        ValueError: If *T_db* is outside [-100, 200] °C.

    Example:
        >>> round(saturation_pressure(25.0), 1)
        3169.9
    """
    _validate_temperature(T_db, "T_db")
    if T_db >= 0.0:
        # Buck (1981) – liquid water
        p_sat = 611.21 * math.exp((18.678 - T_db / 234.5) * (T_db / (257.14 + T_db)))
    else:
        # Buck (1981) – ice surface
        p_sat = 611.15 * math.exp((23.036 - T_db / 333.7) * (T_db / (279.82 + T_db)))
    logger.debug("saturation_pressure(T_db=%.2f °C) = %.4f Pa", T_db, p_sat)
    return p_sat


def humidity_ratio_from_rh(
    T_db: float,
    rh: float,
    P: float = P_STD,
) -> float:
    """Compute humidity ratio from dry-bulb temperature and relative humidity.

    The humidity ratio (specific humidity) is defined as the mass of
    water vapour per unit mass of dry air:

    .. math::
        \\omega = 0.621945 \\, \\frac{\\phi \\, p_{sat}}{P - \\phi \\, p_{sat}}

    Args:
        T_db: Dry-bulb temperature  [°C].
        rh:   Relative humidity as a fraction in [0, 1].
        P:    Total (barometric) pressure  [Pa].  Defaults to standard
              atmosphere (101 325 Pa).

    Returns:
        Humidity ratio ω  [kg_water / kg_dryair].

    Raises:
        ValueError: If *rh* is outside [0, 1] or if vapour pressure
            exceeds total pressure.
    """
    _validate_temperature(T_db, "T_db")
    _validate_pressure(P)
    if not 0.0 <= rh <= 1.0:
        raise ValueError(f"Relative humidity must be in [0, 1], got rh = {rh}.")

    p_sat = saturation_pressure(T_db)
    p_v = rh * p_sat
    if p_v >= P:
        raise ValueError(
            f"Vapour pressure p_v = {p_v:.2f} Pa ≥ total pressure P = {P:.2f} Pa. "
            "Check that T_db and rh are consistent with the given pressure."
        )
    omega = MR_RATIO * p_v / (P - p_v)
    logger.debug("humidity_ratio_from_rh → ω = %.6f kg/kg", omega)
    return omega


def humidity_ratio(
    T_db: float,
    T_wb: float,
    P: float = P_STD,
) -> float:
    """Compute humidity ratio from dry- and wet-bulb temperatures.

    Applies the Sprung psychrometric formula (sling psychrometer):

    .. math::
        \\omega = \\omega_{sat}(T_{wb}) -
        A \\cdot (T_{db} - T_{wb}) \\cdot P / 1000

    where *A* = 6.6×10⁻⁴ K⁻¹ for a ventilated (aspirated) psychrometer.

    Args:
        T_db: Dry-bulb temperature  [°C].
        T_wb: Wet-bulb temperature  [°C].
        P:    Barometric pressure  [Pa].

    Returns:
        Humidity ratio ω  [kg_water / kg_dryair].

    Raises:
        ValueError: If *T_wb* > *T_db* or if temperatures are out of range.
    """
    _validate_temperature(T_db, "T_db")
    _validate_temperature(T_wb, "T_wb")
    _validate_pressure(P)
    if T_wb > T_db + EPSILON:
        raise ValueError(
            f"Wet-bulb temperature T_wb = {T_wb} °C cannot exceed "
            f"dry-bulb temperature T_db = {T_db} °C."
        )

    # Sprung formula in vapour-pressure space:
    #   p_v = p_sat(T_wb) - A * P * (T_db - T_wb)
    # where A = 6.6×10⁻⁴ K⁻¹ for a ventilated (Assmann) psychrometer.
    A = 6.6e-4
    p_sat_wb = saturation_pressure(T_wb)
    p_v = p_sat_wb - A * P * (T_db - T_wb)
    p_v = max(p_v, 0.0)
    omega = MR_RATIO * p_v / max(P - p_v, 1.0)
    logger.debug(
        "humidity_ratio(T_db=%.2f, T_wb=%.2f) → ω = %.6f kg/kg",
        T_db,
        T_wb,
        omega,
    )
    return omega


def specific_enthalpy(T_db: float, omega: float) -> float:
    """Compute the specific enthalpy of moist air.

    Evaluated relative to liquid water at 0 °C and dry air at 0 °C:

    .. math::
        h = c_{p,a} T + \\omega (h_{fg,0} + c_{p,v} T)

    where T is in °C (consistent with the chosen datum).

    Args:
        T_db: Dry-bulb temperature  [°C].
        omega: Humidity ratio  [kg_water / kg_dryair].

    Returns:
        Specific enthalpy  [J / kg_dryair].

    Raises:
        ValueError: If *omega* is negative.
    """
    _validate_temperature(T_db, "T_db")
    if omega < 0.0:
        raise ValueError(f"Humidity ratio cannot be negative, got ω = {omega}.")

    h = CP_AIR * T_db + omega * (H_FG_0 + CP_VAPOUR * T_db)
    logger.debug("specific_enthalpy(T=%.2f, ω=%.6f) → h = %.2f J/kg", T_db, omega, h)
    return h


def relative_humidity(T_db: float, omega: float, P: float = P_STD) -> float:
    """Compute relative humidity from dry-bulb temperature and humidity ratio.

    Args:
        T_db:  Dry-bulb temperature  [°C].
        omega: Humidity ratio  [kg_water / kg_dryair].
        P:     Barometric pressure  [Pa].

    Returns:
        Relative humidity φ in [0, 1].

    Raises:
        ValueError: If *omega* is negative or resulting φ > 1.
    """
    _validate_temperature(T_db, "T_db")
    _validate_pressure(P)
    if omega < 0.0:
        raise ValueError(f"Humidity ratio cannot be negative, got ω = {omega}.")

    p_sat = saturation_pressure(T_db)
    p_v = omega * P / (MR_RATIO + omega)
    rh = p_v / (p_sat + EPSILON)
    if rh > 1.0 + 1e-4:
        raise ValueError(
            f"Computed φ = {rh:.4f} > 1.  The combination T_db = {T_db} °C, "
            f"ω = {omega:.6f} kg/kg, P = {P:.0f} Pa implies supersaturation."
        )
    return min(rh, 1.0)


def dew_point_temperature(omega: float, P: float = P_STD) -> float:
    """Compute dew-point temperature from humidity ratio.

    Inverts the Buck formula for liquid water via Newton–Raphson.

    Args:
        omega: Humidity ratio  [kg_water / kg_dryair].
        P:     Barometric pressure  [Pa].

    Returns:
        Dew-point temperature  [°C].

    Raises:
        ValueError: If *omega* is non-positive.
    """
    _validate_pressure(P)
    if omega <= 0.0:
        raise ValueError(
            f"Humidity ratio must be positive to compute dew point, got ω = {omega}."
        )

    p_v = omega * P / (MR_RATIO + omega)
    # Magnus inversion (Lawrence 2005) – accurate to ±0.35 °C over 0–60 °C
    ln_pv = math.log(p_v / 611.21)
    T_dp = (234.5 * ln_pv) / (17.368 - ln_pv)
    logger.debug("dew_point_temperature(ω=%.6f) → T_dp = %.2f °C", omega, T_dp)
    return T_dp


def wet_bulb_temperature(
    T_db: float,
    omega: float,
    P: float = P_STD,
    tol: float = 1e-4,
    max_iter: int = 50,
) -> float:
    """Estimate wet-bulb temperature via iterative inversion of the Sprung formula.

    Args:
        T_db:     Dry-bulb temperature  [°C].
        omega:    Humidity ratio  [kg_water / kg_dryair].
        P:        Barometric pressure  [Pa].
        tol:      Convergence tolerance on T_wb  [°C].
        max_iter: Maximum number of Newton iterations.

    Returns:
        Wet-bulb temperature  [°C].

    Raises:
        ValueError: If inputs are out of physical range.
        RuntimeError: If the iteration does not converge within *max_iter*.
    """
    _validate_temperature(T_db, "T_db")
    _validate_pressure(P)
    if omega < 0.0:
        raise ValueError(
            f"Humidity ratio cannot be negative, got ω = {omega} kg/kg."
        )

    # Feasibility: ω must not exceed saturation at T_db, otherwise the
    # Newton iteration will wander outside the valid psychrometric region
    # and eventually raise a generic non-convergence error.
    omega_sat = humidity_ratio_from_rh(T_db, 1.0, P)
    if omega > omega_sat * (1.0 + 1e-6):
        raise ValueError(
            f"Humidity ratio ω = {omega:.6f} kg/kg exceeds saturation "
            f"ω_sat = {omega_sat:.6f} kg/kg at T_db = {T_db} °C, "
            f"P = {P:.0f} Pa. The state is supersaturated; no wet-bulb "
            "temperature exists."
        )

    T_wb = T_db  # initial guess
    A = 6.6e-4
    for i in range(max_iter):
        omega_calc = humidity_ratio(T_db, T_wb, P)
        residual = omega_calc - omega
        if abs(residual) < tol:
            logger.debug(
                "wet_bulb_temperature converged in %d iterations → T_wb = %.4f °C",
                i + 1,
                T_wb,
            )
            return T_wb
        # Derivative of omega w.r.t. T_wb (finite difference)
        # Finite-difference Jacobian (robust across all operating points)
        eps = 0.1
        T_wb_hi = min(T_wb + eps, T_db)
        T_wb_lo = max(T_wb - eps, -90.0)
        dOmega_dTwb = (
            humidity_ratio(T_db, T_wb_hi, P) - humidity_ratio(T_db, T_wb_lo, P)
        ) / max(T_wb_hi - T_wb_lo, EPSILON)
        T_wb -= residual / (dOmega_dTwb + EPSILON)

    raise RuntimeError(
        f"wet_bulb_temperature did not converge within {max_iter} iterations "
        f"for T_db={T_db} °C, ω={omega:.6f} kg/kg.  "
        "Check that inputs are physically consistent."
    )


def degree_of_saturation(T_db: float, omega: float, P: float = P_STD) -> float:
    """Compute the degree of saturation μ = ω / ω_sat(T_db).

    Args:
        T_db:  Dry-bulb temperature  [°C].
        omega: Humidity ratio  [kg_water / kg_dryair].
        P:     Barometric pressure  [Pa].

    Returns:
        Degree of saturation μ in [0, 1].
    """
    omega_sat = humidity_ratio_from_rh(T_db, 1.0, P)
    return min(omega / (omega_sat + EPSILON), 1.0)


def evaporation_rate(
    omega_out: float,
    omega_in: float,
    m_air: float,
) -> float:
    """Compute the evaporation (make-up) water mass-flow rate.

    Derived from the steady-flow mass balance on water vapour:

    .. math::
        \\dot{m}_{evap} = (\\omega_2 - \\omega_1) \\, \\dot{m}_a

    Args:
        omega_out: Exit air humidity ratio ω₂  [kg_water / kg_dryair].
        omega_in:  Inlet air humidity ratio ω₁  [kg_water / kg_dryair].
        m_air:     Dry air mass-flow rate ṁₐ  [kg/s].

    Returns:
        Evaporation rate ṁ_evap  [kg/s].  Returns 0.0 if ω_out < ω_in
        (condensation regime — not physically expected in a cooling tower).

    Raises:
        ValueError: If any argument is negative.
    """
    for name, val in [("omega_out", omega_out), ("omega_in", omega_in), ("m_air", m_air)]:
        if val < 0.0:
            raise ValueError(f"{name} must be non-negative, got {val}.")

    m_evap = max(0.0, (omega_out - omega_in) * m_air)
    logger.debug(
        "evaporation_rate(Δω=%.6f, ṁa=%.4f) → ṁ_evap = %.6f kg/s",
        omega_out - omega_in,
        m_air,
        m_evap,
    )
    return m_evap
