"""Physical and thermodynamic constants for cooling tower analysis.

All values are in SI units unless explicitly stated otherwise.
"""

from __future__ import annotations

__all__ = [
    "CP_WATER",
    "CP_AIR",
    "CP_VAPOUR",
    "H_FG_0",
    "R_AIR",
    "R_VAPOUR",
    "P_STD",
    "T_REF",
    "EPSILON",
    "MR_RATIO",
]

# ── Specific heat capacities [J/(kg·K)] ──────────────────────────────────────
CP_WATER: float = 4186.0
"""Specific heat of liquid water at ~25 °C  [J/(kg·K)]."""

CP_AIR: float = 1006.0
"""Specific heat of dry air at constant pressure  [J/(kg·K)]."""

CP_VAPOUR: float = 1805.0
"""Specific heat of water vapour at constant pressure  [J/(kg·K)]."""

# ── Latent heat ───────────────────────────────────────────────────────────────
H_FG_0: float = 2_501_000.0
"""Latent heat of vaporisation of water at 0 °C  [J/kg]."""

# ── Gas constants [J/(kg·K)] ──────────────────────────────────────────────────
R_AIR: float = 287.058
"""Specific gas constant of dry air  [J/(kg·K)]."""

R_VAPOUR: float = 461.522
"""Specific gas constant of water vapour  [J/(kg·K)]."""

# ── Reference conditions ──────────────────────────────────────────────────────
P_STD: float = 101_325.0
"""Standard atmospheric pressure  [Pa]."""

T_REF: float = 273.15
"""Reference temperature offset: 0 °C in Kelvin  [K]."""

# ── Molecular mass ratio ──────────────────────────────────────────────────────
MR_RATIO: float = 0.621_945
"""Ratio of molecular mass of water vapour to dry air (18.015 / 28.964)."""

# ── Numerical tolerance ───────────────────────────────────────────────────────
EPSILON: float = 1e-9
"""Small positive number used to guard against division-by-zero."""
