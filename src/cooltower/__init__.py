"""cooltower — Cooling Tower Thermodynamic & Control Analysis Library.

A pvlib-inspired Python package for engineering analysis of mechanical-
draught cooling towers.  Covers:

- **Psychrometrics**: humidity ratio, specific enthalpy, wet-bulb
  temperature, saturation pressure (Rogers & Mayhew / ASHRAE).
- **Energy balances**: steady-flow mass and energy balances, evaporation
  rates, approach and range temperatures, L/G ratio.
- **Process control**: FOPDT model identification from step tests; lambda
  (IMC), Ziegler–Nichols, and Cohen–Coon PI tuning; closed-loop
  simulation and performance indices.

Quickstart
----------
>>> from cooltower.psychrometrics import humidity_ratio, specific_enthalpy
>>> omega = humidity_ratio(T_db=28.0, T_wb=22.0)
>>> h = specific_enthalpy(T_db=28.0, omega=omega)
>>> print(f"ω = {omega:.4f} kg/kg,  h = {h/1000:.2f} kJ/kg_da")

>>> from cooltower.control import FOPDTModel, tune_lambda
>>> model = FOPDTModel(K_p=0.8, tau_p=120.0, theta=15.0)
>>> pi = tune_lambda(model)
>>> print(pi)
"""

from __future__ import annotations

__version__ = "0.1.0"
__author__ = "Defne Nihal Ertugrul"
__license__ = "MIT"

__all__ = [
    "__version__",
    "__author__",
    "__license__",
    "psychrometrics",
    "energy_balance",
    "control",
    "constants",
]

from cooltower import constants, control, energy_balance, psychrometrics
