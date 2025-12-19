"""PI/PID controller tuning and step-response analysis for cooling tower control.

Implements:
  - First-order-plus-dead-time (FOPDT) model identification from step tests.
  - Lambda (IMC-based) tuning for PI controllers.
  - Ziegler–Nichols and Cohen–Coon rules (for comparison).
  - Step-response simulation for closed-loop performance assessment.

Lambda tuning is preferred for cooling tower temperature control because:
  1. The process has slow, noisy dynamics (temperature lags flow changes).
  2. The closed-loop time constant λ can be tuned conservatively to avoid
     windup in the presence of flow-rate noise.
  3. Unlike Ziegler–Nichols, no sustained oscillations are needed during
     commissioning.

Reference: Seborg, Edgar, Mellichamp & Doyle, *Process Dynamics and Control*,
4th ed., Chapters 11–12.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Sequence

__all__ = [
    "FOPDTModel",
    "PIParameters",
    "identify_fopdt",
    "tune_lambda",
    "tune_ziegler_nichols",
    "tune_cohen_coon",
    "step_response",
    "closed_loop_response",
    "performance_indices",
]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FOPDTModel:
    """First-order-plus-dead-time (FOPDT) process model.

    .. math::

        G_p(s) = \\frac{K_p e^{-\\theta s}}{\\tau_p s + 1}

    Attributes:
        K_p:   Process gain  [output_units / input_units].
        tau_p: Process time constant  [s].
        theta: Dead time (transport delay)  [s].
    """

    K_p: float
    tau_p: float
    theta: float

    def __post_init__(self) -> None:
        if self.tau_p <= 0:
            raise ValueError(f"Process time constant tau_p must be positive, got {self.tau_p}.")
        if self.theta < 0:
            raise ValueError(f"Dead time theta must be non-negative, got {self.theta}.")


@dataclass(frozen=True)
class PIParameters:
    """Tuned PI controller parameters.

    In parallel (ISA) form:

    .. math::

        u(t) = K_c \\left[ e(t) + \\frac{1}{\\tau_I} \\int e \\, dt \\right]

    Attributes:
        K_c:   Controller gain.
        tau_I: Integral time constant  [s].
        method: Tuning method used (e.g., ``"lambda"``, ``"ziegler_nichols"``).
    """

    K_c: float
    tau_I: float
    method: str

    @property
    def K_i(self) -> float:
        """Integral gain K_I = K_c / τ_I."""
        return self.K_c / self.tau_I

    def __str__(self) -> str:
        return (
            f"PI({self.method}): K_c={self.K_c:.4f},  τ_I={self.tau_I:.2f} s  "
            f"(K_I={self.K_i:.4f})"
        )


# ---------------------------------------------------------------------------
# Model identification
# ---------------------------------------------------------------------------


def identify_fopdt(
    time: Sequence[float],
    output: Sequence[float],
    step_time: float,
    step_magnitude: float,
    method: str = "tangent",
) -> FOPDTModel:
    """Identify FOPDT parameters from an open-loop step-test record.

    Two methods are supported:

    * ``"tangent"`` – draws the tangent at the inflection point of the
      step response (classical graphical method; sensitive to noise).
    * ``"two_point"`` – uses the 28 % and 63 % response times to
      estimate τ_p and θ (more robust for noisy field data).

    Args:
        time:            Time vector  [s], monotonically increasing.
        output:          Measured output (e.g., temperature) vector.
        step_time:       Time at which the step was applied  [s].
        step_magnitude:  Magnitude of the step change in the manipulated
                         variable (e.g., % valve opening).
        method:          Identification method: ``"tangent"`` or
                         ``"two_point"``.

    Returns:
        An :class:`FOPDTModel` with identified K_p, τ_p, θ.

    Raises:
        ValueError: If *time* and *output* have different lengths, if
            *step_magnitude* is zero, or if the response does not reach
            28 % of its final value.
    """
    t = list(time)
    y = list(output)

    if len(t) != len(y):
        raise ValueError(
            f"time and output must have the same length, got {len(t)} and {len(y)}."
        )
    if abs(step_magnitude) < 1e-9:
        raise ValueError("step_magnitude must be non-zero.")

    # Trim to post-step data
    pre_step = [yi for ti, yi in zip(t, y) if ti < step_time]
    post_step_t = [ti for ti in t if ti >= step_time]
    post_step_y = [yi for ti, yi in zip(t, y) if ti >= step_time]

    if not post_step_y:
        raise ValueError(
            "step_time must lie within the time vector such that post-step data exist."
        )
    # If step applied at t[0] there is no pre-step baseline; use first sample
    if not pre_step:
        pre_step = [post_step_y[0]]

    y0 = sum(pre_step[-min(5, len(pre_step)) :]) / min(5, len(pre_step))
    y_inf = post_step_y[-1]
    delta_y = y_inf - y0
    K_p = delta_y / step_magnitude

    if method == "two_point":
        t28 = _interpolate_response_time(post_step_t, post_step_y, y0, 0.283 * delta_y)
        t63 = _interpolate_response_time(post_step_t, post_step_y, y0, 0.632 * delta_y)
        tau_p = 1.5 * (t63 - t28)
        theta = t63 - tau_p - (post_step_t[0] - step_time)
    elif method == "tangent":
        # Find inflection point (max slope in post-step data)
        slopes = [
            (post_step_y[i + 1] - post_step_y[i]) / max(post_step_t[i + 1] - post_step_t[i], 1e-9)
            for i in range(len(post_step_y) - 1)
        ]
        idx_max = max(range(len(slopes)), key=lambda i: abs(slopes[i]))
        slope_max = slopes[idx_max]
        t_infl = post_step_t[idx_max]
        y_infl = post_step_y[idx_max]
        # Tangent line: y = y_infl + slope_max * (t - t_infl)
        t_start = t_infl - (y_infl - y0) / (slope_max + 1e-12)
        t_end = t_infl + (y_inf - y_infl) / (slope_max + 1e-12)
        theta = max(t_start - step_time, 0.0)
        tau_p = t_end - t_start
    else:
        raise ValueError(f"Unknown method '{method}'. Choose 'tangent' or 'two_point'.")

    tau_p = max(tau_p, 1e-3)
    theta = max(theta, 0.0)

    model = FOPDTModel(K_p=K_p, tau_p=tau_p, theta=theta)
    logger.info(
        "FOPDT identified (%s): K_p=%.4f, τ_p=%.2f s, θ=%.2f s",
        method,
        K_p,
        tau_p,
        theta,
    )
    return model


def _interpolate_response_time(
    t: list[float],
    y: list[float],
    y0: float,
    target_delta: float,
) -> float:
    """Find the time at which (y − y0) first reaches *target_delta* by linear interpolation."""
    for i in range(len(y) - 1):
        if (y[i] - y0) <= target_delta <= (y[i + 1] - y0):
            frac = (target_delta - (y[i] - y0)) / max(y[i + 1] - y[i], 1e-12)
            return t[i] + frac * (t[i + 1] - t[i])
    raise ValueError(
        f"Response does not reach the target level (Δy = {target_delta:.4f}) "
        "within the supplied data range.  Extend the step test or check inputs."
    )


# ---------------------------------------------------------------------------
# Tuning rules
# ---------------------------------------------------------------------------


def tune_lambda(
    model: FOPDTModel,
    lambda_: float | None = None,
    lambda_factor: float = 3.0,
) -> PIParameters:
    """Compute PI parameters using the IMC-based lambda (closed-loop) tuning rule.

    .. math::

        K_c = \\frac{\\tau_p}{K_p (\\lambda + \\theta)}, \\quad
        \\tau_I = \\tau_p

    *λ* is the desired closed-loop time constant.  Larger λ → slower,
    more robust response.  The default heuristic is λ = 3θ (Seborg et al.
    recommend λ ≥ θ).

    Lambda tuning is preferred over Ziegler–Nichols for cooling tower
    temperature control because:
      - The process is slow-moving; aggressive integral action causes
        windup against the valve envelope.
      - Flow-rate measurements are noisy; Ziegler–Nichols would require
        finding the ultimate gain under steady oscillation, which is
        destructive to operation.
      - Cohen–Coon needs a large dead-time fraction (θ/τ_p ≥ 0.1),
        not typical of temperature loops.

    Args:
        model:         Identified FOPDT model.
        lambda_:       Desired closed-loop time constant  [s].  If
                       ``None``, defaults to ``lambda_factor * theta``
                       (or ``tau_p / 2`` if theta ≈ 0).
        lambda_factor: Multiplier applied to θ when *lambda_* is
                       ``None``.  Ignored if *lambda_* is supplied.

    Returns:
        Tuned :class:`PIParameters`.

    Raises:
        ValueError: If the resulting K_c or τ_I are non-physical.
    """
    if lambda_ is None:
        lambda_ = max(lambda_factor * model.theta, model.tau_p / 2.0)
        logger.debug("Lambda auto-set to %.2f s (= max(%.1f·θ, τ_p/2))", lambda_, lambda_factor)

    if lambda_ <= 0:
        raise ValueError(f"lambda_ must be positive, got {lambda_}.")

    K_c = model.tau_p / (model.K_p * (lambda_ + model.theta))
    tau_I = model.tau_p

    params = PIParameters(K_c=K_c, tau_I=tau_I, method="lambda")
    logger.info("Lambda tuning → %s", params)
    return params


def tune_ziegler_nichols(model: FOPDTModel) -> PIParameters:
    """Compute PI parameters using the Ziegler–Nichols open-loop (reaction-curve) rule.

    .. math::

        K_c = \\frac{0.9 \\tau_p}{K_p \\theta}, \\quad
        \\tau_I = 3.33 \\theta

    .. warning::
        Ziegler–Nichols produces aggressive settings with ~25 % overshoot
        and is **not recommended** for cooling tower temperature control
        (noisy, slow dynamics).  Provided for comparison only.

    Args:
        model: Identified FOPDT model.

    Returns:
        Tuned :class:`PIParameters`.

    Raises:
        ValueError: If dead time θ ≈ 0 (formula is singular).
    """
    if model.theta < 1e-6:
        raise ValueError(
            "Ziegler–Nichols open-loop tuning requires θ > 0. "
            "The process has negligible dead time — use lambda tuning instead."
        )

    K_c = 0.9 * model.tau_p / (model.K_p * model.theta)
    tau_I = 3.33 * model.theta
    params = PIParameters(K_c=K_c, tau_I=tau_I, method="ziegler_nichols")
    logger.info("Ziegler–Nichols tuning → %s", params)
    return params


def tune_cohen_coon(model: FOPDTModel) -> PIParameters:
    """Compute PI parameters using the Cohen–Coon rule.

    .. math::

        K_c = \\frac{\\tau_p}{K_p \\theta}
              \\left(0.9 + \\frac{\\theta}{12 \\tau_p}\\right), \\quad
        \\tau_I = \\theta \\frac{30 + 3(\\theta/\\tau_p)}{9 + 20(\\theta/\\tau_p)}

    .. warning::
        Cohen–Coon is designed for processes where θ/τ_p ∈ [0.1, 1.0].
        It requires quick output responses to step changes, making it
        unsuitable for slow temperature dynamics in cooling towers.

    Args:
        model: Identified FOPDT model.

    Returns:
        Tuned :class:`PIParameters`.
    """
    r = model.theta / model.tau_p
    K_c = (model.tau_p / (model.K_p * model.theta)) * (0.9 + r / 12.0)
    tau_I = model.theta * (30.0 + 3.0 * r) / (9.0 + 20.0 * r)
    params = PIParameters(K_c=K_c, tau_I=tau_I, method="cohen_coon")
    logger.info("Cohen–Coon tuning → %s", params)
    return params


# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------


def step_response(
    model: FOPDTModel,
    t_end: float,
    dt: float = 1.0,
    step_magnitude: float = 1.0,
) -> tuple[list[float], list[float]]:
    """Simulate the open-loop step response of an FOPDT model.

    Args:
        model:           FOPDT process model.
        t_end:           Simulation end time  [s].
        dt:              Time step  [s].
        step_magnitude:  Input step magnitude.

    Returns:
        Tuple ``(time, output)`` as plain Python lists.
    """
    n = int(t_end / dt) + 1
    time = [i * dt for i in range(n)]
    output = []
    for t in time:
        t_eff = t - model.theta
        if t_eff <= 0:
            output.append(0.0)
        else:
            output.append(model.K_p * step_magnitude * (1.0 - math.exp(-t_eff / model.tau_p)))
    return time, output


def closed_loop_response(
    model: FOPDTModel,
    pi: PIParameters,
    setpoint: float,
    t_end: float,
    dt: float = 1.0,
    disturbance: float = 0.0,
    disturbance_time: float | None = None,
) -> tuple[list[float], list[float], list[float]]:
    """Simulate the closed-loop PI step response using Euler integration.

    Implements the velocity (incremental) form to avoid reset windup:

    .. math::

        \\Delta u_k = K_c \\left[ (e_k - e_{k-1})
            + \\frac{\\Delta t}{\\tau_I} e_k \\right]

    Args:
        model:            FOPDT process model.
        pi:               Tuned PI parameters.
        setpoint:         Step setpoint change magnitude.
        t_end:            Simulation end time  [s].
        dt:               Integration step  [s].
        disturbance:      Magnitude of a load disturbance applied at
                          *disturbance_time*  [process output units].
        disturbance_time: Time of load disturbance  [s].  Ignored if
                          *disturbance* == 0.

    Returns:
        Tuple ``(time, output, control_signal)``.
    """
    n = int(t_end / dt) + 1
    time_vec = [i * dt for i in range(n)]

    # Store delayed outputs for dead-time approximation (integer delay)
    delay_steps = max(1, int(model.theta / dt))
    u_history: list[float] = [0.0] * (delay_steps + n)

    y = 0.0
    u = 0.0
    e_prev = 0.0

    out_y: list[float] = []
    out_u: list[float] = []

    for k, t in enumerate(time_vec):
        sp = setpoint if t >= 0 else 0.0
        dist = (
            disturbance
            if (disturbance_time is not None and t >= disturbance_time)
            else 0.0
        )
        e = sp - y - dist
        delta_u = pi.K_c * ((e - e_prev) + dt / pi.tau_I * e)
        u = u + delta_u
        e_prev = e

        u_history[k + delay_steps] = u
        u_delayed = u_history[k]

        # First-order Euler update
        dy = (-y + model.K_p * u_delayed) / model.tau_p * dt
        y += dy

        out_y.append(y)
        out_u.append(u)

    return time_vec, out_y, out_u


def performance_indices(
    time: Sequence[float],
    error: Sequence[float],
) -> dict[str, float]:
    """Compute integral performance indices from a closed-loop error sequence.

    Computes:
      - **ISE** – Integral of Squared Error  ∫ e² dt
      - **IAE** – Integral of Absolute Error  ∫ |e| dt
      - **ITAE** – Integral of Time-weighted Absolute Error  ∫ t|e| dt

    Args:
        time:  Time vector  [s].
        error: Error (setpoint − output) vector.

    Returns:
        Dictionary with keys ``"ISE"``, ``"IAE"``, ``"ITAE"``.

    Raises:
        ValueError: If *time* and *error* have different lengths.
    """
    t = list(time)
    e = list(error)
    if len(t) != len(e):
        raise ValueError(
            f"time and error must have equal length, got {len(t)} and {len(e)}."
        )

    ISE = IAE = ITAE = 0.0
    for i in range(len(t) - 1):
        dt = t[i + 1] - t[i]
        ISE += 0.5 * (e[i] ** 2 + e[i + 1] ** 2) * dt
        IAE += 0.5 * (abs(e[i]) + abs(e[i + 1])) * dt
        ITAE += 0.5 * (t[i] * abs(e[i]) + t[i + 1] * abs(e[i + 1])) * dt

    indices = {"ISE": ISE, "IAE": IAE, "ITAE": ITAE}
    logger.debug("Performance indices: %s", indices)
    return indices
