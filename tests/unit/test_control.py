"""Unit tests for cooltower.control."""

from __future__ import annotations

import pytest

from cooltower.control import (
    FOPDTModel,
    PIParameters,
    closed_loop_response,
    identify_fopdt,
    performance_indices,
    step_response,
    tune_cohen_coon,
    tune_lambda,
    tune_ziegler_nichols,
)

# ---------------------------------------------------------------------------
# FOPDTModel
# ---------------------------------------------------------------------------


class TestFOPDTModel:
    def test_valid_model(self, fopdt_typical: FOPDTModel) -> None:
        assert fopdt_typical.K_p == pytest.approx(0.75)
        assert fopdt_typical.tau_p == pytest.approx(120.0)
        assert fopdt_typical.theta == pytest.approx(15.0)

    def test_negative_tau_raises(self) -> None:
        with pytest.raises(ValueError, match="tau_p"):
            FOPDTModel(K_p=1.0, tau_p=-1.0, theta=5.0)

    def test_negative_theta_raises(self) -> None:
        with pytest.raises(ValueError, match="theta"):
            FOPDTModel(K_p=1.0, tau_p=60.0, theta=-1.0)

    def test_zero_theta_allowed(self) -> None:
        m = FOPDTModel(K_p=1.0, tau_p=60.0, theta=0.0)
        assert m.theta == 0.0


# ---------------------------------------------------------------------------
# tune_lambda
# ---------------------------------------------------------------------------


class TestTuneLambda:
    def test_kc_positive(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical)
        assert pi.K_c > 0.0

    def test_tau_i_equals_tau_p(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical)
        assert pi.tau_I == pytest.approx(fopdt_typical.tau_p, rel=1e-6)

    def test_method_label(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical)
        assert pi.method == "lambda"

    def test_larger_lambda_smaller_kc(self, fopdt_typical: FOPDTModel) -> None:
        pi_fast = tune_lambda(fopdt_typical, lambda_=20.0)
        pi_slow = tune_lambda(fopdt_typical, lambda_=200.0)
        assert pi_slow.K_c < pi_fast.K_c

    def test_invalid_lambda_raises(self, fopdt_typical: FOPDTModel) -> None:
        with pytest.raises(ValueError, match="lambda_"):
            tune_lambda(fopdt_typical, lambda_=-10.0)

    def test_k_i_property(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical)
        assert pi.K_i == pytest.approx(pi.K_c / pi.tau_I, rel=1e-6)


# ---------------------------------------------------------------------------
# tune_ziegler_nichols
# ---------------------------------------------------------------------------


class TestTuneZieglerNichols:
    def test_returns_pi(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_ziegler_nichols(fopdt_typical)
        assert pi.K_c > 0.0
        assert pi.tau_I == pytest.approx(3.33 * fopdt_typical.theta, rel=1e-3)

    def test_zero_theta_raises(self) -> None:
        model = FOPDTModel(K_p=1.0, tau_p=60.0, theta=0.0)
        with pytest.raises(ValueError, match="negligible dead time"):
            tune_ziegler_nichols(model)

    def test_more_aggressive_than_lambda(self, fopdt_typical: FOPDTModel) -> None:
        """Ziegler-Nichols should give higher gain than conservative lambda."""
        pi_zn = tune_ziegler_nichols(fopdt_typical)
        pi_lam = tune_lambda(fopdt_typical, lambda_=3 * fopdt_typical.theta)
        # ZN is typically more aggressive (higher K_c)
        assert pi_zn.K_c > pi_lam.K_c


# ---------------------------------------------------------------------------
# tune_cohen_coon
# ---------------------------------------------------------------------------


class TestTuneCohenCoon:
    def test_returns_positive_params(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_cohen_coon(fopdt_typical)
        assert pi.K_c > 0.0
        assert pi.tau_I > 0.0

    def test_method_label(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_cohen_coon(fopdt_typical)
        assert pi.method == "cohen_coon"


# ---------------------------------------------------------------------------
# step_response
# ---------------------------------------------------------------------------


class TestStepResponse:
    def test_length_consistent(self, fopdt_typical: FOPDTModel) -> None:
        t, y = step_response(fopdt_typical, t_end=300.0, dt=1.0)
        assert len(t) == len(y)

    def test_initial_output_zero_during_dead_time(self, fopdt_typical: FOPDTModel) -> None:
        t, y = step_response(fopdt_typical, t_end=200.0, dt=1.0)
        # During dead time output should be zero
        dead_steps = int(fopdt_typical.theta)
        for i in range(dead_steps):
            assert y[i] == pytest.approx(0.0)

    def test_asymptotic_value(self, fopdt_typical: FOPDTModel) -> None:
        step_mag = 5.0
        t, y = step_response(fopdt_typical, t_end=1000.0, dt=1.0, step_magnitude=step_mag)
        expected_ss = fopdt_typical.K_p * step_mag
        assert y[-1] == pytest.approx(expected_ss, rel=1e-2)


# ---------------------------------------------------------------------------
# closed_loop_response
# ---------------------------------------------------------------------------


class TestClosedLoopResponse:
    def test_output_approaches_setpoint(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical, lambda_=50.0)
        sp = 5.0
        t, y, u = closed_loop_response(fopdt_typical, pi, setpoint=sp, t_end=1200.0, dt=2.0)
        # Final value should be within 5 % of setpoint
        assert abs(y[-1] - sp) / sp < 0.05

    def test_lengths_consistent(self, fopdt_typical: FOPDTModel) -> None:
        pi = tune_lambda(fopdt_typical)
        t, y, u = closed_loop_response(fopdt_typical, pi, setpoint=1.0, t_end=300.0)
        assert len(t) == len(y) == len(u)


# ---------------------------------------------------------------------------
# identify_fopdt
# ---------------------------------------------------------------------------


class TestIdentifyFOPDT:
    def test_two_point_recovers_approximate_model(self, step_test_data: dict) -> None:
        model = identify_fopdt(
            time=step_test_data["time"],
            output=step_test_data["output"],
            step_time=step_test_data["step_time"],
            step_magnitude=step_test_data["step_magnitude"],
            method="two_point",
        )
        true = step_test_data["model"]
        # K_p within 10 %
        assert abs(model.K_p - true.K_p) / true.K_p < 0.10
        # tau_p within 30 % (due to noise)
        assert abs(model.tau_p - true.tau_p) / true.tau_p < 0.30

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="same length"):
            identify_fopdt([0, 1, 2], [0, 1], 0.0, 1.0)

    def test_zero_step_magnitude_raises(self) -> None:
        with pytest.raises(ValueError, match="step_magnitude"):
            identify_fopdt([0, 1, 2], [0, 0.1, 0.2], 0.0, 0.0)

    def test_unknown_method_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown method"):
            identify_fopdt([0, 1, 2, 3], [0, 0, 0.5, 1.0], 1.0, 1.0, method="banana")


# ---------------------------------------------------------------------------
# performance_indices
# ---------------------------------------------------------------------------


class TestPerformanceIndices:
    def test_zero_error(self) -> None:
        t = [0.0, 1.0, 2.0]
        e = [0.0, 0.0, 0.0]
        idx = performance_indices(t, e)
        assert idx["ISE"] == pytest.approx(0.0)
        assert idx["IAE"] == pytest.approx(0.0)
        assert idx["ITAE"] == pytest.approx(0.0)

    def test_ise_iae_itae_ordering(self) -> None:
        t = [float(i) for i in range(100)]
        e = [1.0] * 100
        idx = performance_indices(t, e)
        # ITAE should be much larger than IAE for long time horizons
        assert idx["ITAE"] > idx["IAE"]

    def test_length_mismatch_raises(self) -> None:
        with pytest.raises(ValueError, match="equal length"):
            performance_indices([0, 1, 2], [0, 1])
