"""
Benchmark for the cooltower psychrometrics hot path.

The heaviest operation in cooling-tower analysis is wet_bulb_temperature,
which runs a Newton iteration with a finite-difference Jacobian — each
iteration calls humidity_ratio() three times, and each of those calls
saturation_pressure() plus a full validation stack. A typical design sweep
(hundreds of ambient conditions) spends most of its time here.

# BASELINE (pre-optimization):
#   mean wall time over 5 runs (10,000 wet_bulb_temperature calls): 0.0858 s
#   cProfile (cumulative):
#     wet_bulb_temperature         0.294 s / 10,000 calls
#     humidity_ratio               0.249 s / 70,000 calls
#     logging.Logger.debug         0.118 s / 150,000 calls  ← 40% of runtime
#     saturation_pressure          0.098 s / 70,000 calls
#     logging.isEnabledFor         0.093 s / 150,000 calls
#     _validate_temperature        0.017 s / 220,000 calls
#     math.exp                     0.005 s / 70,000 calls
#
# AFTER optimization (remove logger.debug from 5 hot-path functions):
#   mean wall time over 5 runs: 0.0363 s
#   speedup: 2.36x
#   cProfile (cumulative):
#     wet_bulb_temperature         0.164 s / 10,000 calls
#     humidity_ratio               0.117 s / 70,000 calls
#     saturation_pressure          0.034 s / 70,000 calls
#     logging.Logger.debug         0.010 s / 10,000 calls  (only in wet_bulb now)
#
# The logger.debug calls dispatched through Python-level isEnabledFor even
# when debug was disabled (~0.8 us overhead per call). At 150k calls, that
# was 40% of the sweep's runtime despite producing zero actual log output.
# The remaining wet_bulb_temperature convergence log (1 per outer call) is
# kept — it's useful, low-frequency diagnostic output.
"""

from __future__ import annotations

import cProfile
import pstats
import time
from io import StringIO

from cooltower.psychrometrics import wet_bulb_temperature

N_RUNS = 5
N_POINTS = 10_000


def workload() -> None:
    # Sweep of 10k (T_db, omega) ambient conditions across a realistic envelope.
    for i in range(N_POINTS):
        T_db = 15.0 + 0.002 * i         # 15 ... 35 °C
        omega = 0.005 + 1.5e-6 * i      # 0.005 ... 0.020 kg/kg
        wet_bulb_temperature(T_db, omega)


def main() -> None:
    workload()  # warm up

    times = []
    for _ in range(N_RUNS):
        t0 = time.perf_counter()
        workload()
        times.append(time.perf_counter() - t0)

    mean = sum(times) / len(times)
    print(f"mean wall time over {N_RUNS} runs: {mean:.4f} s")
    print(f"individual runs: {['%.4f' % t for t in times]}")

    prof = cProfile.Profile()
    prof.enable()
    workload()
    prof.disable()
    s = StringIO()
    pstats.Stats(prof, stream=s).sort_stats("cumulative").print_stats(15)
    print("\ncProfile (top 15 by cumulative time):")
    print(s.getvalue())


if __name__ == "__main__":
    main()
