# Agent Handoff Notes

This file is for future coding agents working in this repository. It captures
the current project conventions, module map, and verification commands.

## Working Principles

- Read the current tree before editing. This checkout changes often.
- Do not revert user changes unless explicitly asked.
- Keep SPICE-heavy, loop-heavy, and non-template logic in `.cpp` files.
- Keep `Stations.hpp` as the stable, SPICE-agnostic station data model.
- Reusable library modules should not call `furnsh_c`; tests and demos load
  kernels.
- Prefer Eigen for vector/matrix math.
- Use km, s, km/s, and km/s^2 consistently.
- When adding SPICE calls, use `od::SpiceErrorModeGuard` and
  `od::throwIfSpiceFailed` from `include/utils/CSPICE/SpiceError.hpp`.

## Current Module Map

### Stations

- `include/stations/Stations.hpp`: base station model.
- `include/stations/StationCatalog.hpp`
- `include/stations/ElevationMask.hpp`
- `src/stations/StationCatalog.cpp`
- `src/stations/ElevationMask.cpp`

Station catalog helpers map DSN names such as `DSS-43` to NAIF station IDs and
can build station geometry from CSPICE kernels.
Elevation-mask helpers compute station-target elevation from CSPICE station SPK
and ITRF93 Earth-orientation geometry; this kernel set does not define
`DSS-43_TOPO`.

### Integrator

- `include/RKF45Integrator.hpp`
- `include/dynamics/EphemerisInterpolator.hpp`
- `src/RKF45Integrator.cpp`
- `src/dynamics/EphemerisInterpolator.cpp`
- `test_rkf45.cpp`

Adaptive RKF45 uses Eigen states and preallocated stage storage. The result
history can be loaded into the Hermite ephemeris interpolator for
propagate-once/interpolate-many observation evaluation.

### Filters

- `include/filters/filter.hpp`
- `include/filters/WLS.hpp`
- `src/filters/filter.cpp`
- `src/filters/WLS.cpp`

WLS implements a batch a-priori normal equation update. It uses solves rather
than explicitly forming inverse measurement covariance.

### Synthetic Observations

- `include/observations/synth/obs_synth.hpp`
- `include/observations/synth/RangeSynth.hpp`
- `include/observations/synth/RangeRateSynth.hpp`
- `src/observations/synth/obs_synth.cpp`
- `src/observations/synth/RangeSynth.cpp`
- `src/observations/synth/RangeRateSynth.cpp`
- `tests/test_synth_observations.cpp`
- `tests/test_voyager_position_od.cpp`

The generic `RangeSynth`/`RangeRateSynth` classes use CSPICE geometry, optional
Gaussian noise, and solar Shapiro delay. `RangeRateSynth` models Doppler as
centered differenced range over a configurable count time, defaulting to 60 s,
and scales its sigma by `1 / sqrt(count_time_s)`. The Voyager synthetic test
now seeds Voyager's initial state from CSPICE once, propagates it with RKF45 and
the OD dynamics, applies manual one-way light-time against that propagated
ephemeris, and keeps only samples above the station elevation mask. The Voyager
synthetic report is generated at:

```text
synthetic_observations_dss43_voyager1.txt
```

### Perturbations

- `include/perturbations/Shapiro.hpp`
- `include/perturbations/Gravitational.hpp`
- `include/perturbations/SRP.hpp`
- `src/perturbations/Gravitational.cpp`
- `src/perturbations/SRP.cpp`

Current perturbation classes:

- Solar Shapiro delay helpers.
- Third-body gravity with configurable massive bodies.
- Cannonball solar radiation pressure.

## Build and Test Commands

Targeted builds:

```sh
cmake --build build-clang --target test_rkf45 -j4
cmake --build build-clang --target test_synth_observations -j4
cmake --build build-clang --target test_voyager_position_od -j4
```

Targeted tests:

```sh
ctest --test-dir build-clang -R test_rkf45 --output-on-failure
ctest --test-dir build-clang -R test_synth_observations --output-on-failure
ctest --test-dir build-clang -R test_voyager_position_od --output-on-failure
```

SPICE-backed tests must run from `build-clang`, because `kernels.tm` contains:

```text
PATH_VALUES = ( '../Kernels' )
```

The canonical meta-kernel filename in the repository is lowercase `kernels.tm`.
Older files may still refer to `../Kernels.tm`; check and fix the literal before
trusting an all-test run on a case-sensitive filesystem.

## Known Runtime Setup

Environment:

```sh
export CSPICE_HOME=/path/to/cspice
```

Required kernel inventory is referenced by `kernels.tm`, including:

- `naif0012.tls`
- `pck00011.tpc`
- `de442.bsp`
- `vgr1_jup230.bsp`
- `jup310.bsp`
- `Voyager_1.a54206u_V0.2_merged.bsp`
- `gm_de440.tpc`
- DSN station and Earth orientation kernels

`jup310.bsp` is required by `test_voyager_position_od`; it provides the 1979
coverage needed for Jupiter body `599` and the Galilean moons `501-504`.
Earlier attempts with `jup230`/`jup363` did not provide usable 1979 moon states
for this test.

## Current Verified Values

`test_synth_observations`:

```text
Stations: DSS-43, DSS-63
Target: Voyager 1 (-31)
Requested arc: 1979-01-10T00:00:00 to 1979-01-14T00:00:00
Elevation mask: 10 deg
Cadence: 180 s
VLBI cadence: 20 s
Samples after mask: 1655 total; DSS-43 659, DSS-63 996
VLBI samples after mask: 2902 total; DSS-43/DSS-63 82, DSS-14/DSS-43 2820
Range sigma: 0.010 km
Range-rate base sigma: 1.0e-6 km/s at 1 s
Range-rate count time: 60 s
Range-rate row sigma: 1.290994e-7 km/s
VLBI sigma: 0.001 km
```

`test_voyager_position_od` reads that report and estimates the initial Voyager
state with DP853 propagation and cached SPICE body/station states:
Sun gravity, Jupiter body `599`, Galilean moons `501-504`, Saturn barycenter
`6`, SRP, light-time, Shapiro delay, centered count-time range-rate, and VLBI
differential range. Its debug report is generated at
`tests/voyager_position_estimation_report.txt`. It also exports plotting CSVs:
`tests/voyager_od_postfit_diagnostics_VLBI.csv`,
`tests/voyager_od_trajectory_error_VLBI.csv`, and
`tests/voyager_station_observability_windows_VLBI.csv`.

Latest verified `jup310.bsp` run:

```text
Active third bodies: 599 501 502 503 504 6
Prior position error:      70.710678 km
Posterior position error:  5.766723 km
Solver status:             CONVERGED in 5 iterations
Solver weighted RMS:       0.694371
Prefit range RMS:          150770.580 m
Postfit range RMS:         9.671 m
Prefit range-rate RMS:     528.089 mm/s
Postfit range-rate RMS:    0.179 mm/s
Prefit VLBI RMS:           1.598 m
Postfit VLBI RMS:          1.005 m
OD test runtime:           about 4 s
```

## Worktree Caution

At the time this handoff was written, the worktree had unrelated local changes.
In particular, `include/linearization.hpp` appeared deleted. Do not restore or
remove unrelated files unless the user explicitly asks.

## Next Logical Steps

- Promote the `test_voyager_position_od` computed-observation and dynamics code
  into reusable modules once the interfaces settle.
- Add SVD, covariance, and conditioning diagnostics around the finite-difference
  design matrix.
- Continue checking final-state accuracy against SPICE after post-fit residuals
  improve; a short two-station arc can still leave weakly observed directions.
