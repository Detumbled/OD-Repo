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
- When adding SPICE calls, use `erract_c("RETURN")` with an RAII-style guard and
  restore the previous state.

## Current Module Map

### Stations

- `include/Stations.hpp`: base station model.
- `include/StationCatalog.hpp`
- `src/StationCatalog.cpp`

Station catalog helpers map DSN names such as `DSS-43` to NAIF station IDs and
can build station geometry from CSPICE kernels.

### Integrator

- `include/RKF45Integrator.hpp`
- `src/RKF45Integrator.cpp`
- `test_rkf45.cpp`

Adaptive RKF45 uses Eigen states and preallocated stage storage.

### Filters

- `include/filters/filter.hpp`
- `include/filters/WLS.hpp`
- `src/filters/filter.cpp`
- `src/filters/WLS.cpp`
- `tests/test_wls.cpp`

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

Synthetic range/range-rate uses CSPICE geometry, optional Gaussian noise, and
solar Shapiro delay. The DSS-43/DSS-63 Voyager report is generated at:

```text
synthetic_observations_dss43_voyager1.txt
```

### Perturbations

- `include/perturbations/Shapiro.hpp`
- `include/perturbations/Gravitational.hpp`
- `include/perturbations/SRP.hpp`
- `src/perturbations/Gravitational.cpp`
- `src/perturbations/SRP.cpp`
- `tests/test_perturbations.cpp`

Current perturbation classes:

- Solar Shapiro delay helpers.
- Third-body gravity with configurable massive bodies.
- Cannonball solar radiation pressure.

## Build and Test Commands

Targeted builds:

```sh
cmake --build build-clang --target test_rkf45 -j4
cmake --build build-clang --target test_wls -j4
cmake --build build-clang --target test_synth_observations -j4
cmake --build build-clang --target test_perturbations -j4
cmake --build build-clang --target test_voyager_position_od -j4
```

Targeted tests:

```sh
ctest --test-dir build-clang -R test_rkf45 --output-on-failure
ctest --test-dir build-clang -R test_wls --output-on-failure
ctest --test-dir build-clang -R test_synth_observations --output-on-failure
ctest --test-dir build-clang -R test_perturbations --output-on-failure
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

`test_perturbations` at `1979-03-05T00:00:00`:

```text
Jupiter third-body acceleration norm = 1.463970346623e-04 km/s^2
SRP acceleration norm                = 2.941469913405e-12 km/s^2
```

`test_synth_observations`:

```text
Stations: DSS-43, DSS-63
Target: Voyager 1 (-31)
Arc: 1979-03-05T00:00:00 to 1979-03-05T03:00:00
Cadence: 180 s
Samples: 122 total, 61 per station
Range sigma: 0.010 km
Range-rate sigma: 1.0e-6 km/s
```

`test_voyager_position_od` reads that report and estimates the initial Voyager
state with RKF45, Sun gravity, Jupiter body `599`, Galilean moons `501-504`,
Saturn barycenter `6`, SRP, light-time, and Shapiro delay. Its debug report is
generated at `tests/voyager_position_estimation_report.txt`.

Latest verified `jup310.bsp` run:

```text
Active third bodies: 599 501 502 503 504 6
Prior position error:      70.710678 km
Posterior position error:  51.587303 km
Prefit range RMS:          37308.328 m
Postfit range RMS:         161.648 m
Prefit range-rate RMS:     457.768 mm/s
Postfit range-rate RMS:    13.464 mm/s
OD test runtime:           about 129 s
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
