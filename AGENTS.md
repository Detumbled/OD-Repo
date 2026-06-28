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

Synthetic range/range-rate uses CSPICE geometry, optional Gaussian noise, and
solar Shapiro delay. The DSS-43/Voyager report is generated at:

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
```

Targeted tests:

```sh
ctest --test-dir build-clang -R test_rkf45 --output-on-failure
ctest --test-dir build-clang -R test_wls --output-on-failure
ctest --test-dir build-clang -R test_synth_observations --output-on-failure
ctest --test-dir build-clang -R test_perturbations --output-on-failure
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
- `Voyager_1.a54206u_V0.2_merged.bsp`
- `gm_de440.tpc`
- DSN station and Earth orientation kernels

## Current Verified Values

`test_perturbations` at `1979-03-05T00:00:00`:

```text
Jupiter third-body acceleration norm = 1.463970346623e-04 km/s^2
SRP acceleration norm                = 2.941469913405e-12 km/s^2
```

`test_synth_observations`:

```text
Station: DSS-43
Target: Voyager 1 (-31)
Arc: 1979-03-05T00:00:00 to 1979-03-05T01:00:00
Cadence: 300 s
Range sigma: 0.010 km
Range-rate sigma: 1.0e-6 km/s
```

## Worktree Caution

At the time this handoff was written, the worktree had unrelated local changes.
In particular, `include/linearization.hpp` appeared deleted. Do not restore or
remove unrelated files unless the user explicitly asks.

## Next Logical Steps

- Connect `ThirdBodyGravity` and `SolarRadiationPressure` to a dynamics builder
  used by `RKF45Integrator`.
- Add a computed-observation module with iterative light-time and Shapiro delay
  using the same `perturbations/Shapiro.hpp` functions.
- Extend WLS tests from synthetic linearized range to full computed observation
  residuals.
- Add covariance and conditioning diagnostics when the estimator is connected to
  real observation partials.
