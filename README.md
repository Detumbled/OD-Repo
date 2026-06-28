# OD-Repo

Orbit-determination and flight-dynamics sandbox for DSN/Voyager tracking work.
The current implementation is centered on C++/Eigen numerical utilities,
CSPICE-backed station geometry, synthetic observations, perturbation models, and
batch weighted least-squares filtering.

## Current Scope

Implemented modules:

- CSPICE-backed DSN station catalog.
- Adaptive RKF45 integrator using Eigen vectors and preallocated stage storage.
- Batch Weighted Least Squares filter with a priori information.
- Synthetic one-way range and range-rate observations.
- Solar Shapiro delay applied to synthetic radiometric observables.
- Modular perturbation models for third-body gravity and cannonball SRP.
- Focused tests for RKF45, WLS, synthetic observations, perturbations,
  station-kernel sanity checks, and a Voyager 1 position OD smoke case.

Planned or partial areas:

- Computed-observation light-time iteration integrated with the estimator.
- Full dynamics builder that sums central gravity, third bodies, SRP, and future force models.
- Sequential filters and richer OD diagnostics.
- Attitude/AOCS and visualization work.

## Repository Layout

```text
include/
  RKF45Integrator.hpp
  StationCatalog.hpp
  Stations.hpp
  filters/
    filter.hpp
    WLS.hpp
  observations/synth/
    obs_synth.hpp
    RangeSynth.hpp
    RangeRateSynth.hpp
  perturbations/
    Gravitational.hpp
    SRP.hpp
    Shapiro.hpp

src/
  RKF45Integrator.cpp
  StationCatalog.cpp
  filters/
    filter.cpp
    WLS.cpp
  observations/synth/
    obs_synth.cpp
    RangeSynth.cpp
    RangeRateSynth.cpp
  perturbations/
    Gravitational.cpp
    SRP.cpp

tests/
  test_wls.cpp
  test_synth_observations.cpp
  test_perturbations.cpp
  test_voyager_position_od.cpp

test_rkf45.cpp
test_stations.cpp
station_catalog_demo.cpp
kernels.tm
```

## Dependencies

The CMake project currently uses C++23 and links these dependencies:

- Eigen3
- CSPICE
- Ceres
- GLFW, GLEW, OpenGL
- local `third_party/imgui`, `third_party/implot`, and `third_party/glm`

`CSPICE_HOME` must point to the CSPICE installation:

```sh
export CSPICE_HOME=/path/to/cspice
```

The project expects the Voyager/DSN kernels under `Kernels/` and the meta-kernel
at `kernels.tm`. The Voyager Jupiter-encounter OD test requires `jup310.bsp`
for 1979 coverage of Jupiter body `599` and Galilean moons `501-504`.

Important path convention: `kernels.tm` uses `PATH_VALUES = ( '../Kernels' )`.
For tests or demos that load `../kernels.tm`, run from `build-clang`.

## Build

```sh
cmake --build build-clang --target test_rkf45 -j4
cmake --build build-clang --target test_wls -j4
cmake --build build-clang --target test_synth_observations -j4
cmake --build build-clang --target test_perturbations -j4
cmake --build build-clang --target test_voyager_position_od -j4
cmake --build build-clang --target station_catalog_demo -j4
```

The current `CMakeLists.txt` filters source candidates with `EXISTS` so missing
legacy examples do not break generation.

## Tests

Run focused tests:

```sh
ctest --test-dir build-clang -R test_rkf45 --output-on-failure
ctest --test-dir build-clang -R test_wls --output-on-failure
ctest --test-dir build-clang -R test_synth_observations --output-on-failure
ctest --test-dir build-clang -R test_perturbations --output-on-failure
ctest --test-dir build-clang -R test_voyager_position_od --output-on-failure
```

SPICE-backed tests are registered with `build-clang` as their working directory
so `../kernels.tm` resolves correctly.

`test_stations.cpp` is an older kernel smoke test and currently hard-codes
`../Kernels.tm`. On case-sensitive filesystems, either update that literal to
`../kernels.tm` or provide a matching compatibility file before relying on it in
an all-test run.

## Implemented Components

### Station Catalog

`StationCatalog.hpp/.cpp` provides:

- `defaultDsnStationCatalog()`
- `stationNaifIdFromName(...)`
- `buildStationFromKernel(...)`
- `buildDefaultDsnCatalogFromKernel(...)`

Kernel loading remains the caller's responsibility. The catalog code queries
CSPICE and returns `od::Station` objects without mutating the base station model.

### RKF45 Integrator

`RKF45Integrator` uses:

- `Eigen::VectorXd` state representation.
- `Eigen::Ref` derivative callbacks.
- Fehlberg 4(5) embedded error estimate.
- Adaptive step acceptance/rejection.
- Preallocated trial, fourth-order, fifth-order, and stage matrices.

`test_rkf45` propagates a simple circular Kepler orbit for one period and checks
position closure, velocity closure, energy, angular momentum, and runtime.

### WLS Filter

`fd::filters::WLS` implements a batch Weighted Least Squares update with a priori
information:

```text
Lambda = H^T R^-1 H + P0^-1
N      = H^T R^-1 r + P0^-1 (x_prior - x_nominal)
dx     = Lambda^-1 N
```

The implementation uses Eigen solves rather than explicitly forming `R^-1`.
LDLT is attempted first, with QR fallback for robustness.

`test_wls` creates an 8-hour synthetic range arc with 10 m noise and verifies
that the filter reduces both linearized and nonlinear range residual RMS.

### Synthetic Observations

Synthetic observation classes live under `fd::observations::synth`.

Shared features:

- `GeometryConfig`: target, station, frame, aberration correction.
- `NoiseConfig`: enabled flag, sigma, random seed.
- Configurable time grid via `makeEpochGrid(start, end, step)`.
- CSPICE station names resolved through the DSN station catalog.

`RangeSynth`:

- Gets light-time corrected station-to-target geometry through CSPICE.
- Adds solar Shapiro range delay.
- Adds optional Gaussian range noise.

`RangeRateSynth`:

- Computes geometric one-way range-rate from relative position/velocity.
- Adds a finite-difference derivative of solar Shapiro delay using `dt = 1 s`.
- Adds optional Gaussian range-rate noise.

`test_synth_observations` generates interleaved station observations:

- Stations: DSS-43, DSS-63
- Target: Voyager 1 (`-31`)
- Start: `1979-03-05T00:00:00`
- Duration: 3 hours
- Cadence: 3 minutes
- Samples: 122 total, 61 per station
- Range sigma: 0.010 km
- Range-rate sigma: `1.0e-6 km/s`

The output report keeps the historical filename but now includes a `station`
column:

```text
synthetic_observations_dss43_voyager1.txt
```

### Perturbations

`fd::perturbations::Shapiro` provides inline reusable functions:

- `computeShapiroRangeDelay(...)`
- `computeShapiroTimeDelay(...)`

`fd::perturbations::ThirdBodyGravity` computes:

```text
a_3rd = mu_i * (r_sc_to_i / |r_sc_to_i|^3 - r_central_to_i / |r_central_to_i|^3)
```

`fd::perturbations::SolarRadiationPressure` computes cannonball SRP:

- anti-sunward direction
- configurable `C_R`, area, and mass
- returns `km/s^2`

`test_perturbations` evaluates both at `1979-03-05T00:00:00` for Voyager 1:

```text
Jupiter third-body acceleration norm = 1.463970346623e-04 km/s^2
SRP acceleration norm                = 2.941469913405e-12 km/s^2
```

## Units

Unless otherwise noted:

- distance: km
- time: s
- velocity: km/s
- acceleration: km/s^2
- SRP pressure: N/m^2 internally, converted to km/s^2 at output

## SPICE Notes

- Reusable modules do not call `furnsh_c`.
- Tests and demos load meta-kernels explicitly.
- The canonical repository meta-kernel is lowercase `kernels.tm`; some legacy
  files may still mention `Kernels.tm`.
- SPICE error mode is temporarily set to `RETURN` inside reusable modules that
  query CSPICE, then restored with RAII guards.
- Pure dynamics perturbation geometry uses no aberration correction (`"NONE"`).
- Synthetic observations currently default to light-time correction (`"LT"`) and
  then add Shapiro delay explicitly.

## Generated Files

`test_synth_observations` writes:

```text
synthetic_observations_dss43_voyager1.txt
```

This file is useful for inspecting synthetic observation values and for checking
whether Shapiro/noise settings changed output.

`test_voyager_position_od` consumes that report and estimates the Voyager 1
initial state using RKF45 with Sun gravity, Jupiter body `599`, Galilean moons
`501-504`, Saturn barycenter `6`, SRP, light-time, and Shapiro delay. It writes:

```text
tests/voyager_position_estimation_report.txt
```

Latest verified `jup310.bsp` run:

```text
Active third bodies: 599 501 502 503 504 6
Prior position error:      70.710678 km
Posterior position error:  51.587303 km
Range RMS:                 37308.328 m -> 161.648 m
Range-rate RMS:            457.768 mm/s -> 13.464 mm/s
```

This is a short-arc diagnostic test, not a final high-precision OD solution.
Use the report residuals and future conditioning diagnostics before interpreting
state accuracy from residual reduction alone.
