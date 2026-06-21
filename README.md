# Deep Space Navigation and Attitude Simulation Framework

## Overview

This project aims to develop a high-fidelity simulation framework that combines:

* Orbit propagation
* Orbit determination
* Attitude dynamics
* Attitude determination and control (AOCS)
* Deep-space antenna pointing
* Synthetic measurement generation
* Navigation performance analysis

The primary mission scenario consists of a spacecraft equipped with a steerable high-gain antenna tracking Voyager 1 while simultaneously transmitting navigation information to Earth.

The simulation environment is built around NASA SPICE kernels and provides a complete chain from truth-state generation to navigation-state estimation.

---

# Mission Scenario

## Concept

A deep-space spacecraft acts as a relay and observation platform.

The spacecraft:

1. Tracks Voyager 1.
2. Maintains antenna lock on Voyager.
3. Generates synthetic observations.
4. Downlinks measurements to Earth.
5. Performs onboard orbit determination.

The framework compares:

* Spaceborne measurements
* Ground-based measurements
* Estimated states
* Truth states from SPICE

---

# Main Objectives

## Orbit Determination

Implement and compare:

### Initial Orbit Determination

* Gauss Method
* Gibbs Method
* Herrick-Gibbs Method
* Lambert-based estimation

### Statistical Orbit Determination

* Weighted Least Squares
* Batch Least Squares
* Sequential Least Squares
* Extended Kalman Filter
* Unscented Kalman Filter

Based primarily on:

* Tapley, Statistical Orbit Determination

---

## Attitude Determination and Control

Implement:

### Reference Frames

* ECI
* ECEF
* LVLH
* Hill Frame
* Spacecraft Body Frame

### Attitude Representations

* Direction Cosine Matrices
* Euler Angles
* Modified Rodrigues Parameters
* Quaternions

### Attitude Dynamics

Rigid-body equations:

Iω̇ + ω × (Iω) = τ

### Attitude Controllers

* PD Controller
* Quaternion Feedback Control
* Reaction Wheel Control

Based primarily on:

* Schaub & Junkins

---

# Software Architecture

## Module Structure

```text
include/

    Orbit/
        Propagator.hpp
        ForceModels.hpp
        Lambert.hpp
        Gibbs.hpp
        BatchLS.hpp
        EKF.hpp

    Attitude/
        Quaternion.hpp
        DCM.hpp
        AttitudeDynamics.hpp
        Controllers.hpp

    Sensors/
        StarTracker.hpp
        Gyroscope.hpp
        Antenna.hpp

    Measurements/
        Range.hpp
        Doppler.hpp
        AngleMeasurement.hpp

    Spice/
        SpiceManager.hpp
        FrameTransform.hpp

    Visualization/
        Camera.hpp
        SpacecraftRenderer.hpp

src/

    orbit/
    attitude/
    sensors/
    measurements/
    spice/
    visualization/
```

---

# Phase 1 — Truth Model

## SPICE Integration

Load:

* Planetary ephemerides
* Voyager kernels
* Earth orientation kernels
* Leap second kernels

Primary CSPICE functions:

* furnsh_c()
* spkezr_c()
* sxform_c()
* pxform_c()
* str2et_c()

Outputs:

* Spacecraft state
* Voyager state
* Earth state
* Frame transformations

---

# Phase 2 — Orbit Propagation

## Two-Body Dynamics

State vector:

x = [r,v]

Equation:

r¨ = -μr/r³

---

## Perturbation Models

### Earth Missions

* J2
* J3
* Atmospheric Drag
* Solar Radiation Pressure

### Deep-Space Missions

* Third-body perturbations
* Solar Radiation Pressure
* Relativistic corrections (optional)

---

## Numerical Integrators

Implement:

### Fixed Step

* RK4

### Adaptive

* RKF45
* Dormand-Prince 8(7)

---

# Phase 3 — Attitude Simulation

## Spacecraft Model

State:

x_att = [q,ω]

where:

q = quaternion

ω = body angular velocity

---

## Actuators

### Reaction Wheels

Model:

τ_rw

### Thrusters

Optional future extension.

---

## Sensors

### Star Tracker

Output:

q_measured

Noise:

Gaussian attitude error

---

### Gyroscope

Output:

ω_measured

Include:

* Bias
* Scale factor
* White noise

---

# Phase 4 — Antenna Pointing

## High Gain Antenna

Define:

* Beam width
* Maximum slew rate
* Boresight vector

---

## Pointing Error

Compute:

θ = acos(b · LOS)

where:

LOS is the line-of-sight vector to Voyager.

---

## Performance Metrics

Track:

* Pointing accuracy
* Time inside beam
* Maximum pointing error
* RMS pointing error

---

# Phase 5 — Measurement Generation

## Angular Measurements

Generate:

* Right Ascension
* Declination

Measurement model:

z = h(x) + v

---

## Range Measurements

ρ = ||r_target-r_observer||

Noise model:

σ_range

---

## Doppler Measurements

fD = -(vLOS/c)f0

Noise model:

σ_doppler

---

## Antenna-Based Measurements

Estimate:

* Signal strength
* Gain loss
* Pointing offset

---

# Phase 6 — Orbit Determination

## Batch Least Squares

Implement:

minimize:

J = rᵀWr

where:

r = residual vector

W = weighting matrix

---

## Covariance Analysis

Compute:

P = (HᵀWH)^(-1)

Analyze:

* Position uncertainty
* Velocity uncertainty
* Correlation structure

---

## Extended Kalman Filter

Propagation:

x̂(k+1)

Covariance:

P(k+1)

Measurement update:

K = PHᵀ(HPHᵀ+R)^(-1)

---

# Phase 7 — Earth Comparison

Simulate observations from:

* Goldstone
* Madrid
* Canberra

Compare:

* Ground-only OD
* Spacecraft-only OD
* Combined OD

Metrics:

* RMS position error
* RMS velocity error
* Covariance volume
* Convergence time

---

# Visualization

## 3D Scene

Display:

* Sun
* Earth
* Voyager
* Observer spacecraft

Render:

* Trajectories
* Attitude axes
* Antenna beam cone
* Measurement rays

---

## Analysis Plots

Generate:

* Position errors
* Velocity errors
* Attitude errors
* Covariance evolution
* Residual histories

---

# Long-Term Extensions

## Autonomous Navigation

* Optical navigation
* Beacon tracking
* Relative navigation

## Multi-Spacecraft Scenarios

* Formation flying
* Relay constellations
* Interplanetary navigation networks

## Advanced Filters

* Unscented Kalman Filter
* Square Root Information Filter
* Particle Filter

---

# Expected Outcomes

The project should provide:

1. A reusable astrodynamics library.
2. A reusable AOCS library.
3. A complete orbit determination framework.
4. A deep-space navigation simulator.
5. A SPICE-based mission analysis environment.
6. Validation against real ephemerides and tracking geometries.
7. A portfolio-grade aerospace software project approaching professional mission-analysis tools.
   """
