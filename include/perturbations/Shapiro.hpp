#pragma once

#include <cmath>
#include <stdexcept>

namespace fd::perturbations {

inline constexpr double kSpeedOfLightKmPerSec = 299792.458;
inline constexpr double kSunGravitationalParameterKm3PerSec2 = 132712440041.9394;

inline double computeShapiroRangeDelay(double r_obs_mag,
                                       double r_target_mag,
                                       double rho_mag,
                                       double mu_central) {
    if (!std::isfinite(r_obs_mag) || !std::isfinite(r_target_mag)
        || !std::isfinite(rho_mag) || !std::isfinite(mu_central)) {
        throw std::runtime_error("Shapiro delay inputs must be finite.");
    }
    if (r_obs_mag <= 0.0 || r_target_mag <= 0.0 || rho_mag < 0.0 || mu_central <= 0.0) {
        throw std::runtime_error("Shapiro delay inputs must have physically valid magnitudes.");
    }

    const double numerator = r_obs_mag + r_target_mag + rho_mag;
    const double denominator = r_obs_mag + r_target_mag - rho_mag;

    if (denominator <= 0.0) {
        throw std::runtime_error("Shapiro delay denominator is non-positive; singular central-body geometry.");
    }
    if (numerator <= 0.0) {
        throw std::runtime_error("Shapiro delay numerator is non-positive.");
    }

    return (2.0 * mu_central / (kSpeedOfLightKmPerSec * kSpeedOfLightKmPerSec))
        * std::log(numerator / denominator);
}

inline double computeShapiroTimeDelay(double r_obs_mag,
                                      double r_target_mag,
                                      double rho_mag,
                                      double mu_central) {
    return computeShapiroRangeDelay(r_obs_mag, r_target_mag, rho_mag, mu_central)
        / kSpeedOfLightKmPerSec;
}

} // namespace fd::perturbations
