#include "perturbations/SRP.hpp"

#include <SpiceUsr.h>

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace fd::perturbations {
namespace {

using SpiceErrorActionBuffer = std::array<SpiceChar, 32>;

class SpiceErrorActionGuard {
public:
    SpiceErrorActionGuard() {
        erract_c("GET", static_cast<SpiceInt>(previous_action_.size()), previous_action_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
    }

    SpiceErrorActionGuard(const SpiceErrorActionGuard&) = delete;
    SpiceErrorActionGuard& operator=(const SpiceErrorActionGuard&) = delete;

    ~SpiceErrorActionGuard() {
        erract_c("SET", 0, previous_action_.data());
    }

private:
    SpiceErrorActionBuffer previous_action_ {};
};

void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar short_message[1841] = {0};
    SpiceChar long_message[1841] = {0};
    getmsg_c("SHORT", sizeof(short_message), short_message);
    getmsg_c("LONG", sizeof(long_message), long_message);
    reset_c();

    std::ostringstream message;
    message << context << ": " << short_message << " | " << long_message;
    throw std::runtime_error(message.str());
}

void validateReflectivity(double cr) {
    if (!std::isfinite(cr) || cr < 0.0) {
        throw std::invalid_argument("SRP reflectivity coefficient must be finite and non-negative.");
    }
}

void validatePositiveFinite(double value, const char* label) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(std::string(label) + " must be positive and finite.");
    }
}

[[nodiscard]] Eigen::Vector3d sunPositionRelativeToCentral(double tdb,
                                                           const std::string& frame,
                                                           const std::string& centralBody) {
    if (centralBody == "SUN" || centralBody == "10") {
        return Eigen::Vector3d::Zero();
    }

    SpiceErrorActionGuard action_guard;

    SpiceDouble position[3] = {0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkpos_c("SUN",
             tdb,
             frame.c_str(),
             "NONE",
             centralBody.c_str(),
             position,
             &light_time);

    throwIfSpiceFailed("Failed to compute Sun position for SRP.");

    return {position[0], position[1], position[2]};
}

} // namespace

SolarRadiationPressure::SolarRadiationPressure(std::string centralBody,
                                               std::string frame,
                                               double cr,
                                               double area_m2,
                                               double mass_kg)
    : centralBody_(std::move(centralBody)),
      frame_(std::move(frame)) {
    if (centralBody_.empty()) {
        throw std::invalid_argument("SRP central body cannot be empty.");
    }
    if (frame_.empty()) {
        throw std::invalid_argument("SRP frame cannot be empty.");
    }

    setReflectivity(cr);
    setArea(area_m2);
    setMass(mass_kg);
}

void SolarRadiationPressure::setReflectivity(double cr) {
    validateReflectivity(cr);
    cr_ = cr;
}

void SolarRadiationPressure::setArea(double area_m2) {
    validatePositiveFinite(area_m2, "SRP area");
    area_m2_ = area_m2;
}

void SolarRadiationPressure::setMass(double mass_kg) {
    validatePositiveFinite(mass_kg, "SRP mass");
    mass_kg_ = mass_kg;
}

Eigen::Vector3d SolarRadiationPressure::computeAcceleration(double tdb,
                                                            const Eigen::Vector3d& scPosition) const {
    if (!std::isfinite(tdb) || !scPosition.allFinite()) {
        throw std::invalid_argument("SRP acceleration inputs must be finite.");
    }

    const Eigen::Vector3d central_to_sun = sunPositionRelativeToCentral(tdb, frame_, centralBody_);
    const Eigen::Vector3d sun_to_spacecraft = scPosition - central_to_sun;
    const double sun_spacecraft_distance_km = sun_to_spacecraft.norm();
    if (sun_spacecraft_distance_km <= 0.0) {
        throw std::runtime_error("SRP geometry has zero Sun-spacecraft separation.");
    }

    const Eigen::Vector3d anti_sunward_unit = sun_to_spacecraft / sun_spacecraft_distance_km;
    const double acceleration_m_per_s2 = kSolarPressure1AU_N_m2
        * cr_
        * (area_m2_ / mass_kg_)
        * std::pow(kAu_km / sun_spacecraft_distance_km, 2.0);

    return (acceleration_m_per_s2 / 1000.0) * anti_sunward_unit;
}

} // namespace fd::perturbations
