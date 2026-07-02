#include "observations/synth/obs_synth.hpp"

#include "stations/StationCatalog.hpp"
#include "perturbations/Shapiro.hpp"

#include <SpiceUsr.h>

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
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

[[nodiscard]] bool isPositiveFinite(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] std::string stationTargetName(const std::string& stationName) {
    try {
        return std::to_string(od::stationNaifIdFromName(stationName));
    } catch (const std::invalid_argument&) {
        return stationName;
    }
}

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

void validateGeometryConfig(const GeometryConfig& geometry) {
    if (geometry.target.empty()) {
        throw std::invalid_argument("Synthetic observation target cannot be empty.");
    }
    if (geometry.stationName.empty()) {
        throw std::invalid_argument("Synthetic observation station name cannot be empty.");
    }
    if (geometry.frame.empty()) {
        throw std::invalid_argument("Synthetic observation frame cannot be empty.");
    }
    if (geometry.aberrationCorrection.empty()) {
        throw std::invalid_argument("Synthetic observation aberration correction cannot be empty.");
    }
}

void validateNoiseConfig(const NoiseConfig& noise) {
    if (!std::isfinite(noise.sigma) || noise.sigma < 0.0) {
        throw std::invalid_argument("Synthetic observation noise sigma must be finite and non-negative.");
    }
}

} // namespace

SyntheticObservation::SyntheticObservation(GeometryConfig geometry, NoiseConfig noise)
    : geometry_(std::move(geometry)),
      noise_(noise),
      rng_(noise_.seed) {
    validateGeometryConfig(geometry_);
    validateNoiseConfig(noise_);
}

SyntheticObservation::~SyntheticObservation() = default;

void SyntheticObservation::setNoiseEnabled(bool enabled) noexcept {
    noise_.enabled = enabled;
}

void SyntheticObservation::setNoiseSigma(double sigma) {
    if (!std::isfinite(sigma) || sigma < 0.0) {
        throw std::invalid_argument("Synthetic observation noise sigma must be finite and non-negative.");
    }

    noise_.sigma = sigma;
}

void SyntheticObservation::reseed(std::uint32_t seed) {
    noise_.seed = seed;
    rng_.seed(seed);
}

const GeometryConfig& SyntheticObservation::geometryConfig() const noexcept {
    return geometry_;
}

const NoiseConfig& SyntheticObservation::noiseConfig() const noexcept {
    return noise_;
}

std::vector<double> SyntheticObservation::makeEpochGrid(double startTdb,
                                                        double endTdb,
                                                        double stepSeconds) {
    if (!std::isfinite(startTdb) || !std::isfinite(endTdb)) {
        throw std::invalid_argument("Synthetic observation time bounds must be finite.");
    }
    if (!isPositiveFinite(stepSeconds)) {
        throw std::invalid_argument("Synthetic observation step must be positive and finite.");
    }
    if (endTdb < startTdb) {
        throw std::invalid_argument("Synthetic observation end epoch cannot precede start epoch.");
    }

    std::vector<double> epochs;
    const double duration = endTdb - startTdb;
    const auto interval_count = static_cast<std::size_t>(std::floor(duration / stepSeconds));
    epochs.reserve(interval_count + 2);

    for (std::size_t index = 0; index <= interval_count; ++index) {
        epochs.push_back(startTdb + static_cast<double>(index) * stepSeconds);
    }

    if (epochs.empty() || std::abs(epochs.back() - endTdb) > 1.0e-9) {
        epochs.push_back(endTdb);
    }

    return epochs;
}

Eigen::Matrix<double, 6, 1> SyntheticObservation::relativeTargetState(double epochTdb) const {
    return relativeTargetGeometry(epochTdb).stationToTargetState;
}

RelativeGeometry SyntheticObservation::relativeTargetGeometry(double receiveEpochTdb) const {
    return relativeTargetGeometryFor(geometry_.target, geometry_.stationName, receiveEpochTdb);
}

RelativeGeometry SyntheticObservation::relativeTargetGeometryFor(const std::string& target,
                                                                 const std::string& stationName,
                                                                 double receiveEpochTdb) const {
    if (target.empty()) {
        throw std::invalid_argument("Synthetic observation target cannot be empty.");
    }
    if (stationName.empty()) {
        throw std::invalid_argument("Synthetic observation station name cannot be empty.");
    }
    if (!std::isfinite(receiveEpochTdb)) {
        throw std::invalid_argument("Synthetic observation epoch must be finite.");
    }

    SpiceErrorActionGuard action_guard;

    const std::string observer = stationTargetName(stationName);
    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;

    spkezr_c(target.c_str(),
             receiveEpochTdb,
             geometry_.frame.c_str(),
             geometry_.aberrationCorrection.c_str(),
             observer.c_str(),
             spice_state,
             &light_time);

    throwIfSpiceFailed("Failed to compute synthetic observation geometry for station " + stationName);

    RelativeGeometry geometry;
    geometry.lightTimeSec = light_time;
    geometry.receiveEpochTdb = receiveEpochTdb;
    geometry.emitEpochTdb = receiveEpochTdb - light_time;

    for (Eigen::Index i = 0; i < geometry.stationToTargetState.size(); ++i) {
        geometry.stationToTargetState[i] = spice_state[i];
    }

    return geometry;
}

Eigen::Matrix<double, 6, 1> SyntheticObservation::sunRelativeState(const std::string& body,
                                                                   double epochTdb) const {
    if (body.empty()) {
        throw std::invalid_argument("Synthetic observation Sun-relative body cannot be empty.");
    }
    if (!std::isfinite(epochTdb)) {
        throw std::invalid_argument("Synthetic observation Sun-relative epoch must be finite.");
    }

    SpiceErrorActionGuard action_guard;

    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;

    spkezr_c(body.c_str(),
             epochTdb,
             geometry_.frame.c_str(),
             "NONE",
             "SUN",
             spice_state,
             &light_time);

    throwIfSpiceFailed("Failed to compute Sun-relative synthetic observation state for " + body);

    Eigen::Matrix<double, 6, 1> state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spice_state[i];
    }

    return state;
}

double SyntheticObservation::shapiroRangeDelay(double receiveEpochTdb) const {
    return shapiroRangeDelay(relativeTargetGeometry(receiveEpochTdb));
}

double SyntheticObservation::shapiroRangeDelay(const RelativeGeometry& geometry) const {
    return shapiroRangeDelayFor(geometry_.target, geometry_.stationName, geometry);
}

double SyntheticObservation::shapiroRangeDelayFor(const std::string& target,
                                                  const std::string& stationName,
                                                  const RelativeGeometry& geometry) const {
    if (target.empty()) {
        throw std::invalid_argument("Synthetic observation Shapiro target cannot be empty.");
    }
    if (stationName.empty()) {
        throw std::invalid_argument("Synthetic observation Shapiro station cannot be empty.");
    }

    const std::string observer = stationTargetName(stationName);
    const Eigen::Matrix<double, 6, 1> observer_sun_state = sunRelativeState(observer,
                                                                            geometry.receiveEpochTdb);
    const Eigen::Matrix<double, 6, 1> target_sun_state = sunRelativeState(target,
                                                                          geometry.emitEpochTdb);

    const Eigen::Vector3d observer_from_sun = observer_sun_state.segment<3>(0);
    const Eigen::Vector3d target_from_sun = target_sun_state.segment<3>(0);
    const double r_obs_mag = observer_from_sun.norm();
    const double r_target_mag = target_from_sun.norm();
    const double rho_mag = (target_from_sun - observer_from_sun).norm();

    return fd::perturbations::computeShapiroRangeDelay(
        r_obs_mag,
        r_target_mag,
        rho_mag,
        fd::perturbations::kSunGravitationalParameterKm3PerSec2);
}

double SyntheticObservation::drawNoise() {
    if (!noise_.enabled || noise_.sigma == 0.0) {
        return 0.0;
    }

    std::normal_distribution<double> distribution(0.0, noise_.sigma);
    return distribution(rng_);
}

} // namespace fd::observations::synth
