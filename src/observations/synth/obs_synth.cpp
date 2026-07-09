#include "observations/synth/obs_synth.hpp"

#include "stations/StationCatalog.hpp"
#include "perturbations/Shapiro.hpp"
#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

constexpr double kLightTimeToleranceSec = 1.0e-9;
constexpr int kMaxLightTimeIterations = 12;

using od::throwIfSpiceFailed;

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

[[nodiscard]] Eigen::Vector3d rotateToItrf93(const Eigen::Vector3d& vector,
                                             double epochTdb,
                                             const std::string& frame) {
    SpiceDouble rotation[3][3] = {};
    pxform_c(frame.c_str(), "ITRF93", epochTdb, rotation);

    return Eigen::Vector3d{
        rotation[0][0] * vector.x() + rotation[0][1] * vector.y() + rotation[0][2] * vector.z(),
        rotation[1][0] * vector.x() + rotation[1][1] * vector.y() + rotation[1][2] * vector.z(),
        rotation[2][0] * vector.x() + rotation[2][1] * vector.y() + rotation[2][2] * vector.z()
    };
}

[[nodiscard]] Eigen::Vector3d stationUpVectorItrf(const od::Station& station) {
    const double cos_lat = std::cos(station.latitudeRad());
    return Eigen::Vector3d{
        cos_lat * std::cos(station.longitudeRad()),
        cos_lat * std::sin(station.longitudeRad()),
        std::sin(station.latitudeRad())
    };
}

[[nodiscard]] double shapiroRangeDelayForPositions(const Eigen::Vector3d& observerFromSun,
                                                   const Eigen::Vector3d& targetFromSun) {
    const double rho_mag = (targetFromSun - observerFromSun).norm();
    return fd::perturbations::computeShapiroRangeDelay(
        observerFromSun.norm(),
        targetFromSun.norm(),
        rho_mag,
        fd::perturbations::kSunGravitationalParameterKm3PerSec2);
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

double SyntheticObservation::targetElevationRad(
    double receiveEpochTdb,
    const od::Station& station,
    const TargetStateProvider& targetProvider) const {
    return targetElevationRadFor(geometry_.stationName, station, receiveEpochTdb, targetProvider);
}

bool SyntheticObservation::isVisibleAboveElevation(
    double receiveEpochTdb,
    const od::Station& station,
    double minimumElevationRad,
    const TargetStateProvider& targetProvider) const {
    if (!std::isfinite(minimumElevationRad)) {
        throw std::invalid_argument("Synthetic observation elevation mask must be finite.");
    }

    return targetElevationRad(receiveEpochTdb, station, targetProvider) >= minimumElevationRad;
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

RelativeGeometry SyntheticObservation::relativeTargetGeometry(
    double receiveEpochTdb,
    const TargetStateProvider& targetProvider) const {
    return relativeTargetGeometryFor(geometry_.stationName, receiveEpochTdb, targetProvider);
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

    const SpiceTargetStateProvider target_provider(target, geometry_.frame, "SUN");
    return relativeTargetGeometryFor(stationName, receiveEpochTdb, target_provider);
}

RelativeGeometry SyntheticObservation::relativeTargetGeometryFor(
    const std::string& stationName,
    double receiveEpochTdb,
    const TargetStateProvider& targetProvider) const {
    if (stationName.empty()) {
        throw std::invalid_argument("Synthetic observation station name cannot be empty.");
    }
    if (!std::isfinite(receiveEpochTdb)) {
        throw std::invalid_argument("Synthetic observation epoch must be finite.");
    }

    const std::string observer = stationTargetName(stationName);
    const Eigen::Matrix<double, 6, 1> station_state = sunRelativeState(observer, receiveEpochTdb);
    const Eigen::Vector3d station_position = station_state.segment<3>(0);

    double light_time_sec = 0.0;
    Eigen::Matrix<double, 6, 1> target_state = targetProvider.stateAt(receiveEpochTdb);

    for (int iteration = 0; iteration < kMaxLightTimeIterations; ++iteration) {
        const Eigen::Vector3d target_position = target_state.segment<3>(0);
        const double geometric_range_km = (target_position - station_position).norm();
        const double shapiro_delay_km =
            shapiroRangeDelayForPositions(station_position, target_position);
        const double updated_light_time_sec =
            (geometric_range_km + shapiro_delay_km)
            / fd::perturbations::kSpeedOfLightKmPerSec;

        if (std::abs(updated_light_time_sec - light_time_sec) < kLightTimeToleranceSec) {
            light_time_sec = updated_light_time_sec;
            break;
        }

        light_time_sec = updated_light_time_sec;
        target_state = targetProvider.stateAt(receiveEpochTdb - light_time_sec);
    }

    target_state = targetProvider.stateAt(receiveEpochTdb - light_time_sec);

    RelativeGeometry geometry;
    geometry.stationToTargetState = target_state - station_state;
    geometry.lightTimeSec = light_time_sec;
    geometry.receiveEpochTdb = receiveEpochTdb;
    geometry.emitEpochTdb = receiveEpochTdb - light_time_sec;

    return geometry;
}

double SyntheticObservation::targetElevationRadFor(
    const std::string& stationName,
    const od::Station& station,
    double receiveEpochTdb,
    const TargetStateProvider& targetProvider) const {
    if (stationName.empty()) {
        throw std::invalid_argument("Synthetic observation station name cannot be empty.");
    }
    if (!std::isfinite(receiveEpochTdb)) {
        throw std::invalid_argument("Synthetic observation epoch must be finite.");
    }

    const std::string observer = stationTargetName(stationName);
    const Eigen::Matrix<double, 6, 1> station_state = sunRelativeState(observer, receiveEpochTdb);
    const Eigen::Matrix<double, 6, 1> target_state = targetProvider.stateAt(receiveEpochTdb);
    const Eigen::Vector3d line_of_sight =
        target_state.segment<3>(0) - station_state.segment<3>(0);
    const Eigen::Vector3d line_of_sight_itrf =
        rotateToItrf93(line_of_sight, receiveEpochTdb, geometry_.frame);
    throwIfSpiceFailed("Failed to transform synthetic observation line of sight into ITRF93");

    const double range_km = line_of_sight_itrf.norm();
    if (!(range_km > 0.0) || !std::isfinite(range_km)) {
        throw std::runtime_error("Cannot compute elevation for zero station-target range.");
    }

    const double sine_elevation =
        std::clamp(line_of_sight_itrf.normalized().dot(stationUpVectorItrf(station)), -1.0, 1.0);
    return std::asin(sine_elevation);
}

Eigen::Matrix<double, 6, 1> SyntheticObservation::sunRelativeState(const std::string& body,
                                                                   double epochTdb) const {
    if (body.empty()) {
        throw std::invalid_argument("Synthetic observation Sun-relative body cannot be empty.");
    }
    if (!std::isfinite(epochTdb)) {
        throw std::invalid_argument("Synthetic observation Sun-relative epoch must be finite.");
    }

    od::SpiceErrorModeGuard action_guard;

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
    return shapiroRangeDelayForPositions(observer_from_sun, target_from_sun);
}

double SyntheticObservation::shapiroRangeDelayFor(
    const std::string& stationName,
    const RelativeGeometry& geometry,
    const TargetStateProvider& targetProvider) const {
    if (stationName.empty()) {
        throw std::invalid_argument("Synthetic observation Shapiro station cannot be empty.");
    }

    const std::string observer = stationTargetName(stationName);
    const Eigen::Matrix<double, 6, 1> observer_sun_state = sunRelativeState(observer,
                                                                            geometry.receiveEpochTdb);
    const Eigen::Matrix<double, 6, 1> target_sun_state = targetProvider.stateAt(geometry.emitEpochTdb);

    const Eigen::Vector3d observer_from_sun = observer_sun_state.segment<3>(0);
    const Eigen::Vector3d target_from_sun = target_sun_state.segment<3>(0);
    return shapiroRangeDelayForPositions(observer_from_sun, target_from_sun);
}

double SyntheticObservation::drawNoise() {
    return drawNoise(noise_.sigma);
}

double SyntheticObservation::drawNoise(double sigma) {
    if (!std::isfinite(sigma) || sigma < 0.0) {
        throw std::invalid_argument("Synthetic observation noise sigma must be finite and non-negative.");
    }
    if (!noise_.enabled || sigma == 0.0) {
        return 0.0;
    }

    std::normal_distribution<double> distribution(0.0, sigma);
    return distribution(rng_);
}

} // namespace fd::observations::synth
