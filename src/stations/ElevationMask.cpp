#include "stations/ElevationMask.hpp"

#include "stations/StationCatalog.hpp"

#include <SpiceUsr.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace od {
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

[[nodiscard]] std::string stationObserverName(const std::string& stationName) {
    try {
        return std::to_string(stationNaifIdFromName(stationName));
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

void validateInputs(const std::string& target,
                    const std::string& stationName,
                    double epochTdb,
                    const std::string& aberrationCorrection) {
    if (target.empty()) {
        throw std::invalid_argument("Elevation-mask target cannot be empty.");
    }
    if (stationName.empty()) {
        throw std::invalid_argument("Elevation-mask station name cannot be empty.");
    }
    if (!std::isfinite(epochTdb)) {
        throw std::invalid_argument("Elevation-mask epoch must be finite.");
    }
    if (aberrationCorrection.empty()) {
        throw std::invalid_argument("Elevation-mask aberration correction cannot be empty.");
    }
}

} // namespace

TopocentricTargetGeometry computeTopocentricTargetGeometry(
    const std::string& target,
    const std::string& stationName,
    double epochTdb,
    const std::string& aberrationCorrection) {
    validateInputs(target, stationName, epochTdb, aberrationCorrection);

    SpiceErrorActionGuard action_guard;

    const std::string observer = stationObserverName(stationName);
    SpiceDouble station_position[3] = {0.0, 0.0, 0.0};
    SpiceDouble station_light_time = 0.0;
    SpiceDouble line_of_sight_itrf[3] = {0.0, 0.0, 0.0};
    SpiceDouble target_light_time = 0.0;

    spkpos_c(observer.c_str(),
             epochTdb,
             "ITRF93",
             "NONE",
             "EARTH",
             station_position,
             &station_light_time);

    SpiceInt dimension = 0;
    SpiceDouble radii[3] = {0.0, 0.0, 0.0};
    bodvrd_c("EARTH", "RADII", 3, &dimension, radii);

    const SpiceDouble equatorial_radius_km = radii[0];
    const SpiceDouble polar_radius_km = radii[2];
    const SpiceDouble flattening = (equatorial_radius_km - polar_radius_km) / equatorial_radius_km;

    SpiceDouble longitude_rad = 0.0;
    SpiceDouble latitude_rad = 0.0;
    SpiceDouble altitude_km = 0.0;
    recgeo_c(station_position,
             equatorial_radius_km,
             flattening,
             &longitude_rad,
             &latitude_rad,
             &altitude_km);

    spkpos_c(target.c_str(),
             epochTdb,
             "ITRF93",
             aberrationCorrection.c_str(),
             observer.c_str(),
             line_of_sight_itrf,
             &target_light_time);

    throwIfSpiceFailed("Failed to compute topocentric target geometry for " + stationName);

    const double range = std::sqrt(line_of_sight_itrf[0] * line_of_sight_itrf[0]
                                   + line_of_sight_itrf[1] * line_of_sight_itrf[1]
                                   + line_of_sight_itrf[2] * line_of_sight_itrf[2]);
    if (!(range > 0.0) || !std::isfinite(range)) {
        throw std::runtime_error("Topocentric target geometry has invalid range.");
    }

    const double cos_lat = std::cos(latitude_rad);
    const double up_x = cos_lat * std::cos(longitude_rad);
    const double up_y = cos_lat * std::sin(longitude_rad);
    const double up_z = std::sin(latitude_rad);
    const double sine_elevation = std::clamp((line_of_sight_itrf[0] * up_x
                                              + line_of_sight_itrf[1] * up_y
                                              + line_of_sight_itrf[2] * up_z) / range,
                                             -1.0,
                                             1.0);
    return TopocentricTargetGeometry{
        range,
        std::asin(sine_elevation),
        target_light_time
    };
}

double computeElevationRad(const std::string& target,
                           const std::string& stationName,
                           double epochTdb,
                           const std::string& aberrationCorrection) {
    return computeTopocentricTargetGeometry(target,
                                            stationName,
                                            epochTdb,
                                            aberrationCorrection).elevationRad;
}

double computeElevationDeg(const std::string& target,
                           const std::string& stationName,
                           double epochTdb,
                           const std::string& aberrationCorrection) {
    return computeElevationRad(target, stationName, epochTdb, aberrationCorrection) * dpr_c();
}

bool isVisibleAboveElevation(const std::string& target,
                             const std::string& stationName,
                             double epochTdb,
                             double minimumElevationRad,
                             const std::string& aberrationCorrection) {
    if (!std::isfinite(minimumElevationRad)) {
        throw std::invalid_argument("Elevation mask must be finite.");
    }

    return computeElevationRad(target, stationName, epochTdb, aberrationCorrection) >= minimumElevationRad;
}

} // namespace od
