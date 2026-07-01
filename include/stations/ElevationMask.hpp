#pragma once

#include <string>

namespace od {

struct TopocentricTargetGeometry {
    double rangeKm {0.0};
    double elevationRad {0.0};
    double lightTimeSec {0.0};
};

[[nodiscard]] TopocentricTargetGeometry computeTopocentricTargetGeometry(
    const std::string& target,
    const std::string& stationName,
    double epochTdb,
    const std::string& aberrationCorrection = "LT");

[[nodiscard]] double computeElevationRad(const std::string& target,
                                         const std::string& stationName,
                                         double epochTdb,
                                         const std::string& aberrationCorrection = "LT");

[[nodiscard]] double computeElevationDeg(const std::string& target,
                                         const std::string& stationName,
                                         double epochTdb,
                                         const std::string& aberrationCorrection = "LT");

[[nodiscard]] bool isVisibleAboveElevation(const std::string& target,
                                           const std::string& stationName,
                                           double epochTdb,
                                           double minimumElevationRad,
                                           const std::string& aberrationCorrection = "LT");

} // namespace od
