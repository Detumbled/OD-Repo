#ifndef STATION_CATALOG_HPP
#define STATION_CATALOG_HPP

#include "stations/Stations.hpp"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace od {

using StationCatalogEntry = std::pair<std::string, int>;
using StationBiasMap = std::unordered_map<std::string, Station::Vector3>;
using StationNoiseMap = std::unordered_map<std::string, Station::Vector3>;
using StationElevationMaskMap = std::unordered_map<std::string, double>;

[[nodiscard]] const std::vector<StationCatalogEntry>& defaultDsnStationCatalog();

[[nodiscard]] int stationNaifIdFromName(const std::string& commonName);

[[nodiscard]] Station buildStationFromKernel(
    const std::string& commonName,
    int naifId,
    double referenceEt,
    const Station::Vector3& biasKm = {0.0, 0.0, 0.0},
    const Station::Vector3& noiseSigmaKm = {0.0, 0.0, 0.0},
    double elevationMaskRad = 0.0);

[[nodiscard]] std::vector<Station> buildStationCatalogFromKernel(
    const std::vector<StationCatalogEntry>& stations,
    double referenceEt,
    const StationBiasMap& biasByName = {},
    const StationNoiseMap& noiseSigmaByName = {},
    const StationElevationMaskMap& elevationMaskByName = {});

[[nodiscard]] std::vector<Station> buildDefaultDsnCatalogFromKernel(
    double referenceEt,
    const StationBiasMap& biasByName = {},
    const StationNoiseMap& noiseSigmaByName = {},
    const StationElevationMaskMap& elevationMaskByName = {});

} // namespace od

#endif // STATION_CATALOG_HPP
