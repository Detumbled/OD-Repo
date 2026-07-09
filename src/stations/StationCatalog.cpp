#include "stations/StationCatalog.hpp"

#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <stdexcept>

namespace od {
namespace {

[[nodiscard]] const std::vector<StationCatalogEntry>& catalogStorage() {
    static const std::vector<StationCatalogEntry> kCatalog = {
        {"DSS-13", 399013},
        {"DSS-14", 399014},
        {"DSS-15", 399015},
        {"DSS-24", 399024},
        {"DSS-25", 399025},
        {"DSS-26", 399026},
        {"DSS-34", 399034},
        {"DSS-43", 399043},
        {"DSS-45", 399045},
        {"DSS-53", 399053},
        {"DSS-54", 399054},
        {"DSS-55", 399055},
        {"DSS-63", 399063},
        {"DSS-65", 399065}
    };

    return kCatalog;
}

[[nodiscard]] Station::Vector3 lookupVectorOrDefault(
    const std::unordered_map<std::string, Station::Vector3>& valuesByName,
    const std::string& stationName) {
    const auto iterator = valuesByName.find(stationName);
    return iterator == valuesByName.end() ? Station::Vector3{0.0, 0.0, 0.0} : iterator->second;
}

[[nodiscard]] double lookupMaskOrDefault(
    const StationElevationMaskMap& valuesByName,
    const std::string& stationName) {
    const auto iterator = valuesByName.find(stationName);
    return iterator == valuesByName.end() ? 0.0 : iterator->second;
}

template <typename Callable>
void withSpiceReturnMode(const std::string& context, Callable&& callable) {
    SpiceErrorModeGuard action_guard;
    callable();
    throwIfSpiceFailed(context);
}

} // namespace

const std::vector<StationCatalogEntry>& defaultDsnStationCatalog() {
    return catalogStorage();
}

int stationNaifIdFromName(const std::string& commonName) {
    for (const auto& [station_name, naif_id] : catalogStorage()) {
        if (station_name == commonName) {
            return naif_id;
        }
    }

    throw std::invalid_argument("Unknown station name: " + commonName);
}

Station buildStationFromKernel(const std::string& commonName,
                               int naifId,
                               double referenceEt,
                               const Station::Vector3& biasKm,
                               const Station::Vector3& noiseSigmaKm,
                               double elevationMaskRad) {
    Station::Vector3 ecef_position_km {0.0, 0.0, 0.0};
    double latitude_rad = 0.0;
    double longitude_rad = 0.0;
    double altitude_km = 0.0;

    const std::string naif_id_string = std::to_string(naifId);
    const std::string context = "Failed to build station " + commonName
        + " (NAIF ID " + naif_id_string + ")";

    withSpiceReturnMode(context, [&]() {
        SpiceDouble pos[3] = {0.0, 0.0, 0.0};
        SpiceDouble lt = 0.0;
        spkpos_c(naif_id_string.c_str(), referenceEt, "ITRF93", "NONE", "EARTH", pos, &lt);

        ecef_position_km = {pos[0], pos[1], pos[2]};

        SpiceInt dimension = 0;
        SpiceDouble radii[3] = {0.0, 0.0, 0.0};
        bodvrd_c("EARTH", "RADII", 3, &dimension, radii);

        const SpiceDouble equatorial_radius_km = radii[0];
        const SpiceDouble polar_radius_km = radii[2];
        const SpiceDouble flattening = (equatorial_radius_km - polar_radius_km) / equatorial_radius_km;

        SpiceDouble lon = 0.0;
        SpiceDouble lat = 0.0;
        SpiceDouble alt = 0.0;
        recgeo_c(pos, equatorial_radius_km, flattening, &lon, &lat, &alt);

        latitude_rad = lat;
        longitude_rad = lon;
        altitude_km = alt;
    });

    return Station(commonName,
                   ecef_position_km,
                   biasKm,
                   noiseSigmaKm,
                   latitude_rad,
                   longitude_rad,
                   altitude_km,
                   elevationMaskRad,
                   true);
}

std::vector<Station> buildStationCatalogFromKernel(
    const std::vector<StationCatalogEntry>& stations,
    double referenceEt,
    const StationBiasMap& biasByName,
    const StationNoiseMap& noiseSigmaByName,
    const StationElevationMaskMap& elevationMaskByName) {
    std::vector<Station> built_stations;
    built_stations.reserve(stations.size());

    for (const auto& [station_name, naif_id] : stations) {
        built_stations.push_back(buildStationFromKernel(station_name,
                                                        naif_id,
                                                        referenceEt,
                                                        lookupVectorOrDefault(biasByName, station_name),
                                                        lookupVectorOrDefault(noiseSigmaByName, station_name),
                                                        lookupMaskOrDefault(elevationMaskByName, station_name)));
    }

    return built_stations;
}

std::vector<Station> buildDefaultDsnCatalogFromKernel(
    double referenceEt,
    const StationBiasMap& biasByName,
    const StationNoiseMap& noiseSigmaByName,
    const StationElevationMaskMap& elevationMaskByName) {
    return buildStationCatalogFromKernel(catalogStorage(),
                                         referenceEt,
                                         biasByName,
                                         noiseSigmaByName,
                                         elevationMaskByName);
}

} // namespace od
