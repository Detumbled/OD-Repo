#include "StationCatalog.hpp"

#include <SpiceUsr.h>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

constexpr double kDegreeTolerance = 0.05;
constexpr double kAltitudeToleranceKm = 0.05;

void configureSpiceToReturn() {
    erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
    errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
}

void loadKernelPool(const char* metakernelPath) {
    furnsh_c(metakernelPath);

    if (failed_c()) {
        SpiceChar short_message[1841] = {0};
        SpiceChar long_message[1841] = {0};
        getmsg_c("SHORT", sizeof(short_message), short_message);
        getmsg_c("LONG", sizeof(long_message), long_message);
        reset_c();

        throw std::runtime_error(std::string("Failed to load meta-kernel: ")
                                 + short_message + " | " + long_message);
    }
}

double degreesToRadians(double degrees) {
    return degrees * rpd_c();
}

void validateDss43(const od::Station& station) {
    const double latitude_deg = station.latitudeRad() * dpr_c();
    const double longitude_deg = station.longitudeRad() * dpr_c();
    const double altitude_km = station.altitudeKm();

    const bool latitude_ok = std::abs(latitude_deg - (-35.402)) <= kDegreeTolerance;
    const bool longitude_ok = std::abs(longitude_deg - 148.981) <= kDegreeTolerance;
    const bool altitude_ok = std::abs(altitude_km - 0.689) <= kAltitudeToleranceKm;

    if (!latitude_ok || !longitude_ok || !altitude_ok) {
        throw std::runtime_error("DSS-43 validation failed: geodetic coordinates are outside tolerance.");
    }
}

void printStationSummary(const od::Station& station) {
    std::cout << std::fixed << std::setprecision(3)
              << station.name()
              << "  lat=" << station.latitudeRad() * dpr_c() << " deg"
              << "  lon=" << station.longitudeRad() * dpr_c() << " deg E"
              << "  alt=" << station.altitudeKm() << " km\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        configureSpiceToReturn();

        const char* metakernel = argc > 1 ? argv[1] : "../Kernels.tm";
        loadKernelPool(metakernel);

        SpiceDouble reference_et = 0.0;
        str2et_c("2026-06-21T00:00:00", &reference_et);
        if (failed_c()) {
            SpiceChar short_message[1841] = {0};
            SpiceChar long_message[1841] = {0};
            getmsg_c("SHORT", sizeof(short_message), short_message);
            getmsg_c("LONG", sizeof(long_message), long_message);
            reset_c();

            throw std::runtime_error(std::string("Failed to convert reference epoch: ")
                                     + short_message + " | " + long_message);
        }

        const std::vector<od::StationCatalogEntry> network = {
            {"DSS-14", od::stationNaifIdFromName("DSS-14")},
            {"DSS-43", od::stationNaifIdFromName("DSS-43")},
            {"DSS-63", od::stationNaifIdFromName("DSS-63")}
        };

        const auto stations = od::buildStationCatalogFromKernel(network, reference_et);

        for (const auto& station : stations) {
            printStationSummary(station);
        }

        const auto dss43 = od::buildStationFromKernel("DSS-43",
                                                      od::stationNaifIdFromName("DSS-43"),
                                                      reference_et,
                                                      {0.0, 0.0, 0.0},
                                                      {0.0, 0.0, 0.0},
                                                      degreesToRadians(10.0));
        validateDss43(dss43);

        std::cout << "DSS-43 validation passed.\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& exception) {
        std::cerr << exception.what() << '\n';
        return EXIT_FAILURE;
    }
}
