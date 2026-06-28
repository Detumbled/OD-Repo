#include "StationCatalog.hpp"
#include "observations/synth/RangeRateSynth.hpp"
#include "observations/synth/RangeSynth.hpp"

#include <SpiceUsr.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-03-05T00:00:00";
constexpr const char* kStationName = "DSS-43";
constexpr const char* kTarget = "-31";
constexpr double kArcDurationSec = 3600.0;
constexpr double kCadenceSec = 300.0;
constexpr double kRangeSigmaKm = 0.010;
constexpr double kRangeRateSigmaKmPerSec = 1.0e-6;
constexpr const char* kReportPath = "../synthetic_observations_dss43_voyager1.txt";

void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar short_message[1841] = {0};
    SpiceChar long_message[1841] = {0};
    getmsg_c("SHORT", sizeof(short_message), short_message);
    getmsg_c("LONG", sizeof(long_message), long_message);
    reset_c();

    throw std::runtime_error(context + ": " + short_message + " | " + long_message);
}

[[nodiscard]] std::string utcFromEt(double et) {
    SpiceChar utc[128] = {0};
    et2utc_c(et, "ISOC", 3, static_cast<SpiceInt>(sizeof(utc)), utc);
    throwIfSpiceFailed("Failed to convert ET to UTC");
    return utc;
}

void writeReport(
    const std::string& path,
    const od::Station& station,
    const fd::observations::synth::GeometryConfig& geometry,
    const std::vector<fd::observations::synth::SyntheticObservationSample>& ranges,
    const std::vector<fd::observations::synth::SyntheticObservationSample>& range_rates) {
    if (ranges.size() != range_rates.size()) {
        throw std::runtime_error("Synthetic range and range-rate sample counts do not match.");
    }

    std::ofstream report(path);
    if (!report) {
        throw std::runtime_error("Failed to open synthetic observation report: " + path);
    }

    report << std::fixed << std::setprecision(9);
    report << "# Synthetic Voyager 1 observations\n"
           << "# target_naif_id        : " << kTarget << '\n'
           << "# station               : " << station.name() << '\n'
           << "# station_lat_deg       : " << station.latitudeRad() * dpr_c() << '\n'
           << "# station_lon_deg_east  : " << station.longitudeRad() * dpr_c() << '\n'
           << "# station_alt_km        : " << station.altitudeKm() << '\n'
           << "# frame                 : " << geometry.frame << '\n'
           << "# aberration_correction : " << geometry.aberrationCorrection << '\n'
           << "# relativistic_delay    : Solar Shapiro range delay applied\n"
           << "# range_sigma_km        : " << kRangeSigmaKm << '\n'
           << "# range_rate_sigma_km_s : " << kRangeRateSigmaKmPerSec << '\n'
           << "# cadence_seconds       : " << kCadenceSec << '\n'
           << "# columns: utc epoch_tdb range_truth_shapiro_km range_noise_km range_observed_km "
              "range_sigma_km range_rate_truth_shapiro_km_s range_rate_noise_km_s "
              "range_rate_observed_km_s range_rate_sigma_km_s\n";

    for (std::size_t i = 0; i < ranges.size(); ++i) {
        report << utcFromEt(ranges[i].epochTdb) << ' '
               << ranges[i].epochTdb << ' '
               << ranges[i].truth << ' '
               << ranges[i].noise << ' '
               << ranges[i].observed << ' '
               << ranges[i].sigma << ' '
               << range_rates[i].truth << ' '
               << range_rates[i].noise << ' '
               << range_rates[i].observed << ' '
               << range_rates[i].sigma << '\n';
    }
}

} // namespace

int main() {
    try {
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
        errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load synthetic observation meta-kernel");

        SpiceDouble start_et = 0.0;
        str2et_c(kStartUtc, &start_et);
        throwIfSpiceFailed("Failed to convert synthetic observation start epoch");

        const int station_naif_id = od::stationNaifIdFromName(kStationName);
        const od::Station station = od::buildStationFromKernel(kStationName, station_naif_id, start_et);

        fd::observations::synth::GeometryConfig geometry;
        geometry.target = kTarget;
        geometry.stationName = kStationName;
        geometry.frame = "J2000";
        geometry.aberrationCorrection = "LT";

        fd::observations::synth::NoiseConfig range_noise;
        range_noise.enabled = true;
        range_noise.sigma = kRangeSigmaKm;
        range_noise.seed = 4301U;

        fd::observations::synth::NoiseConfig range_rate_noise;
        range_rate_noise.enabled = true;
        range_rate_noise.sigma = kRangeRateSigmaKmPerSec;
        range_rate_noise.seed = 4302U;

        fd::observations::synth::RangeSynth range_synth(geometry, range_noise);
        fd::observations::synth::RangeRateSynth range_rate_synth(geometry, range_rate_noise);

        const double end_et = start_et + kArcDurationSec;
        const auto ranges = range_synth.generate(start_et, end_et, kCadenceSec);
        const auto range_rates = range_rate_synth.generate(start_et, end_et, kCadenceSec);

        writeReport(kReportPath, station, geometry, ranges, range_rates);

        std::cout << "Generated " << ranges.size()
                  << " synthetic DSS-43 Voyager 1 observations in " << kReportPath << '\n'
                  << "First UTC: " << utcFromEt(ranges.front().epochTdb) << '\n'
                  << "Last UTC : " << utcFromEt(ranges.back().epochTdb) << '\n'
                  << "Range sigma: " << kRangeSigmaKm * 1000.0 << " m\n"
                  << "Range-rate sigma: " << kRangeRateSigmaKmPerSec * 1000.0 << " m/s\n";

        kclear_c();
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Synthetic observation test failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
