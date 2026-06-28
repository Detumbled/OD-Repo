#include "StationCatalog.hpp"
#include "observations/synth/RangeRateSynth.hpp"
#include "observations/synth/RangeSynth.hpp"

#include <SpiceUsr.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-03-05T00:00:00";
constexpr std::array<const char*, 2> kStationNames {"DSS-43", "DSS-63"};
constexpr const char* kTarget = "-31";
constexpr double kArcDurationSec = 3.0 * 3600.0;
constexpr double kCadenceSec = 180.0;
constexpr double kRangeSigmaKm = 0.010;
constexpr double kRangeRateSigmaKmPerSec = 1.0e-6;
constexpr const char* kReportPath = "../synthetic_observations_dss43_voyager1.txt";

struct StationObservationSet {
    od::Station station;
    fd::observations::synth::GeometryConfig geometry;
    std::vector<fd::observations::synth::SyntheticObservationSample> ranges;
    std::vector<fd::observations::synth::SyntheticObservationSample> rangeRates;
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
    const std::vector<StationObservationSet>& stationSets) {
    if (stationSets.empty()) {
        throw std::runtime_error("Synthetic observation report requires at least one station.");
    }

    const std::size_t sample_count = stationSets.front().ranges.size();
    if (sample_count == 0) {
        throw std::runtime_error("Synthetic observation report cannot be empty.");
    }

    for (const StationObservationSet& station_set : stationSets) {
        if (station_set.ranges.size() != sample_count || station_set.rangeRates.size() != sample_count) {
            throw std::runtime_error("Synthetic station sample counts do not match.");
        }
    }

    std::ofstream report(path);
    if (!report) {
        throw std::runtime_error("Failed to open synthetic observation report: " + path);
    }

    report << std::fixed << std::setprecision(9);
    report << "# Synthetic Voyager 1 observations\n"
           << "# target_naif_id        : " << kTarget << '\n'
           << "# station_count         : " << stationSets.size() << '\n';

    for (const StationObservationSet& station_set : stationSets) {
        const od::Station& station = station_set.station;
        report << "# station               : " << station.name()
               << " lat_deg=" << station.latitudeRad() * dpr_c()
               << " lon_deg_east=" << station.longitudeRad() * dpr_c()
               << " alt_km=" << station.altitudeKm() << '\n';
    }

    report << "# frame                 : " << stationSets.front().geometry.frame << '\n'
           << "# aberration_correction : " << stationSets.front().geometry.aberrationCorrection << '\n'
           << "# relativistic_delay    : Solar Shapiro range delay applied\n"
           << "# range_sigma_km        : " << kRangeSigmaKm << '\n'
           << "# range_rate_sigma_km_s : " << kRangeRateSigmaKmPerSec << '\n'
           << "# cadence_seconds       : " << kCadenceSec << '\n'
           << "# columns: utc station epoch_tdb range_truth_shapiro_km range_noise_km range_observed_km "
              "range_sigma_km range_rate_truth_shapiro_km_s range_rate_noise_km_s "
              "range_rate_observed_km_s range_rate_sigma_km_s\n";

    for (std::size_t i = 0; i < sample_count; ++i) {
        for (const StationObservationSet& station_set : stationSets) {
            const auto& range = station_set.ranges[i];
            const auto& range_rate = station_set.rangeRates[i];
            report << utcFromEt(range.epochTdb) << ' '
                   << station_set.station.name() << ' '
                   << range.epochTdb << ' '
                   << range.truth << ' '
                   << range.noise << ' '
                   << range.observed << ' '
                   << range.sigma << ' '
                   << range_rate.truth << ' '
                   << range_rate.noise << ' '
                   << range_rate.observed << ' '
                   << range_rate.sigma << '\n';
        }
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

        const double end_et = start_et + kArcDurationSec;
        std::vector<StationObservationSet> station_sets;
        station_sets.reserve(kStationNames.size());

        for (std::size_t station_index = 0; station_index < kStationNames.size(); ++station_index) {
            const std::string station_name = kStationNames[station_index];
            const int station_naif_id = od::stationNaifIdFromName(station_name);
            od::Station station = od::buildStationFromKernel(station_name, station_naif_id, start_et);

            fd::observations::synth::GeometryConfig geometry;
            geometry.target = kTarget;
            geometry.stationName = station_name;
            geometry.frame = "J2000";
            geometry.aberrationCorrection = "LT";

            fd::observations::synth::NoiseConfig range_noise;
            range_noise.enabled = true;
            range_noise.sigma = kRangeSigmaKm;
            range_noise.seed = 4301U + static_cast<std::uint32_t>(100U * station_index);

            fd::observations::synth::NoiseConfig range_rate_noise;
            range_rate_noise.enabled = true;
            range_rate_noise.sigma = kRangeRateSigmaKmPerSec;
            range_rate_noise.seed = 4302U + static_cast<std::uint32_t>(100U * station_index);

            fd::observations::synth::RangeSynth range_synth(geometry, range_noise);
            fd::observations::synth::RangeRateSynth range_rate_synth(geometry, range_rate_noise);

            station_sets.push_back(StationObservationSet{
                std::move(station),
                std::move(geometry),
                range_synth.generate(start_et, end_et, kCadenceSec),
                range_rate_synth.generate(start_et, end_et, kCadenceSec)
            });
        }

        writeReport(kReportPath, station_sets);

        const std::size_t sample_count = station_sets.front().ranges.size();
        const std::size_t total_observation_count = sample_count * station_sets.size();

        std::cout << "Generated " << total_observation_count
                  << " synthetic DSS-43/DSS-63 Voyager 1 observations in " << kReportPath << '\n'
                  << "First UTC: " << utcFromEt(station_sets.front().ranges.front().epochTdb) << '\n'
                  << "Last UTC : " << utcFromEt(station_sets.front().ranges.back().epochTdb) << '\n'
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
