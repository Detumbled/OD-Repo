#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "stations/StationCatalog.hpp"
#include "observations/synth/RangeRateSynth.hpp"
#include "observations/synth/RangeSynth.hpp"
#include "observations/synth/VLBISynth.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "perturbations/Shapiro.hpp"
#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using State6 = Eigen::Matrix<double, 6, 1>;
using od::throwIfSpiceFailed;

struct VLBIBaselineConfig {
    const char* stationOneName;
    const char* stationTwoName;
    std::uint32_t noiseSeed;
};

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-01-10T00:00:00";
constexpr std::array<const char*, 2> kStationNames {"DSS-43", "DSS-63"};
constexpr std::array<VLBIBaselineConfig, 2> kVLBIBaselines {{
    {"DSS-43", "DSS-63", 5301U},
    {"DSS-14", "DSS-43", 5401U}
}};
constexpr const char* kTarget = "-31";
constexpr const char* kCentralBody = "SUN";
constexpr const char* kFrame = "J2000";
constexpr const char* kJupiter = "599";
constexpr const char* kIo = "501";
constexpr const char* kEuropa = "502";
constexpr const char* kGanymede = "503";
constexpr const char* kCallisto = "504";
constexpr const char* kSaturnBarycenter = "6";
constexpr double kArcDurationSec = 2 * 48.0 * 3600.0;
constexpr double kCadenceSec = 180.0;
constexpr double kVLBICadenceSec = 20.0;
constexpr std::size_t kMinVLBISamplesPerBaseline = 20;
constexpr double kElevationMaskDeg = 10.0;
constexpr double kRangeSigmaKm = 0.010;
constexpr double kRangeRateSigmaKmPerSec = 1.0e-6;
constexpr double kVLBISigmaKm = 0.001;
constexpr double kSunMuKm3PerSec2 = fd::perturbations::kSunGravitationalParameterKm3PerSec2;
constexpr double kJupiterBodyMuKm3PerSec2 = 126686531.9003704;
constexpr double kIoMuKm3PerSec2 = 5959.916;
constexpr double kEuropaMuKm3PerSec2 = 3202.739;
constexpr double kGanymedeMuKm3PerSec2 = 9887.834;
constexpr double kCallistoMuKm3PerSec2 = 7179.289;
constexpr double kSaturnSystemMuKm3PerSec2 = 37940585.2;
constexpr double kVoyagerReflectivity = 1.3;
constexpr double kVoyagerAreaM2 = 10.0;
constexpr double kVoyagerMassKg = 721.9;
constexpr double kRangeRateCountTimeSec = 60.0;
constexpr double kSyntheticEphemerisLookbackSec = 3.0 * 3600.0;
constexpr const char* kReportPath = "../synthetic_observations_dss43_voyager1.txt";

struct StationObservationSet {
    od::Station station;
    std::string stationNaif;
    fd::observations::synth::GeometryConfig geometry;
    std::vector<fd::observations::synth::SyntheticObservationSample> ranges;
    std::vector<fd::observations::synth::SyntheticObservationSample> rangeRates;
};

struct VLBIObservationSet {
    std::string stationOneName;
    std::string stationTwoName;
    std::string stationOneNaif;
    std::string stationTwoNaif;
    std::vector<fd::observations::synth::SyntheticObservationSample> delays;
};

struct StationInfo {
    od::Station station;
    std::string stationNaif;
};

struct ReportRow {
    const StationObservationSet* stationSet {nullptr};
    std::size_t sampleIndex {0};
};

struct ComparisonStats {
    std::size_t count {0};
    double sumSquared {0.0};
    double maxAbs {0.0};

    void add(double delta) {
        ++count;
        sumSquared += delta * delta;
        maxAbs = std::max(maxAbs, std::abs(delta));
    }

    [[nodiscard]] double rms() const {
        return count == 0U ? 0.0 : std::sqrt(sumSquared / static_cast<double>(count));
    }
};

[[nodiscard]] State6 spiceState(const std::string& target,
                                double tdb,
                                const std::string& observer,
                                const std::string& aberrationCorrection = "LT") {
    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkezr_c(target.c_str(),
             tdb,
             kFrame,
             aberrationCorrection.c_str(),
             observer.c_str(),
             spice_state,
             &light_time);
    throwIfSpiceFailed("Failed to fetch SPICE state for " + target + " relative to " + observer);

    State6 state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spice_state[i];
    }
    return state;
}

[[nodiscard]] bool spiceStateAvailable(const std::string& target,
                                       double tdb,
                                       const std::string& observer) {
    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkezr_c(target.c_str(), tdb, kFrame, "NONE", observer.c_str(), spice_state, &light_time);
    if (failed_c()) {
        reset_c();
        return false;
    }
    return true;
}

void addThirdBodyIfCovered(fd::perturbations::ThirdBodyGravity& thirdBody,
                           const fd::perturbations::MassiveBody& body,
                           double startEpoch,
                           double endEpoch,
                           bool required,
                           std::vector<std::string>& activeBodies,
                           std::vector<std::string>& skippedBodies) {
    const bool has_coverage = spiceStateAvailable(body.name, startEpoch, kCentralBody)
        && spiceStateAvailable(body.name, endEpoch, kCentralBody);

    if (has_coverage) {
        thirdBody.addBody(body);
        activeBodies.push_back(body.name);
        return;
    }

    if (required) {
        throw std::runtime_error("Required third-body SPK coverage is missing for " + body.name);
    }

    skippedBodies.push_back(body.name);
}

[[nodiscard]] od::RKF45Integrator::DerivativeFunction makeDynamics(
    const fd::perturbations::ThirdBodyGravity& thirdBody,
    const fd::perturbations::SolarRadiationPressure& srp) {
    return [&thirdBody, &srp](double tdb,
                              od::RKF45Integrator::ConstStateRef state,
                              od::RKF45Integrator::StateRef derivative) {
        const Eigen::Vector3d position = state.segment<3>(0);
        const Eigen::Vector3d velocity = state.segment<3>(3);
        const double radius = position.norm();
        if (radius <= 0.0) {
            throw std::runtime_error("Dynamics received zero Sun-spacecraft radius.");
        }

        derivative.segment<3>(0) = velocity;
        derivative.segment<3>(3) =
            -kSunMuKm3PerSec2 * position / (radius * radius * radius)
            + thirdBody.computeAcceleration(tdb, position)
            + srp.computeAcceleration(tdb, position);
    };
}

void appendHistory(od::EphemerisInterpolator& ephemeris,
                   const od::RKF45Integrator::Result& result) {
    if (result.history.empty()) {
        throw std::runtime_error("RKF45 propagation returned an empty ephemeris history.");
    }

    for (const od::RKF45Integrator::Result::EphemerisNode& node : result.history) {
        ephemeris.addNode(node.tdb, node.state, node.derivative);
    }
}

[[nodiscard]] State6 asState6(const Eigen::VectorXd& state, const std::string& context) {
    if (state.size() != 6 || !state.allFinite()) {
        throw std::runtime_error(context + " produced an invalid Voyager state.");
    }

    State6 fixed_state;
    fixed_state = state;
    return fixed_state;
}

[[nodiscard]] od::EphemerisInterpolator propagateEphemerisSpan(
    const od::RKF45Integrator& integrator,
    const od::RKF45Integrator::DerivativeFunction& dynamics,
    double referenceEpoch,
    const State6& referenceState,
    double startEpoch,
    double finalEpoch) {
    if (!std::isfinite(startEpoch) || !std::isfinite(finalEpoch) || startEpoch > finalEpoch) {
        throw std::invalid_argument("Ephemeris span bounds must be finite and ordered.");
    }

    od::EphemerisInterpolator ephemeris;
    if (startEpoch < referenceEpoch) {
        const od::RKF45Integrator::State dynamic_initial = referenceState;
        const auto backward = integrator.integrate(referenceEpoch, dynamic_initial, startEpoch, dynamics);
        (void) asState6(backward.state, "Backward propagation");
        appendHistory(ephemeris, backward);
    }

    if (finalEpoch > referenceEpoch) {
        const od::RKF45Integrator::State dynamic_initial = referenceState;
        const auto forward = integrator.integrate(referenceEpoch, dynamic_initial, finalEpoch, dynamics);
        (void) asState6(forward.state, "Forward propagation");
        appendHistory(ephemeris, forward);
    }

    if (ephemeris.history().empty()) {
        const od::RKF45Integrator::State dynamic_initial = referenceState;
        const auto result = integrator.integrate(referenceEpoch, dynamic_initial, referenceEpoch, dynamics);
        appendHistory(ephemeris, result);
    }

    return ephemeris;
}

[[nodiscard]] double rangeRateSigmaForCountTime() {
    return kRangeRateSigmaKmPerSec / std::sqrt(kRangeRateCountTimeSec);
}

[[nodiscard]] std::string utcFromEt(double et) {
    SpiceChar utc[128] = {0};
    et2utc_c(et, "ISOC", 3, static_cast<SpiceInt>(sizeof(utc)), utc);
    throwIfSpiceFailed("Failed to convert ET to UTC");
    return utc;
}

[[nodiscard]] StationInfo stationInfoForName(const std::string& stationName, double epochTdb) {
    const int station_naif_id = od::stationNaifIdFromName(stationName);
    return StationInfo{
        od::buildStationFromKernel(stationName, station_naif_id, epochTdb),
        std::to_string(station_naif_id)
    };
}

[[nodiscard]] fd::observations::synth::NoiseConfig makeNoiseConfig(double sigma,
                                                                   std::uint32_t seed,
                                                                   bool enabled = true) {
    fd::observations::synth::NoiseConfig noise;
    noise.enabled = enabled;
    noise.sigma = sigma;
    noise.seed = seed;
    return noise;
}

[[nodiscard]] const fd::observations::synth::SyntheticObservationSample& singleSample(
    const std::vector<fd::observations::synth::SyntheticObservationSample>& samples,
    const std::string& context) {
    if (samples.size() != 1U) {
        throw std::runtime_error(context + " did not generate exactly one sample.");
    }
    return samples.front();
}

[[nodiscard]] VLBIObservationSet generateVlbiObservationSet(
    const VLBIBaselineConfig& baseline,
    const fd::observations::synth::TargetStateProvider& targetProvider,
    double startEpoch,
    double endEpoch,
    double stationGeometryEpoch) {
    const StationInfo station_one = stationInfoForName(baseline.stationOneName, stationGeometryEpoch);
    const StationInfo station_two = stationInfoForName(baseline.stationTwoName, stationGeometryEpoch);

    VLBIObservationSet vlbi_set;
    vlbi_set.stationOneName = station_one.station.name();
    vlbi_set.stationTwoName = station_two.station.name();
    vlbi_set.stationOneNaif = station_one.stationNaif;
    vlbi_set.stationTwoNaif = station_two.stationNaif;

    fd::observations::synth::VLBIConfig vlbi_config;
    vlbi_config.target = kTarget;
    vlbi_config.stationOneName = vlbi_set.stationOneName;
    vlbi_config.stationTwoName = vlbi_set.stationTwoName;
    vlbi_config.frame = kFrame;
    vlbi_config.aberrationCorrection = "LT";
    vlbi_config.minimumElevationRad = kElevationMaskDeg / dpr_c();

    fd::observations::synth::VLBISynth vlbi_synth(
        vlbi_config,
        makeNoiseConfig(kVLBISigmaKm, baseline.noiseSeed));
    vlbi_set.delays = vlbi_synth.generate(startEpoch,
                                          endEpoch,
                                          kVLBICadenceSec,
                                          targetProvider);

    if (vlbi_set.delays.size() < kMinVLBISamplesPerBaseline) {
        std::ostringstream message;
        message << "Elevation mask produced only " << vlbi_set.delays.size()
                << " synthetic VLBI observations for "
                << vlbi_set.stationOneName << '/' << vlbi_set.stationTwoName
                << "; need at least " << kMinVLBISamplesPerBaseline << '.';
        throw std::runtime_error(message.str());
    }

    return vlbi_set;
}

[[nodiscard]] std::vector<ReportRow> makeSortedReportRows(
    const std::vector<StationObservationSet>& stationSets) {
    std::vector<ReportRow> rows;
    for (const StationObservationSet& station_set : stationSets) {
        if (station_set.ranges.size() != station_set.rangeRates.size()) {
            throw std::runtime_error("Synthetic station range and range-rate sample counts do not match.");
        }

        rows.reserve(rows.size() + station_set.ranges.size());
        for (std::size_t i = 0; i < station_set.ranges.size(); ++i) {
            rows.push_back(ReportRow{&station_set, i});
        }
    }

    std::sort(rows.begin(),
              rows.end(),
              [](const ReportRow& lhs, const ReportRow& rhs) {
                  const auto& lhs_range = lhs.stationSet->ranges[lhs.sampleIndex];
                  const auto& rhs_range = rhs.stationSet->ranges[rhs.sampleIndex];
                  if (lhs_range.epochTdb != rhs_range.epochTdb) {
                      return lhs_range.epochTdb < rhs_range.epochTdb;
                  }
                  return lhs.stationSet->station.name() < rhs.stationSet->station.name();
              });

    return rows;
}

void writeReport(
    const std::string& path,
    const std::vector<StationObservationSet>& stationSets,
    const std::vector<VLBIObservationSet>& vlbiSets) {
    if (stationSets.empty()) {
        throw std::runtime_error("Synthetic observation report requires at least one station.");
    }

    const std::vector<ReportRow> rows = makeSortedReportRows(stationSets);
    if (rows.empty()) {
        throw std::runtime_error("Synthetic observation report cannot be empty.");
    }

    for (const StationObservationSet& station_set : stationSets) {
        if (station_set.ranges.empty()) {
            throw std::runtime_error("Synthetic observation report cannot include an empty station set.");
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
               << " alt_km=" << station.altitudeKm()
               << " visible_samples=" << station_set.ranges.size() << '\n';
    }
    for (const VLBIObservationSet& vlbi_set : vlbiSets) {
        report << "# vlbi_pair             : " << vlbi_set.stationOneName
               << ' ' << vlbi_set.stationTwoName
               << " visible_samples=" << vlbi_set.delays.size() << '\n';
    }

    report << "# frame                 : " << stationSets.front().geometry.frame << '\n'
           << "# aberration_correction : " << stationSets.front().geometry.aberrationCorrection << '\n'
           << "# relativistic_delay    : Solar Shapiro range delay applied\n"
           << "# elevation_mask_deg    : " << kElevationMaskDeg << '\n'
           << "# range_sigma_km        : " << kRangeSigmaKm << '\n'
           << "# range_rate_base_sigma_1s_km_s : " << kRangeRateSigmaKmPerSec << '\n'
           << "# range_rate_count_time_s : " << kRangeRateCountTimeSec << '\n'
           << "# range_rate_sigma_km_s : " << rangeRateSigmaForCountTime() << '\n'
           << "# vlbi_sigma_km         : " << kVLBISigmaKm << '\n'
           << "# cadence_seconds       : " << kCadenceSec << '\n'
           << "# vlbi_cadence_seconds  : " << kVLBICadenceSec << '\n'
           << "# vlbi_min_samples_per_pair : " << kMinVLBISamplesPerBaseline << '\n'
           << "# vlbi_model            : differential one-way range, station2_minus_station1_km\n"
           << "# columns: utc station epoch_tdb range_truth_shapiro_km range_noise_km range_observed_km "
              "range_sigma_km range_rate_truth_shapiro_km_s range_rate_noise_km_s "
              "range_rate_observed_km_s range_rate_sigma_km_s\n"
           << "# vlbi_columns: utc VLBI station1 station2 epoch_tdb vlbi_truth_km vlbi_noise_km "
              "vlbi_observed_km vlbi_sigma_km\n";

    for (const ReportRow& row : rows) {
        const StationObservationSet& station_set = *row.stationSet;
        const auto& range = station_set.ranges[row.sampleIndex];
        const auto& range_rate = station_set.rangeRates[row.sampleIndex];
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

    for (const VLBIObservationSet& vlbi_set : vlbiSets) {
        for (const auto& delay : vlbi_set.delays) {
            report << utcFromEt(delay.epochTdb) << " VLBI "
                   << vlbi_set.stationOneName << ' '
                   << vlbi_set.stationTwoName << ' '
                   << delay.epochTdb << ' '
                   << delay.truth << ' '
                   << delay.noise << ' '
                   << delay.observed << ' '
                   << delay.sigma << '\n';
        }
    }
}

void addCspiceComparisonForStation(const StationObservationSet& stationSet,
                                   ComparisonStats& rangeStats,
                                   ComparisonStats& rangeRateStats) {
    const fd::observations::synth::NoiseConfig no_noise = makeNoiseConfig(0.0, 0U, false);
    fd::observations::synth::RangeSynth cspice_range_synth(stationSet.geometry, no_noise);
    fd::observations::synth::RangeRateSynth cspice_rate_synth(stationSet.geometry, no_noise);

    for (std::size_t i = 0; i < stationSet.ranges.size(); ++i) {
        const auto& simulated_range = stationSet.ranges[i];
        const auto& simulated_rate = stationSet.rangeRates[i];
        const auto cspice_range =
            cspice_range_synth.generate(simulated_range.epochTdb,
                                        simulated_range.epochTdb,
                                        kCadenceSec);
        const auto cspice_rate =
            cspice_rate_synth.generate(simulated_rate.epochTdb,
                                       simulated_rate.epochTdb,
                                       kCadenceSec,
                                       kRangeRateCountTimeSec);

        rangeStats.add(simulated_range.truth
                       - singleSample(cspice_range, "CSPICE range comparison").truth);
        rangeRateStats.add(simulated_rate.truth
                           - singleSample(cspice_rate, "CSPICE range-rate comparison").truth);
    }
}

void addCspiceComparisonForVlbi(const VLBIObservationSet& vlbiSet,
                                ComparisonStats& vlbiStats) {
    const fd::observations::synth::NoiseConfig no_noise = makeNoiseConfig(0.0, 0U, false);
    fd::observations::synth::VLBIConfig vlbi_config;
    vlbi_config.target = kTarget;
    vlbi_config.stationOneName = vlbiSet.stationOneName;
    vlbi_config.stationTwoName = vlbiSet.stationTwoName;
    vlbi_config.frame = kFrame;
    vlbi_config.aberrationCorrection = "LT";
    vlbi_config.minimumElevationRad = -1.57079632679489661923;

    fd::observations::synth::VLBISynth cspice_vlbi_synth(vlbi_config, no_noise);
    for (const auto& simulated_delay : vlbiSet.delays) {
        const auto cspice_vlbi = cspice_vlbi_synth.generate(simulated_delay.epochTdb,
                                                            simulated_delay.epochTdb,
                                                            kVLBICadenceSec);
        vlbiStats.add(simulated_delay.truth
                      - singleSample(cspice_vlbi, "CSPICE VLBI comparison").truth);
    }
}

} // namespace

int main() {
    try {
        const od::SpiceErrorModeGuard spice_error_mode("RETURN", "NONE");

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load synthetic observation meta-kernel");

        SpiceDouble start_et = 0.0;
        str2et_c(kStartUtc, &start_et);
        throwIfSpiceFailed("Failed to convert synthetic observation start epoch");

        const double end_et = start_et + kArcDurationSec;
        const double ephemeris_start_et = start_et - kSyntheticEphemerisLookbackSec;
        const double ephemeris_end_et = end_et + 0.5 * kRangeRateCountTimeSec;

        od::RKF45Integrator::Options options;
        options.absoluteTolerance = 1.0e-5;
        options.relativeTolerance = 1.0e-12;
        options.initialStep = 30.0;
        options.minimumStep = 1.0e-4;
        options.maximumStep = 300.0;
        options.maximumSteps = 200000;
        const od::RKF45Integrator integrator(options);

        fd::perturbations::ThirdBodyGravity third_body(kCentralBody, kFrame);
        std::vector<std::string> active_third_bodies;
        std::vector<std::string> skipped_third_bodies;
        addThirdBodyIfCovered(third_body,
                              {kJupiter, kJupiterBodyMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              true,
                              active_third_bodies,
                              skipped_third_bodies);
        addThirdBodyIfCovered(third_body,
                              {kIo, kIoMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              false,
                              active_third_bodies,
                              skipped_third_bodies);
        addThirdBodyIfCovered(third_body,
                              {kEuropa, kEuropaMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              false,
                              active_third_bodies,
                              skipped_third_bodies);
        addThirdBodyIfCovered(third_body,
                              {kGanymede, kGanymedeMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              false,
                              active_third_bodies,
                              skipped_third_bodies);
        addThirdBodyIfCovered(third_body,
                              {kCallisto, kCallistoMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              false,
                              active_third_bodies,
                              skipped_third_bodies);
        addThirdBodyIfCovered(third_body,
                              {kSaturnBarycenter, kSaturnSystemMuKm3PerSec2},
                              ephemeris_start_et,
                              ephemeris_end_et,
                              true,
                              active_third_bodies,
                              skipped_third_bodies);

        fd::perturbations::SolarRadiationPressure srp(kCentralBody,
                                                      kFrame,
                                                      kVoyagerReflectivity,
                                                      kVoyagerAreaM2,
                                                      kVoyagerMassKg);
        const auto dynamics = makeDynamics(third_body, srp);
        const State6 truth_initial_state = spiceState(kTarget, start_et, kCentralBody, "NONE");
        const od::EphemerisInterpolator truth_ephemeris =
            propagateEphemerisSpan(integrator,
                                   dynamics,
                                   start_et,
                                   truth_initial_state,
                                   ephemeris_start_et,
                                   ephemeris_end_et);
        const fd::observations::synth::InterpolatedTargetStateProvider truth_provider(
            truth_ephemeris,
            "Synthetic truth ephemeris interpolation");

        const std::vector<double> epochs =
            fd::observations::synth::SyntheticObservation::makeEpochGrid(start_et, end_et, kCadenceSec);

        std::vector<StationObservationSet> station_sets;
        station_sets.reserve(kStationNames.size());

        for (std::size_t station_index = 0; station_index < kStationNames.size(); ++station_index) {
            const std::string station_name = kStationNames[station_index];
            StationInfo station_info = stationInfoForName(station_name, start_et);

            fd::observations::synth::GeometryConfig geometry;
            geometry.target = kTarget;
            geometry.stationName = station_name;
            geometry.frame = "J2000";
            geometry.aberrationCorrection = "LT";

            fd::observations::synth::RangeSynth range_synth(
                geometry,
                makeNoiseConfig(kRangeSigmaKm,
                                4301U + static_cast<std::uint32_t>(100U * station_index)));
            fd::observations::synth::RangeRateSynth range_rate_synth(
                geometry,
                makeNoiseConfig(kRangeRateSigmaKmPerSec,
                                4302U + static_cast<std::uint32_t>(100U * station_index)));

            std::vector<fd::observations::synth::SyntheticObservationSample> valid_ranges;
            std::vector<fd::observations::synth::SyntheticObservationSample> valid_rates;
            valid_ranges.reserve(epochs.size());
            valid_rates.reserve(epochs.size());

            for (const double epoch : epochs) {
                if (range_synth.isVisibleAboveElevation(epoch,
                                                        station_info.station,
                                                        kElevationMaskDeg / dpr_c(),
                                                        truth_provider)) {
                    const auto range_sample = range_synth.generate(epoch,
                                                                   epoch,
                                                                   kCadenceSec,
                                                                   truth_provider);
                    const auto rate_sample = range_rate_synth.generate(epoch,
                                                                       epoch,
                                                                       kCadenceSec,
                                                                       kRangeRateCountTimeSec,
                                                                       truth_provider);
                    valid_ranges.push_back(singleSample(range_sample,
                                                        "Simulated range generation"));
                    valid_rates.push_back(singleSample(rate_sample,
                                                       "Simulated range-rate generation"));
                }
            }

            if (valid_ranges.empty()) {
                throw std::runtime_error("Elevation mask rejected every synthetic observation for " + station_name);
            }

            station_sets.push_back(StationObservationSet{
                std::move(station_info.station),
                station_info.stationNaif,
                std::move(geometry),
                std::move(valid_ranges),
                std::move(valid_rates)
            });
        }

        std::vector<VLBIObservationSet> vlbi_sets;
        vlbi_sets.reserve(kVLBIBaselines.size());
        for (const VLBIBaselineConfig& baseline : kVLBIBaselines) {
            vlbi_sets.push_back(generateVlbiObservationSet(baseline,
                                                           truth_provider,
                                                           start_et,
                                                           end_et,
                                                           start_et));
        }

        writeReport(kReportPath, station_sets, vlbi_sets);

        const std::vector<ReportRow> report_rows = makeSortedReportRows(station_sets);
        const std::size_t total_observation_count = report_rows.size();
        std::size_t total_vlbi_count = 0;
        for (const VLBIObservationSet& vlbi_set : vlbi_sets) {
            total_vlbi_count += vlbi_set.delays.size();
        }
        const ReportRow& first_row = report_rows.front();
        const StationObservationSet& first_station_set = *first_row.stationSet;
        const double first_epoch = first_station_set.ranges[first_row.sampleIndex].epochTdb;

        ComparisonStats range_comparison;
        ComparisonStats range_rate_comparison;
        ComparisonStats vlbi_comparison;
        for (const StationObservationSet& station_set : station_sets) {
            addCspiceComparisonForStation(station_set,
                                          range_comparison,
                                          range_rate_comparison);
        }
        for (const VLBIObservationSet& vlbi_set : vlbi_sets) {
            addCspiceComparisonForVlbi(vlbi_set, vlbi_comparison);
        }

        const fd::observations::synth::NoiseConfig no_noise = makeNoiseConfig(0.0, 0U, false);
        fd::observations::synth::RangeSynth cspice_range_synth(first_station_set.geometry, no_noise);
        fd::observations::synth::RangeRateSynth cspice_rate_synth(first_station_set.geometry, no_noise);
        const auto cspice_range = cspice_range_synth.generate(first_epoch, first_epoch, kCadenceSec);
        const auto cspice_rate =
            cspice_rate_synth.generate(first_epoch, first_epoch, kCadenceSec, kRangeRateCountTimeSec);
        const auto& first_cspice_range = singleSample(cspice_range, "First CSPICE range comparison");
        const auto& first_cspice_rate = singleSample(cspice_rate, "First CSPICE range-rate comparison");

        std::cout << "Generated " << total_observation_count
                  << " visible DSS-43/DSS-63 Voyager 1 station observations in " << kReportPath << '\n'
                  << "Generated " << total_vlbi_count
                  << " VLBI observations across " << vlbi_sets.size()
                  << " baselines in " << kReportPath << '\n';
        for (const VLBIObservationSet& vlbi_set : vlbi_sets) {
            std::cout << "  " << vlbi_set.stationOneName << '/' << vlbi_set.stationTwoName
                      << ": " << vlbi_set.delays.size() << " VLBI observations\n";
        }
        std::cout << "First UTC: " << utcFromEt(report_rows.front().stationSet->ranges[report_rows.front().sampleIndex].epochTdb) << '\n'
                  << "Last UTC : " << utcFromEt(report_rows.back().stationSet->ranges[report_rows.back().sampleIndex].epochTdb) << '\n'
                  << "Elevation mask: " << kElevationMaskDeg << " deg\n"
                  << "Range sigma: " << kRangeSigmaKm * 1000.0 << " m\n"
                  << "Range-rate count time: " << kRangeRateCountTimeSec << " s\n"
                  << "Range-rate sigma: " << rangeRateSigmaForCountTime() * 1000.0 << " m/s\n"
                  << "VLBI sigma: " << kVLBISigmaKm * 1000.0 << " m\n";

        std::cout << std::scientific << std::setprecision(12)
                  << "First observation CSPICE/RKF45 LT comparison\n"
                  << "  station                      : " << first_station_set.station.name() << '\n'
                  << "  utc                          : " << utcFromEt(first_epoch) << '\n'
                  << "  cspice range truth km        : " << first_cspice_range.truth << '\n'
                  << "  rkf45 range truth km         : "
                  << first_station_set.ranges[first_row.sampleIndex].truth << '\n'
                  << "  range delta m                : "
                  << (first_station_set.ranges[first_row.sampleIndex].truth
                      - first_cspice_range.truth)
                      * 1000.0 << '\n'
                  << "  cspice range-rate truth km/s : " << first_cspice_rate.truth << '\n'
                  << "  rkf45 range-rate truth km/s  : "
                  << first_station_set.rangeRates[first_row.sampleIndex].truth << '\n'
                  << "  range-rate delta mm/s        : "
                  << (first_station_set.rangeRates[first_row.sampleIndex].truth
                      - first_cspice_rate.truth)
                      * 1.0e6 << '\n'
                  << "Aggregate CSPICE/RKF45 source comparison\n"
                  << "  range samples                : " << range_comparison.count << '\n'
                  << "  range RMS delta m            : " << range_comparison.rms() * 1000.0 << '\n'
                  << "  range max abs delta m        : " << range_comparison.maxAbs * 1000.0 << '\n'
                  << "  range-rate samples           : " << range_rate_comparison.count << '\n'
                  << "  range-rate RMS delta mm/s    : "
                  << range_rate_comparison.rms() * 1.0e6 << '\n'
                  << "  range-rate max abs delta mm/s: "
                  << range_rate_comparison.maxAbs * 1.0e6 << '\n'
                  << "  VLBI samples                 : " << vlbi_comparison.count << '\n'
                  << "  VLBI RMS delta m             : " << vlbi_comparison.rms() * 1000.0 << '\n'
                  << "  VLBI max abs delta m         : " << vlbi_comparison.maxAbs * 1000.0 << '\n'
                  << "  truth ephemeris nodes        : " << truth_ephemeris.history().size() << '\n';

        kclear_c();
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Synthetic observation test failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
