#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "stations/StationCatalog.hpp"
#include "observations/synth/RangeRateSynth.hpp"
#include "observations/synth/RangeSynth.hpp"
#include "observations/synth/VLBISynth.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "perturbations/Shapiro.hpp"

#include <SpiceUsr.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using State6 = Eigen::Matrix<double, 6, 1>;

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-01-10T00:00:00";
constexpr std::array<const char*, 2> kStationNames {"DSS-43", "DSS-63"};
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
constexpr double kLightTimeToleranceSec = 1.0e-9;
constexpr int kMaxLightTimeIterations = 12;
constexpr double kRangeRateStepSec = 1.0;
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

struct ReportRow {
    const StationObservationSet* stationSet {nullptr};
    std::size_t sampleIndex {0};
};

struct ComputedRange {
    double rangeKm {0.0};
    double geometricRangeKm {0.0};
    double shapiroDelayKm {0.0};
    double lightTimeSec {0.0};
    double emitEpochTdb {0.0};
};

class SpiceErrorModeGuard {
public:
    SpiceErrorModeGuard() {
        erract_c("GET", static_cast<SpiceInt>(previous_action_.size()), previous_action_.data());
        errprt_c("GET", static_cast<SpiceInt>(previous_report_.size()), previous_report_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
        errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
    }

    SpiceErrorModeGuard(const SpiceErrorModeGuard&) = delete;
    SpiceErrorModeGuard& operator=(const SpiceErrorModeGuard&) = delete;

    ~SpiceErrorModeGuard() {
        erract_c("SET", 0, previous_action_.data());
        errprt_c("SET", 0, previous_report_.data());
    }

private:
    std::array<SpiceChar, 32> previous_action_ {};
    std::array<SpiceChar, 32> previous_report_ {};
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

[[nodiscard]] State6 interpolateState(const od::EphemerisInterpolator& ephemeris,
                                      double tdb,
                                      const std::string& context) {
    return asState6(ephemeris.interpolate(tdb), context);
}

[[nodiscard]] double shapiroDelayForGeometry(const Eigen::Vector3d& stationFromSun,
                                             const Eigen::Vector3d& targetFromSun) {
    const double geometric_range_km = (targetFromSun - stationFromSun).norm();
    return fd::perturbations::computeShapiroRangeDelay(
        stationFromSun.norm(),
        targetFromSun.norm(),
        geometric_range_km,
        fd::perturbations::kSunGravitationalParameterKm3PerSec2);
}

[[nodiscard]] ComputedRange computeRangeAtReceiveEpoch(const od::EphemerisInterpolator& spacecraftEphemeris,
                                                       const std::string& stationNaif,
                                                       double receiveEpoch) {
    const State6 station_state = spiceState(stationNaif, receiveEpoch, kCentralBody, "NONE");
    const Eigen::Vector3d station_position = station_state.segment<3>(0);

    double light_time_sec = 0.0;
    State6 spacecraft_state = interpolateState(spacecraftEphemeris,
                                               receiveEpoch,
                                               "Receive-epoch synthetic ephemeris interpolation");

    for (int iteration = 0; iteration < kMaxLightTimeIterations; ++iteration) {
        const double geometric_range_km = (spacecraft_state.segment<3>(0) - station_position).norm();
        const double updated_light_time_sec =
            geometric_range_km / fd::perturbations::kSpeedOfLightKmPerSec;

        if (std::abs(updated_light_time_sec - light_time_sec) < kLightTimeToleranceSec) {
            light_time_sec = updated_light_time_sec;
            break;
        }

        light_time_sec = updated_light_time_sec;
        spacecraft_state = interpolateState(spacecraftEphemeris,
                                            receiveEpoch - light_time_sec,
                                            "Emit-epoch synthetic ephemeris interpolation");
    }

    const Eigen::Vector3d target_position = spacecraft_state.segment<3>(0);
    const double geometric_range_km = (target_position - station_position).norm();
    const double shapiro_delay_km = shapiroDelayForGeometry(station_position, target_position);

    return ComputedRange{
        geometric_range_km + shapiro_delay_km,
        geometric_range_km,
        shapiro_delay_km,
        light_time_sec,
        receiveEpoch - light_time_sec
    };
}

[[nodiscard]] double computeVlbiDelayAtReceiveEpoch(const od::EphemerisInterpolator& spacecraftEphemeris,
                                                    const std::string& stationOneNaif,
                                                    const std::string& stationTwoNaif,
                                                    double receiveEpoch) {
    const ComputedRange station_one_range =
        computeRangeAtReceiveEpoch(spacecraftEphemeris, stationOneNaif, receiveEpoch);
    const ComputedRange station_two_range =
        computeRangeAtReceiveEpoch(spacecraftEphemeris, stationTwoNaif, receiveEpoch);
    return station_two_range.rangeKm - station_one_range.rangeKm;
}

[[nodiscard]] Eigen::Vector3d rotateJ2000ToItrf93(const Eigen::Vector3d& vector,
                                                  double epochTdb) {
    SpiceDouble rotation[3][3] = {};
    pxform_c(kFrame, "ITRF93", epochTdb, rotation);
    throwIfSpiceFailed("Failed to transform propagated line-of-sight into ITRF93");

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

[[nodiscard]] double computePropagatedElevationDeg(const od::EphemerisInterpolator& spacecraftEphemeris,
                                                   const od::Station& station,
                                                   const std::string& stationNaif,
                                                   double epochTdb) {
    const State6 spacecraft_state = interpolateState(spacecraftEphemeris,
                                                     epochTdb,
                                                     "Elevation-mask synthetic ephemeris interpolation");
    const State6 station_state = spiceState(stationNaif, epochTdb, kCentralBody, "NONE");
    const Eigen::Vector3d los_j2000 =
        spacecraft_state.segment<3>(0) - station_state.segment<3>(0);
    const Eigen::Vector3d los_itrf = rotateJ2000ToItrf93(los_j2000, epochTdb);
    const double range_km = los_itrf.norm();
    if (range_km <= 0.0) {
        throw std::runtime_error("Cannot compute elevation for zero station-spacecraft range.");
    }

    const double sin_elevation =
        std::clamp(los_itrf.normalized().dot(stationUpVectorItrf(station)), -1.0, 1.0);
    return std::asin(sin_elevation) * dpr_c();
}

[[nodiscard]] double drawGaussian(std::mt19937& rng, double sigma) {
    if (sigma == 0.0) {
        return 0.0;
    }
    std::normal_distribution<double> distribution(0.0, sigma);
    return distribution(rng);
}

[[nodiscard]] std::string utcFromEt(double et) {
    SpiceChar utc[128] = {0};
    et2utc_c(et, "ISOC", 3, static_cast<SpiceInt>(sizeof(utc)), utc);
    throwIfSpiceFailed("Failed to convert ET to UTC");
    return utc;
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
           << "# range_rate_sigma_km_s : " << kRangeRateSigmaKmPerSec << '\n'
           << "# vlbi_sigma_km         : " << kVLBISigmaKm << '\n'
           << "# cadence_seconds       : " << kCadenceSec << '\n'
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

} // namespace

int main() {
    try {
        const SpiceErrorModeGuard spice_error_mode;

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load synthetic observation meta-kernel");

        SpiceDouble start_et = 0.0;
        str2et_c(kStartUtc, &start_et);
        throwIfSpiceFailed("Failed to convert synthetic observation start epoch");

        const double end_et = start_et + kArcDurationSec;
        const double ephemeris_start_et = start_et - kSyntheticEphemerisLookbackSec;
        const double ephemeris_end_et = end_et + kRangeRateStepSec;

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

        const std::vector<double> epochs =
            fd::observations::synth::SyntheticObservation::makeEpochGrid(start_et, end_et, kCadenceSec);

        std::vector<StationObservationSet> station_sets;
        station_sets.reserve(kStationNames.size());

        for (std::size_t station_index = 0; station_index < kStationNames.size(); ++station_index) {
            const std::string station_name = kStationNames[station_index];
            const int station_naif_id = od::stationNaifIdFromName(station_name);
            const std::string station_naif = std::to_string(station_naif_id);
            od::Station station = od::buildStationFromKernel(station_name, station_naif_id, start_et);

            fd::observations::synth::GeometryConfig geometry;
            geometry.target = kTarget;
            geometry.stationName = station_name;
            geometry.frame = "J2000";
            geometry.aberrationCorrection = "LT";

            std::mt19937 range_rng(4301U + static_cast<std::uint32_t>(100U * station_index));
            std::mt19937 range_rate_rng(4302U + static_cast<std::uint32_t>(100U * station_index));

            std::vector<fd::observations::synth::SyntheticObservationSample> valid_ranges;
            std::vector<fd::observations::synth::SyntheticObservationSample> valid_rates;
            valid_ranges.reserve(epochs.size());
            valid_rates.reserve(epochs.size());

            for (const double epoch : epochs) {
                if (computePropagatedElevationDeg(truth_ephemeris,
                                                  station,
                                                  station_naif,
                                                  epoch) >= kElevationMaskDeg) {
                    const ComputedRange range =
                        computeRangeAtReceiveEpoch(truth_ephemeris, station_naif, epoch);
                    const ComputedRange next_range =
                        computeRangeAtReceiveEpoch(truth_ephemeris,
                                                   station_naif,
                                                   epoch + kRangeRateStepSec);
                    const double range_rate_truth =
                        (next_range.rangeKm - range.rangeKm) / kRangeRateStepSec;

                    const double range_noise = drawGaussian(range_rng, kRangeSigmaKm);
                    const double range_rate_noise =
                        drawGaussian(range_rate_rng, kRangeRateSigmaKmPerSec);

                    valid_ranges.push_back(fd::observations::synth::SyntheticObservationSample{
                        epoch,
                        range.rangeKm,
                        range_noise,
                        range.rangeKm + range_noise,
                        kRangeSigmaKm
                    });
                    valid_rates.push_back(fd::observations::synth::SyntheticObservationSample{
                        epoch,
                        range_rate_truth,
                        range_rate_noise,
                        range_rate_truth + range_rate_noise,
                        kRangeRateSigmaKmPerSec
                    });
                }
            }

            if (valid_ranges.empty()) {
                throw std::runtime_error("Elevation mask rejected every synthetic observation for " + station_name);
            }

            station_sets.push_back(StationObservationSet{
                std::move(station),
                station_naif,
                std::move(geometry),
                std::move(valid_ranges),
                std::move(valid_rates)
            });
        }

        if (station_sets.size() != 2U) {
            throw std::runtime_error("VLBI synthetic generation currently expects exactly two stations.");
        }

        std::mt19937 vlbi_rng(5301U);
        VLBIObservationSet vlbi_set;
        vlbi_set.stationOneName = station_sets[0].station.name();
        vlbi_set.stationTwoName = station_sets[1].station.name();
        vlbi_set.stationOneNaif = station_sets[0].stationNaif;
        vlbi_set.stationTwoNaif = station_sets[1].stationNaif;
        vlbi_set.delays.reserve(epochs.size());

        for (const double epoch : epochs) {
            const bool first_visible =
                computePropagatedElevationDeg(truth_ephemeris,
                                              station_sets[0].station,
                                              station_sets[0].stationNaif,
                                              epoch) >= kElevationMaskDeg;
            const bool second_visible =
                computePropagatedElevationDeg(truth_ephemeris,
                                              station_sets[1].station,
                                              station_sets[1].stationNaif,
                                              epoch) >= kElevationMaskDeg;
            if (!first_visible || !second_visible) {
                continue;
            }

            const double truth = computeVlbiDelayAtReceiveEpoch(truth_ephemeris,
                                                                vlbi_set.stationOneNaif,
                                                                vlbi_set.stationTwoNaif,
                                                                epoch);
            const double noise = drawGaussian(vlbi_rng, kVLBISigmaKm);
            vlbi_set.delays.push_back(fd::observations::synth::SyntheticObservationSample{
                epoch,
                truth,
                noise,
                truth + noise,
                kVLBISigmaKm
            });
        }

        if (vlbi_set.delays.empty()) {
            throw std::runtime_error("Elevation mask rejected every synthetic VLBI observation.");
        }

        std::vector<VLBIObservationSet> vlbi_sets;
        vlbi_sets.push_back(std::move(vlbi_set));

        writeReport(kReportPath, station_sets, vlbi_sets);

        const std::vector<ReportRow> report_rows = makeSortedReportRows(station_sets);
        const std::size_t total_observation_count = report_rows.size();
        const std::size_t total_vlbi_count = vlbi_sets.front().delays.size();
        const ReportRow& first_row = report_rows.front();
        const StationObservationSet& first_station_set = *first_row.stationSet;
        const double first_epoch = first_station_set.ranges[first_row.sampleIndex].epochTdb;
        const double first_vlbi_epoch = vlbi_sets.front().delays.front().epochTdb;

        fd::observations::synth::NoiseConfig no_noise;
        no_noise.enabled = false;
        no_noise.sigma = 0.0;
        fd::observations::synth::RangeSynth cspice_range_synth(first_station_set.geometry, no_noise);
        fd::observations::synth::RangeRateSynth cspice_rate_synth(first_station_set.geometry, no_noise);
        fd::observations::synth::VLBIConfig vlbi_config;
        vlbi_config.target = kTarget;
        vlbi_config.stationOneName = vlbi_sets.front().stationOneName;
        vlbi_config.stationTwoName = vlbi_sets.front().stationTwoName;
        vlbi_config.frame = kFrame;
        vlbi_config.aberrationCorrection = "LT";
        vlbi_config.minimumElevationRad = kElevationMaskDeg / dpr_c();
        fd::observations::synth::VLBISynth cspice_vlbi_synth(vlbi_config, no_noise);
        const auto cspice_range = cspice_range_synth.generate(first_epoch, first_epoch, kCadenceSec);
        const auto cspice_rate = cspice_rate_synth.generate(first_epoch, first_epoch, kCadenceSec);
        const auto cspice_vlbi = cspice_vlbi_synth.generate(first_vlbi_epoch,
                                                            first_vlbi_epoch,
                                                            kCadenceSec);
        if (cspice_range.size() != 1U || cspice_rate.size() != 1U || cspice_vlbi.size() != 1U) {
            throw std::runtime_error("Failed to generate one-sample CSPICE comparison observation.");
        }

        std::cout << "Generated " << total_observation_count
                  << " visible DSS-43/DSS-63 Voyager 1 station observations in " << kReportPath << '\n'
                  << "Generated " << total_vlbi_count
                  << " DSS-43/DSS-63 VLBI observations in " << kReportPath << '\n'
                  << "First UTC: " << utcFromEt(report_rows.front().stationSet->ranges[report_rows.front().sampleIndex].epochTdb) << '\n'
                  << "Last UTC : " << utcFromEt(report_rows.back().stationSet->ranges[report_rows.back().sampleIndex].epochTdb) << '\n'
                  << "Elevation mask: " << kElevationMaskDeg << " deg\n"
                  << "Range sigma: " << kRangeSigmaKm * 1000.0 << " m\n"
                  << "Range-rate sigma: " << kRangeRateSigmaKmPerSec * 1000.0 << " m/s\n"
                  << "VLBI sigma: " << kVLBISigmaKm * 1000.0 << " m\n";

        std::cout << std::scientific << std::setprecision(12)
                  << "First observation CSPICE/RKF45 LT comparison\n"
                  << "  station                      : " << first_station_set.station.name() << '\n'
                  << "  utc                          : " << utcFromEt(first_epoch) << '\n'
                  << "  cspice range truth km        : " << cspice_range.front().truth << '\n'
                  << "  rkf45 range truth km         : "
                  << first_station_set.ranges[first_row.sampleIndex].truth << '\n'
                  << "  range delta m                : "
                  << (first_station_set.ranges[first_row.sampleIndex].truth - cspice_range.front().truth)
                      * 1000.0 << '\n'
                  << "  cspice range-rate truth km/s : " << cspice_rate.front().truth << '\n'
                  << "  rkf45 range-rate truth km/s  : "
                  << first_station_set.rangeRates[first_row.sampleIndex].truth << '\n'
                  << "  range-rate delta mm/s        : "
                  << (first_station_set.rangeRates[first_row.sampleIndex].truth - cspice_rate.front().truth)
                      * 1.0e6 << '\n'
                  << "  cspice VLBI truth km         : " << cspice_vlbi.front().truth << '\n'
                  << "  rkf45 VLBI truth km          : " << vlbi_sets.front().delays.front().truth << '\n'
                  << "  VLBI delta m                 : "
                  << (vlbi_sets.front().delays.front().truth - cspice_vlbi.front().truth) * 1000.0 << '\n'
                  << "  truth ephemeris nodes        : " << truth_ephemeris.history().size() << '\n';

        kclear_c();
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Synthetic observation test failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
