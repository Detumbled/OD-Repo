#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "observations/synth/RangeRateSynth.hpp"
#include "observations/synth/RangeSynth.hpp"
#include "observations/synth/VLBISynth.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "perturbations/Shapiro.hpp"
#include "stations/StationCatalog.hpp"
#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using State6 = Eigen::Matrix<double, 6, 1>;
using od::throwIfSpiceFailed;

struct VLBIBaselineConfig {
    const char* stationOneName;
    const char* stationTwoName;
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

struct SourceComparisonStats {
    ComparisonStats range;
    ComparisonStats rangeRate;
    ComparisonStats vlbi;
};

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-01-10T00:00:00";
constexpr std::array<const char*, 2> kStationNames {"DSS-43", "DSS-63"};
constexpr std::array<VLBIBaselineConfig, 2> kVLBIBaselines {{
    {"DSS-43", "DSS-63"},
    {"DSS-14", "DSS-43"}
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
constexpr double kElevationMaskDeg = 10.0;
constexpr double kRangeRateCountTimeSec = 60.0;
constexpr double kSyntheticEphemerisLookbackSec = 3.0 * 3600.0;
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
                           bool required) {
    const bool has_coverage = spiceStateAvailable(body.name, startEpoch, kCentralBody)
        && spiceStateAvailable(body.name, endEpoch, kCentralBody);

    if (has_coverage) {
        thirdBody.addBody(body);
        return;
    }

    if (required) {
        throw std::runtime_error("Required third-body SPK coverage is missing for " + body.name);
    }
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

    return ephemeris;
}

[[nodiscard]] fd::observations::synth::NoiseConfig noNoise() {
    fd::observations::synth::NoiseConfig noise;
    noise.enabled = false;
    noise.sigma = 0.0;
    return noise;
}

[[nodiscard]] fd::observations::synth::GeometryConfig makeGeometry(const std::string& stationName) {
    fd::observations::synth::GeometryConfig geometry;
    geometry.target = kTarget;
    geometry.stationName = stationName;
    geometry.frame = kFrame;
    geometry.aberrationCorrection = "NONE";
    return geometry;
}

[[nodiscard]] const fd::observations::synth::SyntheticObservationSample& singleSample(
    const std::vector<fd::observations::synth::SyntheticObservationSample>& samples,
    const std::string& context) {
    if (samples.size() != 1U) {
        throw std::runtime_error(context + " did not generate exactly one sample.");
    }
    return samples.front();
}

void addStationComparisons(const fd::observations::synth::TargetStateProvider& targetProvider,
                           double startEpoch,
                           double endEpoch,
                           SourceComparisonStats& stats) {
    const std::vector<double> epochs =
        fd::observations::synth::SyntheticObservation::makeEpochGrid(startEpoch,
                                                                     endEpoch,
                                                                     kCadenceSec);

    for (const char* station_name : kStationNames) {
        const int station_naif_id = od::stationNaifIdFromName(station_name);
        const od::Station station = od::buildStationFromKernel(station_name,
                                                               station_naif_id,
                                                               startEpoch);

        fd::observations::synth::RangeSynth simulated_range_synth(makeGeometry(station_name),
                                                                  noNoise());
        fd::observations::synth::RangeRateSynth simulated_rate_synth(makeGeometry(station_name),
                                                                     noNoise());
        fd::observations::synth::RangeSynth cspice_range(makeGeometry(station_name), noNoise());
        fd::observations::synth::RangeRateSynth cspice_rate(makeGeometry(station_name), noNoise());

        for (const double epoch : epochs) {
            if (!simulated_range_synth.isVisibleAboveElevation(epoch,
                                                               station,
                                                               kElevationMaskDeg / dpr_c(),
                                                               targetProvider)) {
                continue;
            }

            const double simulated_range =
                singleSample(simulated_range_synth.generate(epoch, epoch, kCadenceSec, targetProvider),
                             "Simulated range").truth;
            const double simulated_rate =
                singleSample(simulated_rate_synth.generate(epoch,
                                                           epoch,
                                                           kCadenceSec,
                                                           kRangeRateCountTimeSec,
                                                           targetProvider),
                             "Simulated range-rate").truth;

            stats.range.add(simulated_range
                            - singleSample(cspice_range.generate(epoch, epoch, kCadenceSec),
                                           "CSPICE relativistic range").truth);
            stats.rangeRate.add(simulated_rate
                                - singleSample(cspice_rate.generate(epoch,
                                                                    epoch,
                                                                    kCadenceSec,
                                                                    kRangeRateCountTimeSec),
                                               "CSPICE relativistic range-rate").truth);
        }
    }
}

[[nodiscard]] fd::observations::synth::VLBIConfig makeVlbiConfig(
    const VLBIBaselineConfig& baseline,
    double minimumElevationRad) {
    fd::observations::synth::VLBIConfig vlbi;
    vlbi.target = kTarget;
    vlbi.stationOneName = baseline.stationOneName;
    vlbi.stationTwoName = baseline.stationTwoName;
    vlbi.frame = kFrame;
    vlbi.aberrationCorrection = "NONE";
    vlbi.minimumElevationRad = minimumElevationRad;
    return vlbi;
}

void addVlbiComparisons(const fd::observations::synth::TargetStateProvider& targetProvider,
                        double startEpoch,
                        double endEpoch,
    SourceComparisonStats& stats) {
    for (const VLBIBaselineConfig& baseline : kVLBIBaselines) {
        fd::observations::synth::VLBISynth simulated_vlbi(
            makeVlbiConfig(baseline, kElevationMaskDeg / dpr_c()),
            noNoise());
        const auto simulated_samples = simulated_vlbi.generate(startEpoch,
                                                               endEpoch,
                                                               kVLBICadenceSec,
                                                               targetProvider);

        fd::observations::synth::VLBISynth cspice_vlbi(
            makeVlbiConfig(baseline, -1.57079632679489661923),
            noNoise());

        for (const auto& simulated : simulated_samples) {
            stats.vlbi.add(simulated.truth
                           - singleSample(cspice_vlbi.generate(simulated.epochTdb,
                                                               simulated.epochTdb,
                                                               kVLBICadenceSec),
                                          "CSPICE relativistic VLBI").truth);
        }
    }
}

[[nodiscard]] od::EphemerisInterpolator buildTruthEphemeris(double startEpoch, double endEpoch) {
    od::RKF45Integrator::Options options;
    options.absoluteTolerance = 1.0e-5;
    options.relativeTolerance = 1.0e-12;
    options.initialStep = 30.0;
    options.minimumStep = 1.0e-4;
    options.maximumStep = 300.0;
    options.maximumSteps = 200000;
    const od::RKF45Integrator integrator(options);

    fd::perturbations::ThirdBodyGravity third_body(kCentralBody, kFrame);
    addThirdBodyIfCovered(third_body, {kJupiter, kJupiterBodyMuKm3PerSec2}, startEpoch, endEpoch, true);
    addThirdBodyIfCovered(third_body, {kIo, kIoMuKm3PerSec2}, startEpoch, endEpoch, false);
    addThirdBodyIfCovered(third_body, {kEuropa, kEuropaMuKm3PerSec2}, startEpoch, endEpoch, false);
    addThirdBodyIfCovered(third_body, {kGanymede, kGanymedeMuKm3PerSec2}, startEpoch, endEpoch, false);
    addThirdBodyIfCovered(third_body, {kCallisto, kCallistoMuKm3PerSec2}, startEpoch, endEpoch, false);
    addThirdBodyIfCovered(third_body,
                          {kSaturnBarycenter, kSaturnSystemMuKm3PerSec2},
                          startEpoch,
                          endEpoch,
                          true);

    fd::perturbations::SolarRadiationPressure srp(kCentralBody,
                                                  kFrame,
                                                  kVoyagerReflectivity,
                                                  kVoyagerAreaM2,
                                                  kVoyagerMassKg);
    const auto dynamics = makeDynamics(third_body, srp);
    const State6 truth_initial_state = spiceState(kTarget, startEpoch, kCentralBody, "NONE");
    return propagateEphemerisSpan(integrator,
                                  dynamics,
                                  startEpoch,
                                  truth_initial_state,
                                  startEpoch - kSyntheticEphemerisLookbackSec,
                                  endEpoch + 0.5 * kRangeRateCountTimeSec);
}

void printStats(const SourceComparisonStats& stats) {
    std::cout << std::scientific << std::setprecision(12)
              << "Synthetic source comparison: simulated provider minus CSPICE custom relativistic provider\n"
              << "  range samples              : " << stats.range.count << '\n'
              << "  range RMS delta m          : " << stats.range.rms() * 1000.0 << '\n'
              << "  range max abs delta m      : " << stats.range.maxAbs * 1000.0 << '\n'
              << "  range-rate samples         : " << stats.rangeRate.count << '\n'
              << "  range-rate RMS mm/s        : " << stats.rangeRate.rms() * 1.0e6 << '\n'
              << "  range-rate max mm/s        : " << stats.rangeRate.maxAbs * 1.0e6 << '\n'
              << "  VLBI samples               : " << stats.vlbi.count << '\n'
              << "  VLBI RMS delta m           : " << stats.vlbi.rms() * 1000.0 << '\n'
              << "  VLBI max abs delta m       : " << stats.vlbi.maxAbs * 1000.0 << '\n';
}

[[nodiscard]] bool validateStats(const SourceComparisonStats& stats) {
    bool passed = true;
    if (stats.range.count == 0U || stats.rangeRate.count == 0U || stats.vlbi.count == 0U) {
        std::cerr << "FAIL: source comparison generated an empty comparison set.\n";
        passed = false;
    }
    return passed;
}

} // namespace

int main() {
    try {
        const od::SpiceErrorModeGuard spice_error_mode("RETURN", "NONE");

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load synthetic source comparison meta-kernel");

        SpiceDouble start_et = 0.0;
        str2et_c(kStartUtc, &start_et);
        throwIfSpiceFailed("Failed to convert synthetic source comparison start epoch");
        const double end_et = start_et + kArcDurationSec;

        const od::EphemerisInterpolator truth_ephemeris = buildTruthEphemeris(start_et, end_et);
        const fd::observations::synth::InterpolatedTargetStateProvider truth_provider(
            truth_ephemeris,
            "Synthetic source comparison ephemeris interpolation");

        SourceComparisonStats stats;
        addStationComparisons(truth_provider, start_et, end_et, stats);
        addVlbiComparisons(truth_provider, start_et, end_et, stats);
        printStats(stats);

        kclear_c();
        return validateStats(stats) ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "Synthetic source comparison failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
