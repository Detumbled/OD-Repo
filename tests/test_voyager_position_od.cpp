#include "stations/StationCatalog.hpp"
#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "filters/WLS.hpp"
#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"
#include "perturbations/Shapiro.hpp"

#include <SpiceUsr.h>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using State6 = Eigen::Matrix<double, 6, 1>;

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kTarget = "-31";
constexpr const char* kDefaultStationName = "DSS-43";
constexpr const char* kDefaultStationNaif = "399043";
constexpr const char* kCentralBody = "SUN";
constexpr const char* kFrame = "J2000";
constexpr const char* kJupiter = "599";
constexpr const char* kIo = "501";
constexpr const char* kEuropa = "502";
constexpr const char* kGanymede = "503";
constexpr const char* kCallisto = "504";
constexpr const char* kSaturnBarycenter = "6";

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
constexpr double kEphemerisLightTimePaddingSec = 600.0;
constexpr int kMaxSolverIterations = 5;
constexpr double kAdoptedRangeSigmaFloorKm = 0.05; //50 meters
constexpr double kAdoptedRangeRateSigmaFloorKmPerSec = 1.0e-4; //

constexpr const char* kReportPath = "../tests/voyager_position_estimation_report.txt";
constexpr const char* kPostfitDiagnosticsCsvPath = "../tests/voyager_od_postfit_diagnostics.csv";
constexpr const char* kTrajectoryErrorCsvPath = "../tests/voyager_od_trajectory_error.csv";
constexpr const char* kObservabilityWindowsCsvPath = "../tests/voyager_station_observability_windows.csv";

struct Observation {
    std::string utc;
    std::string stationName {kDefaultStationName};
    std::string stationNaif {kDefaultStationNaif};
    double epochTdb {0.0};
    double rangeTruthKm {0.0};
    double rangeNoiseKm {0.0};
    double rangeObservedKm {0.0};
    double rangeSigmaKm {0.0};
    double rangeRateTruthKmPerSec {0.0};
    double rangeRateNoiseKmPerSec {0.0};
    double rangeRateObservedKmPerSec {0.0};
    double rangeRateSigmaKmPerSec {0.0};
};

struct ComputedRange {
    double rangeKm {0.0};
    double geometricRangeKm {0.0};
    double shapiroDelayKm {0.0};
    double lightTimeSec {0.0};
    double emitEpochTdb {0.0};
    State6 emittedSpacecraftState {State6::Zero()};
    State6 stationReceiveState {State6::Zero()};
};

struct Prediction {
    Eigen::VectorXd values;
    std::vector<ComputedRange> ranges;
};

struct ResidualStats {
    double rangeRmsKm {0.0};
    double rangeRateRmsKmPerSec {0.0};
    double weightedRms {0.0};
};

struct IterationRecord {
    int iteration {0};
    double weightedRms {0.0};
    double rangeRmsKm {0.0};
    double rangeRateRmsKmPerSec {0.0};
    double correctionPositionNormKm {0.0};
    double correctionVelocityNormKmPerSec {0.0};
    double positionErrorKm {0.0};
    double velocityErrorKmPerSec {0.0};
};

class SpiceErrorModeGuard {
public:
    SpiceErrorModeGuard() {
        erract_c("GET", static_cast<SpiceInt>(previousAction_.size()), previousAction_.data());
        errprt_c("GET", static_cast<SpiceInt>(previousReport_.size()), previousReport_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
        errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
    }

    SpiceErrorModeGuard(const SpiceErrorModeGuard&) = delete;
    SpiceErrorModeGuard& operator=(const SpiceErrorModeGuard&) = delete;

    ~SpiceErrorModeGuard() {
        erract_c("SET", 0, previousAction_.data());
        errprt_c("SET", 0, previousReport_.data());
    }

private:
    std::array<SpiceChar, 32> previousAction_ {};
    std::array<SpiceChar, 32> previousReport_ {};
};

void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar shortMessage[1841] = {0};
    SpiceChar longMessage[1841] = {0};
    getmsg_c("SHORT", sizeof(shortMessage), shortMessage);
    getmsg_c("LONG", sizeof(longMessage), longMessage);
    reset_c();

    throw std::runtime_error(context + ": " + shortMessage + " | " + longMessage);
}

[[nodiscard]] std::filesystem::path firstExistingPath(const std::vector<std::filesystem::path>& candidates,
                                                      const std::string& label) {
    for (const auto& candidate : candidates) {
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
    }

    std::ostringstream message;
    message << "Could not find " << label << ". Tried:";
    for (const auto& candidate : candidates) {
        message << " " << candidate.string();
    }
    throw std::runtime_error(message.str());
}

[[nodiscard]] bool parseDoubleToken(const std::string& token, double& value) {
    char* end = nullptr;
    value = std::strtod(token.c_str(), &end);
    return end != token.c_str() && *end == '\0';
}

[[nodiscard]] std::string stationNaifFromName(const std::string& stationName) {
    try {
        return std::to_string(od::stationNaifIdFromName(stationName));
    } catch (const std::invalid_argument&) {
        double numericStationId = 0.0;
        if (parseDoubleToken(stationName, numericStationId)) {
            return stationName;
        }
        throw;
    }
}

[[nodiscard]] std::string utcFromEt(double et) {
    SpiceChar utc[128] = {0};
    et2utc_c(et, "ISOC", 3, static_cast<SpiceInt>(sizeof(utc)), utc);
    throwIfSpiceFailed("Failed to convert ET to UTC");
    return utc;
}

[[nodiscard]] std::vector<Observation> readObservationReport(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open observation report: " + path.string());
    }

    std::vector<Observation> observations;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') {
            continue;
        }

        std::istringstream parser(line);
        Observation observation;
        std::string epoch_or_station;
        parser >> observation.utc
               >> epoch_or_station;

        if (!parser) {
            throw std::runtime_error("Failed to parse observation row prefix: " + line);
        }

        double parsed_epoch = 0.0;
        if (parseDoubleToken(epoch_or_station, parsed_epoch)) {
            observation.stationName = kDefaultStationName;
            observation.stationNaif = kDefaultStationNaif;
            observation.epochTdb = parsed_epoch;
        } else {
            observation.stationName = epoch_or_station;
            observation.stationNaif = stationNaifFromName(observation.stationName);
            parser >> observation.epochTdb;
        }

        parser >> observation.rangeTruthKm
               >> observation.rangeNoiseKm
               >> observation.rangeObservedKm
               >> observation.rangeSigmaKm
               >> observation.rangeRateTruthKmPerSec
               >> observation.rangeRateNoiseKmPerSec
               >> observation.rangeRateObservedKmPerSec
               >> observation.rangeRateSigmaKmPerSec;

        if (!parser) {
            throw std::runtime_error("Failed to parse observation row: " + line);
        }
        if (observation.rangeSigmaKm <= 0.0 || observation.rangeRateSigmaKmPerSec <= 0.0) {
            throw std::runtime_error("Observation row has non-positive sigma: " + line);
        }

        observations.push_back(std::move(observation));
    }

    if (observations.size() < 4) {
        throw std::runtime_error("Observation report does not contain enough samples for OD.");
    }

    return observations;
}

[[nodiscard]] State6 spiceState(const std::string& target,
                                double tdb,
                                const std::string& observer,
                                const std::string& aberrationCorrection = "LT") {
    SpiceDouble spiceStateVector[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble lightTime = 0.0;
    spkezr_c(target.c_str(),
             tdb,
             kFrame,
             aberrationCorrection.c_str(),
             observer.c_str(),
             spiceStateVector,
             &lightTime);
    throwIfSpiceFailed("Failed to fetch SPICE state for " + target + " relative to " + observer);

    State6 state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spiceStateVector[i];
    }
    return state;
}

[[nodiscard]] bool spiceStateAvailable(const std::string& target,
                                       double tdb,
                                       const std::string& observer) {
    SpiceDouble spiceStateVector[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble lightTime = 0.0;
    spkezr_c(target.c_str(), tdb, kFrame, "NONE", observer.c_str(), spiceStateVector, &lightTime);
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
    const bool hasCoverage = spiceStateAvailable(body.name, startEpoch, kCentralBody)
        && spiceStateAvailable(body.name, endEpoch, kCentralBody);

    if (hasCoverage) {
        thirdBody.addBody(body);
        activeBodies.push_back(body.name);
        return;
    }

    if (required) {
        throw std::runtime_error("Required third-body SPK coverage is missing for " + body.name);
    }

    skippedBodies.push_back(body.name);
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

    State6 fixedState;
    fixedState = state;
    return fixedState;
}

[[nodiscard]] od::EphemerisInterpolator propagateEphemeris(
    const od::RKF45Integrator& integrator,
    const od::RKF45Integrator::DerivativeFunction& dynamics,
    double referenceEpoch,
    const State6& referenceState,
    double finalEpoch) {
    const od::RKF45Integrator::State dynamicInitial = referenceState;
    const auto result = integrator.integrate(referenceEpoch, dynamicInitial, finalEpoch, dynamics);
    (void) asState6(result.state, "Propagation");

    od::EphemerisInterpolator ephemeris;
    appendHistory(ephemeris, result);
    return ephemeris;
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
    bool addedSegment = false;

    if (startEpoch < referenceEpoch) {
        const od::RKF45Integrator::State dynamicInitial = referenceState;
        const auto backward = integrator.integrate(referenceEpoch, dynamicInitial, startEpoch, dynamics);
        (void) asState6(backward.state, "Backward propagation");
        appendHistory(ephemeris, backward);
        addedSegment = true;
    }

    if (finalEpoch > referenceEpoch) {
        const od::RKF45Integrator::State dynamicInitial = referenceState;
        const auto forward = integrator.integrate(referenceEpoch, dynamicInitial, finalEpoch, dynamics);
        (void) asState6(forward.state, "Forward propagation");
        appendHistory(ephemeris, forward);
        addedSegment = true;
    }

    if (!addedSegment) {
        return propagateEphemeris(integrator, dynamics, referenceEpoch, referenceState, referenceEpoch);
    }

    return ephemeris;
}

[[nodiscard]] State6 interpolateState(const od::EphemerisInterpolator& ephemeris,
                                      double tdb,
                                      const std::string& context) {
    return asState6(ephemeris.interpolate(tdb), context);
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

[[nodiscard]] double shapiroDelayForGeometry(const Eigen::Vector3d& stationFromSun,
                                             const Eigen::Vector3d& targetFromSun) {
    const double geometricRangeKm = (targetFromSun - stationFromSun).norm();
    return fd::perturbations::computeShapiroRangeDelay(
        stationFromSun.norm(),
        targetFromSun.norm(),
        geometricRangeKm,
        fd::perturbations::kSunGravitationalParameterKm3PerSec2);
}

[[nodiscard]] double maximumObservedLightTimeSec(const std::vector<Observation>& observations) {
    double maximumRangeKm = 0.0;
    for (const Observation& observation : observations) {
        maximumRangeKm = std::max(maximumRangeKm, std::abs(observation.rangeObservedKm));
    }
    if (!(maximumRangeKm > 0.0) || !std::isfinite(maximumRangeKm)) {
        throw std::runtime_error("Observations do not contain a usable range for ephemeris lookback.");
    }
    return maximumRangeKm / fd::perturbations::kSpeedOfLightKmPerSec;
}

[[nodiscard]] std::pair<double, double> ephemerisBoundsForObservations(
    const std::vector<Observation>& observations) {
    if (observations.empty()) {
        throw std::invalid_argument("Cannot build an ephemeris span for an empty observation set.");
    }

    double firstEpoch = observations.front().epochTdb;
    double lastEpoch = observations.front().epochTdb;
    for (const Observation& observation : observations) {
        firstEpoch = std::min(firstEpoch, observation.epochTdb);
        lastEpoch = std::max(lastEpoch, observation.epochTdb);
    }

    const double lookbackSec = maximumObservedLightTimeSec(observations)
        + kRangeRateStepSec
        + kEphemerisLightTimePaddingSec;
    return {firstEpoch - lookbackSec, lastEpoch + kRangeRateStepSec};
}

[[nodiscard]] od::EphemerisInterpolator propagateObservationEphemeris(
    const std::vector<Observation>& observations,
    const od::RKF45Integrator& integrator,
    const od::RKF45Integrator::DerivativeFunction& dynamics,
    double referenceEpoch,
    const State6& referenceState) {
    const auto [startEpoch, finalEpoch] = ephemerisBoundsForObservations(observations);
    return propagateEphemerisSpan(integrator,
                                  dynamics,
                                  referenceEpoch,
                                  referenceState,
                                  startEpoch,
                                  finalEpoch);
}

[[nodiscard]] ComputedRange computeRangeAtReceiveEpoch(const od::EphemerisInterpolator& spacecraftEphemeris,
                                                       const std::string& stationNaif,
                                                       double receiveEpoch) {
    const State6 stationState = spiceState(stationNaif, receiveEpoch, kCentralBody, "NONE");
    const Eigen::Vector3d stationPosition = stationState.segment<3>(0);

    double lightTimeSec = 0.0;
    State6 spacecraftState = interpolateState(spacecraftEphemeris,
                                              receiveEpoch,
                                              "Receive-epoch ephemeris interpolation");

    for (int iteration = 0; iteration < kMaxLightTimeIterations; ++iteration) {
        const double geometricRangeKm = (spacecraftState.segment<3>(0) - stationPosition).norm();
        const double updatedLightTimeSec =
            geometricRangeKm / fd::perturbations::kSpeedOfLightKmPerSec;

        if (std::abs(updatedLightTimeSec - lightTimeSec) < kLightTimeToleranceSec) {
            lightTimeSec = updatedLightTimeSec;
            break;
        }

        lightTimeSec = updatedLightTimeSec;
        spacecraftState = interpolateState(spacecraftEphemeris,
                                           receiveEpoch - lightTimeSec,
                                           "Transmit-epoch ephemeris interpolation");
    }

    const Eigen::Vector3d targetPosition = spacecraftState.segment<3>(0);
    const double geometricRangeKm = (targetPosition - stationPosition).norm();
    const double shapiroDelayKm = shapiroDelayForGeometry(stationPosition, targetPosition);

    ComputedRange computed;
    computed.rangeKm = geometricRangeKm + shapiroDelayKm;
    computed.geometricRangeKm = geometricRangeKm;
    computed.shapiroDelayKm = shapiroDelayKm;
    computed.lightTimeSec = lightTimeSec;
    computed.emitEpochTdb = receiveEpoch - lightTimeSec;
    computed.emittedSpacecraftState = spacecraftState;
    computed.stationReceiveState = stationState;
    return computed;
}

[[nodiscard]] Prediction predictObservations(const std::vector<Observation>& observations,
                                             const od::RKF45Integrator& integrator,
                                             const od::RKF45Integrator::DerivativeFunction& dynamics,
                                             double referenceEpoch,
                                             const State6& referenceState) {
    const od::EphemerisInterpolator spacecraftEphemeris =
        propagateObservationEphemeris(observations,
                                      integrator,
                                      dynamics,
                                      referenceEpoch,
                                      referenceState);

    Prediction prediction;
    prediction.values.resize(static_cast<Eigen::Index>(2 * observations.size()));
    prediction.ranges.reserve(observations.size());

    for (std::size_t i = 0; i < observations.size(); ++i) {
        const ComputedRange range = computeRangeAtReceiveEpoch(spacecraftEphemeris,
                                                               observations[i].stationNaif,
                                                               observations[i].epochTdb);
        const ComputedRange nextRange = computeRangeAtReceiveEpoch(spacecraftEphemeris,
                                                                   observations[i].stationNaif,
                                                                   observations[i].epochTdb
                                                                       + kRangeRateStepSec);

        prediction.values[static_cast<Eigen::Index>(2 * i)] = range.rangeKm;
        prediction.values[static_cast<Eigen::Index>(2 * i + 1)] =
            (nextRange.rangeKm - range.rangeKm) / kRangeRateStepSec;
        prediction.ranges.push_back(range);
    }

    return prediction;
}

[[nodiscard]] Eigen::VectorXd measuredVector(const std::vector<Observation>& observations) {
    Eigen::VectorXd measured(static_cast<Eigen::Index>(2 * observations.size()));
    for (std::size_t i = 0; i < observations.size(); ++i) {
        measured[static_cast<Eigen::Index>(2 * i)] = observations[i].rangeObservedKm;
        measured[static_cast<Eigen::Index>(2 * i + 1)] = observations[i].rangeRateObservedKmPerSec;
    }
    return measured;
}

[[nodiscard]] Eigen::MatrixXd measurementCovariance(const std::vector<Observation>& observations) {
    Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(2 * observations.size()),
                                                       static_cast<Eigen::Index>(2 * observations.size()));
    for (std::size_t i = 0; i < observations.size(); ++i) {
        const double adoptedRangeSigmaKm = std::max(observations[i].rangeSigmaKm,
                                                    kAdoptedRangeSigmaFloorKm);
        const double adoptedRangeRateSigmaKmPerSec = std::max(observations[i].rangeRateSigmaKmPerSec,
                                                              kAdoptedRangeRateSigmaFloorKmPerSec);

        covariance(static_cast<Eigen::Index>(2 * i), static_cast<Eigen::Index>(2 * i)) =
            adoptedRangeSigmaKm * adoptedRangeSigmaKm;
        covariance(static_cast<Eigen::Index>(2 * i + 1), static_cast<Eigen::Index>(2 * i + 1)) =
            adoptedRangeRateSigmaKmPerSec * adoptedRangeRateSigmaKmPerSec;
    }
    return covariance;
}

[[nodiscard]] Eigen::MatrixXd finiteDifferenceDesignMatrix(
    const std::vector<Observation>& observations,
    const od::RKF45Integrator& integrator,
    const od::RKF45Integrator::DerivativeFunction& dynamics,
    double referenceEpoch,
    const State6& nominalState) {
    Eigen::Matrix<double, 6, 1> perturbationSteps;
    perturbationSteps << 1.0, 1.0, 1.0, 1.0e-4, 1.0e-4, 1.0e-4;

    Eigen::MatrixXd design(static_cast<Eigen::Index>(2 * observations.size()), 6);

    for (Eigen::Index column = 0; column < 6; ++column) {
        State6 plusState = nominalState;
        State6 minusState = nominalState;
        plusState[column] += perturbationSteps[column];
        minusState[column] -= perturbationSteps[column];

        const Eigen::VectorXd plus =
            predictObservations(observations, integrator, dynamics, referenceEpoch, plusState).values;
        const Eigen::VectorXd minus =
            predictObservations(observations, integrator, dynamics, referenceEpoch, minusState).values;

        design.col(column) = (plus - minus) / (2.0 * perturbationSteps[column]);
    }

    if (!design.allFinite()) {
        throw std::runtime_error("Finite-difference design matrix contains non-finite values.");
    }

    return design;
}

[[nodiscard]] ResidualStats computeStats(const Eigen::VectorXd& residuals,
                                         const std::vector<Observation>& observations) {
    if (residuals.size() != static_cast<Eigen::Index>(2 * observations.size())) {
        throw std::runtime_error("Residual vector has inconsistent size.");
    }

    double rangeSum = 0.0;
    double rangeRateSum = 0.0;
    double weightedSum = 0.0;

    for (std::size_t i = 0; i < observations.size(); ++i) {
        const double rangeResidual = residuals[static_cast<Eigen::Index>(2 * i)];
        const double rangeRateResidual = residuals[static_cast<Eigen::Index>(2 * i + 1)];

        rangeSum += rangeResidual * rangeResidual;
        rangeRateSum += rangeRateResidual * rangeRateResidual;
        weightedSum += std::pow(rangeResidual / observations[i].rangeSigmaKm, 2.0)
            + std::pow(rangeRateResidual / observations[i].rangeRateSigmaKmPerSec, 2.0);
    }

    const double count = static_cast<double>(observations.size());
    return ResidualStats{
        std::sqrt(rangeSum / count),
        std::sqrt(rangeRateSum / count),
        std::sqrt(weightedSum / static_cast<double>(residuals.size()))
    };
}

[[nodiscard]] State6 stateAtEpoch(const fd::filters::WLS& filter) {
    if (filter.state().size() != 6) {
        throw std::runtime_error("Filter state dimension is not 6.");
    }
    return filter.state();
}

[[nodiscard]] std::vector<std::string> uniqueStationLabels(const std::vector<Observation>& observations) {
    std::vector<std::string> labels;
    for (const Observation& observation : observations) {
        const std::string label = observation.stationName + "(" + observation.stationNaif + ")";
        if (std::find(labels.begin(), labels.end(), label) == labels.end()) {
            labels.push_back(label);
        }
    }
    return labels;
}

[[nodiscard]] double inferObservationCadenceSec(const std::vector<Observation>& observations) {
    std::vector<double> epochs;
    epochs.reserve(observations.size());
    for (const Observation& observation : observations) {
        epochs.push_back(observation.epochTdb);
    }

    std::sort(epochs.begin(), epochs.end());
    double cadenceSec = std::numeric_limits<double>::infinity();
    for (std::size_t i = 1; i < epochs.size(); ++i) {
        const double delta = epochs[i] - epochs[i - 1];
        if (delta > 1.0e-6) {
            cadenceSec = std::min(cadenceSec, delta);
        }
    }

    if (!std::isfinite(cadenceSec)) {
        throw std::runtime_error("Cannot infer observation cadence from the report.");
    }
    return cadenceSec;
}

void writePostfitDiagnosticsCsv(const std::filesystem::path& path,
                                const std::vector<Observation>& observations,
                                const Prediction& prefitPrediction,
                                const Prediction& postfitPrediction,
                                const od::EphemerisInterpolator& truthEphemeris,
                                const od::EphemerisInterpolator& estimatedEphemeris,
                                double referenceEpoch) {
    std::ofstream csv(path);
    if (!csv) {
        throw std::runtime_error("Failed to open post-fit diagnostics CSV: " + path.string());
    }

    csv << std::scientific << std::setprecision(12);
    csv << "utc,station,station_naif,epoch_tdb,dt_since_reference_s,"
           "range_observed_km,range_truth_km,range_prefit_km,range_postfit_km,"
           "prefit_range_residual_m,postfit_range_residual_m,range_sigma_km,"
           "range_rate_observed_km_s,range_rate_truth_km_s,range_rate_prefit_km_s,"
           "range_rate_postfit_km_s,prefit_range_rate_residual_mm_s,"
           "postfit_range_rate_residual_mm_s,range_rate_sigma_km_s,"
           "shapiro_delay_m,light_time_s,emit_epoch_tdb,"
           "true_x_km,true_y_km,true_z_km,true_vx_km_s,true_vy_km_s,true_vz_km_s,"
           "estimated_x_km,estimated_y_km,estimated_z_km,estimated_vx_km_s,"
           "estimated_vy_km_s,estimated_vz_km_s,"
           "error_x_km,error_y_km,error_z_km,error_vx_km_s,error_vy_km_s,error_vz_km_s,"
           "position_error_norm_km,velocity_error_norm_km_s\n";

    for (std::size_t i = 0; i < observations.size(); ++i) {
        const Observation& observation = observations[i];
        const Eigen::Index rangeRow = static_cast<Eigen::Index>(2 * i);
        const Eigen::Index rangeRateRow = static_cast<Eigen::Index>(2 * i + 1);
        const ComputedRange& postRange = postfitPrediction.ranges[i];
        const State6 truthAtReceive = interpolateState(truthEphemeris,
                                                       observation.epochTdb,
                                                       "Diagnostics truth interpolation");
        const State6 estimatedAtReceive = interpolateState(estimatedEphemeris,
                                                           observation.epochTdb,
                                                           "Diagnostics estimate interpolation");
        const State6 error = estimatedAtReceive - truthAtReceive;

        csv << observation.utc << ','
            << observation.stationName << ','
            << observation.stationNaif << ','
            << observation.epochTdb << ','
            << observation.epochTdb - referenceEpoch << ','
            << observation.rangeObservedKm << ','
            << observation.rangeTruthKm << ','
            << prefitPrediction.values[rangeRow] << ','
            << postfitPrediction.values[rangeRow] << ','
            << (observation.rangeObservedKm - prefitPrediction.values[rangeRow]) * 1000.0 << ','
            << (observation.rangeObservedKm - postfitPrediction.values[rangeRow]) * 1000.0 << ','
            << observation.rangeSigmaKm << ','
            << observation.rangeRateObservedKmPerSec << ','
            << observation.rangeRateTruthKmPerSec << ','
            << prefitPrediction.values[rangeRateRow] << ','
            << postfitPrediction.values[rangeRateRow] << ','
            << (observation.rangeRateObservedKmPerSec - prefitPrediction.values[rangeRateRow]) * 1.0e6 << ','
            << (observation.rangeRateObservedKmPerSec - postfitPrediction.values[rangeRateRow]) * 1.0e6 << ','
            << observation.rangeRateSigmaKmPerSec << ','
            << postRange.shapiroDelayKm * 1000.0 << ','
            << postRange.lightTimeSec << ','
            << postRange.emitEpochTdb;

        for (Eigen::Index component = 0; component < truthAtReceive.size(); ++component) {
            csv << ',' << truthAtReceive[component];
        }
        for (Eigen::Index component = 0; component < estimatedAtReceive.size(); ++component) {
            csv << ',' << estimatedAtReceive[component];
        }
        for (Eigen::Index component = 0; component < error.size(); ++component) {
            csv << ',' << error[component];
        }

        csv << ',' << error.segment<3>(0).norm()
            << ',' << error.segment<3>(3).norm()
            << '\n';
    }
}

void writeTrajectoryErrorCsv(const std::filesystem::path& path,
                             const std::vector<Observation>& observations,
                             const od::EphemerisInterpolator& truthEphemeris,
                             const od::EphemerisInterpolator& estimatedEphemeris,
                             double referenceEpoch) {
    if (observations.empty()) {
        throw std::invalid_argument("Cannot write trajectory diagnostics for an empty observation set.");
    }

    double startEpoch = observations.front().epochTdb;
    double finalEpoch = observations.front().epochTdb;
    for (const Observation& observation : observations) {
        startEpoch = std::min(startEpoch, observation.epochTdb);
        finalEpoch = std::max(finalEpoch, observation.epochTdb);
    }
    const double cadenceSec = inferObservationCadenceSec(observations);

    std::ofstream csv(path);
    if (!csv) {
        throw std::runtime_error("Failed to open trajectory error CSV: " + path.string());
    }

    csv << std::scientific << std::setprecision(12);
    csv << "utc,epoch_tdb,dt_since_reference_s,"
           "true_x_km,true_y_km,true_z_km,true_vx_km_s,true_vy_km_s,true_vz_km_s,"
           "estimated_x_km,estimated_y_km,estimated_z_km,estimated_vx_km_s,"
           "estimated_vy_km_s,estimated_vz_km_s,"
           "error_x_km,error_y_km,error_z_km,error_vx_km_s,error_vy_km_s,error_vz_km_s,"
           "position_error_norm_km,velocity_error_norm_km_s\n";

    for (double epoch = startEpoch; epoch <= finalEpoch + 1.0e-6; epoch += cadenceSec) {
        const State6 truth = interpolateState(truthEphemeris, epoch, "Trajectory truth interpolation");
        const State6 estimate = interpolateState(estimatedEphemeris, epoch, "Trajectory estimate interpolation");
        const State6 error = estimate - truth;

        csv << utcFromEt(epoch) << ','
            << epoch << ','
            << epoch - referenceEpoch;
        for (Eigen::Index component = 0; component < truth.size(); ++component) {
            csv << ',' << truth[component];
        }
        for (Eigen::Index component = 0; component < estimate.size(); ++component) {
            csv << ',' << estimate[component];
        }
        for (Eigen::Index component = 0; component < error.size(); ++component) {
            csv << ',' << error[component];
        }
        csv << ',' << error.segment<3>(0).norm()
            << ',' << error.segment<3>(3).norm()
            << '\n';
    }
}

void writeObservabilityWindowsCsv(const std::filesystem::path& path,
                                  const std::vector<Observation>& observations,
                                  double referenceEpoch) {
    if (observations.empty()) {
        throw std::invalid_argument("Cannot write observability windows for an empty observation set.");
    }

    const double cadenceSec = inferObservationCadenceSec(observations);
    const double windowBreakSec = 1.5 * cadenceSec;

    std::vector<Observation> sortedObservations = observations;
    std::sort(sortedObservations.begin(),
              sortedObservations.end(),
              [](const Observation& lhs, const Observation& rhs) {
                  if (lhs.stationName != rhs.stationName) {
                      return lhs.stationName < rhs.stationName;
                  }
                  return lhs.epochTdb < rhs.epochTdb;
              });

    std::ofstream csv(path);
    if (!csv) {
        throw std::runtime_error("Failed to open station observability windows CSV: " + path.string());
    }

    csv << std::scientific << std::setprecision(12);
    csv << "station,station_naif,window_index,start_utc,end_utc,start_epoch_tdb,end_epoch_tdb,"
           "start_dt_since_reference_s,end_dt_since_reference_s,duration_s,sample_count\n";

    std::size_t i = 0;
    while (i < sortedObservations.size()) {
        const std::string stationName = sortedObservations[i].stationName;
        const std::string stationNaif = sortedObservations[i].stationNaif;
        std::size_t windowIndex = 0;

        while (i < sortedObservations.size() && sortedObservations[i].stationName == stationName) {
            const double startEpoch = sortedObservations[i].epochTdb;
            double endEpoch = startEpoch;
            std::size_t sampleCount = 1;
            ++i;

            while (i < sortedObservations.size()
                   && sortedObservations[i].stationName == stationName
                   && sortedObservations[i].epochTdb - endEpoch <= windowBreakSec) {
                endEpoch = sortedObservations[i].epochTdb;
                ++sampleCount;
                ++i;
            }

            ++windowIndex;
            csv << stationName << ','
                << stationNaif << ','
                << windowIndex << ','
                << utcFromEt(startEpoch) << ','
                << utcFromEt(endEpoch) << ','
                << startEpoch << ','
                << endEpoch << ','
                << startEpoch - referenceEpoch << ','
                << endEpoch - referenceEpoch << ','
                << endEpoch - startEpoch << ','
                << sampleCount << '\n';
        }
    }
}

void writeReport(const std::filesystem::path& path,
                 const std::filesystem::path& observationPath,
                 const std::vector<Observation>& observations,
                 const std::vector<IterationRecord>& iterations,
                 const Prediction& prefitPrediction,
                 const Prediction& postfitPrediction,
                 const State6& truthState,
                 const od::EphemerisInterpolator& truthEphemeris,
                 const State6& priorState,
                 const State6& estimatedState,
                 const Eigen::MatrixXd& covariance,
                 const std::vector<std::string>& activeThirdBodies,
                 const std::vector<std::string>& skippedThirdBodies,
                 const od::RKF45Integrator& integrator,
                 const od::EphemerisInterpolator& estimatedEphemeris,
                 double referenceEpoch) {
    std::ofstream report(path);
    if (!report) {
        throw std::runtime_error("Failed to open Voyager OD report for writing: " + path.string());
    }

    const State6 priorError = priorState - truthState;
    const State6 estimatedError = estimatedState - truthState;
    const Eigen::VectorXd measured = measuredVector(observations);
    const ResidualStats prefitStats = computeStats(measured - prefitPrediction.values, observations);
    const ResidualStats postfitStats = computeStats(measured - postfitPrediction.values, observations);
    const std::vector<std::string> stationLabels = uniqueStationLabels(observations);

    report << std::scientific << std::setprecision(12);
    report << "# Voyager 1 multi-station OD test report\n"
           << "# input_observations: " << observationPath.string() << '\n'
           << "# stations:";
    for (const std::string& stationLabel : stationLabels) {
        report << ' ' << stationLabel;
    }
    report << '\n'
           << "# target: Voyager 1 (" << kTarget << ")\n"
           << "# reference_epoch_tdb: " << referenceEpoch << '\n'
           << "# dynamics: Sun central gravity + covered third bodies + cannonball SRP\n"
           << "# requested_third_body_ids: " << kJupiter << ' ' << kIo << ' ' << kEuropa
           << ' ' << kGanymede << ' ' << kCallisto << ' ' << kSaturnBarycenter << '\n'
           << "# active_third_body_ids:";
    for (const std::string& body : activeThirdBodies) {
        report << ' ' << body;
    }
    report << '\n'
           << "# skipped_third_body_ids_missing_spk:";
    for (const std::string& body : skippedThirdBodies) {
        report << ' ' << body;
    }
    report << '\n'
           << "# truth_model: initial CSPICE Voyager state propagated with the same RKF45 dynamics used by OD\n"
           << "# observation_model: Hermite-interpolated one-way light-time + solar Shapiro + finite-difference range-rate\n"
           << "# covariance_note: WLS uses sigma floors for simplified-dynamics model error; raw residuals below still use the input observations\n"
           << "# adopted_range_sigma_floor_km: " << kAdoptedRangeSigmaFloorKm << '\n'
           << "# adopted_range_rate_sigma_floor_km_s: " << kAdoptedRangeRateSigmaFloorKmPerSec << '\n'
           << "# SRP: CR=" << kVoyagerReflectivity
           << " area_m2=" << kVoyagerAreaM2
           << " mass_kg=" << kVoyagerMassKg << '\n'
           << "# integrator_abs_tol: " << integrator.options().absoluteTolerance << '\n'
           << "# integrator_rel_tol: " << integrator.options().relativeTolerance << '\n'
           << "# integrator_max_step_s: " << integrator.options().maximumStep << '\n'
           << "# prior_state_sun_j2000_km_km_s: " << priorState.transpose() << '\n'
           << "# estimated_state_sun_j2000_km_km_s: " << estimatedState.transpose() << '\n'
           << "# rkf45_truth_initial_state_sun_j2000_km_km_s: " << truthState.transpose() << '\n'
           << "# prior_position_error_km: " << priorError.segment<3>(0).norm() << '\n'
           << "# prior_velocity_error_km_s: " << priorError.segment<3>(3).norm() << '\n'
           << "# posterior_position_error_km: " << estimatedError.segment<3>(0).norm() << '\n'
           << "# posterior_velocity_error_km_s: " << estimatedError.segment<3>(3).norm() << '\n'
           << "# prefit_range_rms_m: " << prefitStats.rangeRmsKm * 1000.0 << '\n'
           << "# postfit_range_rms_m: " << postfitStats.rangeRmsKm * 1000.0 << '\n'
           << "# prefit_range_rate_rms_mm_s: " << prefitStats.rangeRateRmsKmPerSec * 1.0e6 << '\n'
           << "# postfit_range_rate_rms_mm_s: " << postfitStats.rangeRateRmsKmPerSec * 1.0e6 << '\n'
           << "# prefit_weighted_rms: " << prefitStats.weightedRms << '\n'
           << "# postfit_weighted_rms: " << postfitStats.weightedRms << '\n'
           << "# posterior_covariance_diag: " << covariance.diagonal().transpose() << '\n'
           << "# iterations: iteration weighted_rms range_rms_m range_rate_rms_mm_s "
              "correction_pos_m correction_vel_mm_s position_error_km velocity_error_km_s\n";

    for (const IterationRecord& record : iterations) {
        report << "# iter "
               << record.iteration << ' '
               << record.weightedRms << ' '
               << record.rangeRmsKm * 1000.0 << ' '
               << record.rangeRateRmsKmPerSec * 1.0e6 << ' '
               << record.correctionPositionNormKm * 1000.0 << ' '
               << record.correctionVelocityNormKmPerSec * 1.0e6 << ' '
               << record.positionErrorKm << ' '
               << record.velocityErrorKmPerSec << '\n';
    }

    report << "# columns: utc station epoch_tdb obs_range_km truth_range_km prefit_range_km "
              "postfit_range_km prefit_range_residual_m postfit_range_residual_m "
              "obs_range_rate_km_s truth_range_rate_km_s prefit_range_rate_km_s "
              "postfit_range_rate_km_s prefit_range_rate_residual_mm_s "
              "postfit_range_rate_residual_mm_s shapiro_delay_m light_time_s "
              "emit_epoch_tdb propagated_position_error_km propagated_velocity_error_km_s\n";

    for (std::size_t i = 0; i < observations.size(); ++i) {
        const Observation& observation = observations[i];
        const Eigen::Index rangeRow = static_cast<Eigen::Index>(2 * i);
        const Eigen::Index rangeRateRow = static_cast<Eigen::Index>(2 * i + 1);
        const ComputedRange& postRange = postfitPrediction.ranges[i];

        const State6 truthAtReceive = interpolateState(truthEphemeris,
                                                       observation.epochTdb,
                                                       "Truth report interpolation");
        const State6 estimatedAtReceive = interpolateState(estimatedEphemeris,
                                                           observation.epochTdb,
                                                           "Post-fit report interpolation");
        const State6 propagatedError = estimatedAtReceive - truthAtReceive;

        report << observation.utc << ' '
               << observation.stationName << ' '
               << observation.epochTdb << ' '
               << observation.rangeObservedKm << ' '
               << observation.rangeTruthKm << ' '
               << prefitPrediction.values[rangeRow] << ' '
               << postfitPrediction.values[rangeRow] << ' '
               << (observation.rangeObservedKm - prefitPrediction.values[rangeRow]) * 1000.0 << ' '
               << (observation.rangeObservedKm - postfitPrediction.values[rangeRow]) * 1000.0 << ' '
               << observation.rangeRateObservedKmPerSec << ' '
               << observation.rangeRateTruthKmPerSec << ' '
               << prefitPrediction.values[rangeRateRow] << ' '
               << postfitPrediction.values[rangeRateRow] << ' '
               << (observation.rangeRateObservedKmPerSec - prefitPrediction.values[rangeRateRow]) * 1.0e6 << ' '
               << (observation.rangeRateObservedKmPerSec - postfitPrediction.values[rangeRateRow]) * 1.0e6 << ' '
               << postRange.shapiroDelayKm * 1000.0 << ' '
               << postRange.lightTimeSec << ' '
               << postRange.emitEpochTdb << ' '
               << propagatedError.segment<3>(0).norm() << ' '
               << propagatedError.segment<3>(3).norm() << '\n';
    }
}

} // namespace

int main() {
    try {
        const SpiceErrorModeGuard spiceErrorMode;

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load Voyager OD meta-kernel");

        const std::filesystem::path observationPath = firstExistingPath(
            {
                "../synthetic_observations_dss43_voyager1.txt",
                "synthetic_observations_dss43_voyager1.txt",
                "../../synthetic_observations_dss43_voyager1.txt"
            },
            "synthetic Voyager observation report");
        const std::vector<Observation> observations = readObservationReport(observationPath);
        const double referenceEpoch = observations.front().epochTdb;
        const State6 truthState = spiceState(kTarget, referenceEpoch, kCentralBody, "NONE");

        State6 priorState = truthState;
        priorState << truthState[0] + 50.0,
                      truthState[1] - 40.0,
                      truthState[2] + 30.0,
                      truthState[3] + 5.0e-4,
                      truthState[4] - 4.0e-4,
                      truthState[5] + 2.0e-4;

        od::RKF45Integrator::Options options;
        options.absoluteTolerance = 1.0e-5;
        options.relativeTolerance = 1.0e-12;
        options.initialStep = 30.0;
        options.minimumStep = 1.0e-4;
        options.maximumStep = 300.0;
        options.maximumSteps = 200000;
        const od::RKF45Integrator integrator(options);

        fd::perturbations::ThirdBodyGravity thirdBody(kCentralBody, kFrame);
        std::vector<std::string> activeThirdBodies;
        std::vector<std::string> skippedThirdBodies;
        const auto [dynamicsStartEpoch, dynamicsFinalEpoch] =
            ephemerisBoundsForObservations(observations);

        addThirdBodyIfCovered(thirdBody,
                              {kJupiter, kJupiterBodyMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              true,
                              activeThirdBodies,
                              skippedThirdBodies);
        addThirdBodyIfCovered(thirdBody,
                              {kIo, kIoMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              false,
                              activeThirdBodies,
                              skippedThirdBodies);
        addThirdBodyIfCovered(thirdBody,
                              {kEuropa, kEuropaMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              false,
                              activeThirdBodies,
                              skippedThirdBodies);
        addThirdBodyIfCovered(thirdBody,
                              {kGanymede, kGanymedeMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              false,
                              activeThirdBodies,
                              skippedThirdBodies);
        addThirdBodyIfCovered(thirdBody,
                              {kCallisto, kCallistoMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              false,
                              activeThirdBodies,
                              skippedThirdBodies);
        addThirdBodyIfCovered(thirdBody,
                              {kSaturnBarycenter, kSaturnSystemMuKm3PerSec2},
                              dynamicsStartEpoch,
                              dynamicsFinalEpoch,
                              true,
                              activeThirdBodies,
                              skippedThirdBodies);

        fd::perturbations::SolarRadiationPressure srp(kCentralBody,
                                                      kFrame,
                                                      kVoyagerReflectivity,
                                                      kVoyagerAreaM2,
                                                      kVoyagerMassKg);
        const auto dynamics = makeDynamics(thirdBody, srp);
        const od::EphemerisInterpolator truthEphemeris =
            propagateObservationEphemeris(observations,
                                          integrator,
                                          dynamics,
                                          referenceEpoch,
                                          truthState);

        Eigen::MatrixXd priorCovariance = Eigen::MatrixXd::Zero(6, 6);
        priorCovariance.diagonal() << 100.0, 100.0, 100.0,
                                      1.0e-8, 1.0e-8, 1.0e-8;

        fd::filters::WLS filter(kCentralBody);
        filter.setInitialState(priorState, priorCovariance, referenceEpoch);

        const Eigen::VectorXd measured = measuredVector(observations);
        const Eigen::MatrixXd covariance = measurementCovariance(observations);
        const Prediction prefitPrediction = predictObservations(observations,
                                                                integrator,
                                                                dynamics,
                                                                referenceEpoch,
                                                                priorState);

        std::vector<IterationRecord> iterationRecords;
        for (int iteration = 0; iteration < kMaxSolverIterations; ++iteration) {
            const State6 currentState = stateAtEpoch(filter);
            const Prediction prediction = predictObservations(observations,
                                                              integrator,
                                                              dynamics,
                                                              referenceEpoch,
                                                              currentState);
            const Eigen::VectorXd residuals = measured - prediction.values;
            const ResidualStats stats = computeStats(residuals, observations);
            const Eigen::MatrixXd design = finiteDifferenceDesignMatrix(observations,
                                                                        integrator,
                                                                        dynamics,
                                                                        referenceEpoch,
                                                                        currentState);

            filter.processBatch(residuals, design, covariance);

            const State6 updatedState = stateAtEpoch(filter);
            const State6 correction = updatedState - currentState;
            const State6 updatedError = updatedState - truthState;

            iterationRecords.push_back(IterationRecord{
                iteration + 1,
                stats.weightedRms,
                stats.rangeRmsKm,
                stats.rangeRateRmsKmPerSec,
                correction.segment<3>(0).norm(),
                correction.segment<3>(3).norm(),
                updatedError.segment<3>(0).norm(),
                updatedError.segment<3>(3).norm()
            });

            if (correction.segment<3>(0).norm() < 1.0e-3
                && correction.segment<3>(3).norm() < 1.0e-8) {
                break;
            }
        }

        const State6 estimatedState = stateAtEpoch(filter);
        const Prediction postfitPrediction = predictObservations(observations,
                                                                 integrator,
                                                                 dynamics,
                                                                 referenceEpoch,
                                                                 estimatedState);
        const od::EphemerisInterpolator estimatedEphemeris =
            propagateObservationEphemeris(observations,
                                          integrator,
                                          dynamics,
                                          referenceEpoch,
                                          estimatedState);

        const Eigen::VectorXd prefitResiduals = measured - prefitPrediction.values;
        const Eigen::VectorXd postfitResiduals = measured - postfitPrediction.values;
        const ResidualStats prefitStats = computeStats(prefitResiduals, observations);
        const ResidualStats postfitStats = computeStats(postfitResiduals, observations);

        const double priorPositionErrorKm = (priorState - truthState).segment<3>(0).norm();
        const double posteriorPositionErrorKm = (estimatedState - truthState).segment<3>(0).norm();

        writeReport(kReportPath,
                    observationPath,
                    observations,
                    iterationRecords,
                    prefitPrediction,
                    postfitPrediction,
                    truthState,
                    truthEphemeris,
                    priorState,
                    estimatedState,
                    filter.covariance(),
                    activeThirdBodies,
                    skippedThirdBodies,
                    integrator,
                    estimatedEphemeris,
                    referenceEpoch);

        writePostfitDiagnosticsCsv(kPostfitDiagnosticsCsvPath,
                                   observations,
                                   prefitPrediction,
                                   postfitPrediction,
                                   truthEphemeris,
                                   estimatedEphemeris,
                                   referenceEpoch);
        writeTrajectoryErrorCsv(kTrajectoryErrorCsvPath,
                                observations,
                                truthEphemeris,
                                estimatedEphemeris,
                                referenceEpoch);
        writeObservabilityWindowsCsv(kObservabilityWindowsCsvPath,
                                     observations,
                                     referenceEpoch);

        std::cout << std::scientific << std::setprecision(9)
                  << "Voyager multi-station OD from synthetic report\n"
                  << "  observations                  : " << observations.size() << '\n'
                  << "  input report                  : " << observationPath.string() << '\n'
                  << "  output report                 : " << kReportPath << '\n'
                  << "  post-fit diagnostics CSV      : " << kPostfitDiagnosticsCsvPath << '\n'
                  << "  trajectory error CSV          : " << kTrajectoryErrorCsvPath << '\n'
                  << "  observability windows CSV     : " << kObservabilityWindowsCsvPath << '\n'
                  << "  active third bodies           : ";
        for (const std::string& body : activeThirdBodies) {
            std::cout << body << ' ';
        }
        std::cout << '\n'
                  << "  skipped third bodies          : ";
        for (const std::string& body : skippedThirdBodies) {
            std::cout << body << ' ';
        }
        std::cout << '\n'
                  << "  pre-fit range RMS             : " << prefitStats.rangeRmsKm * 1000.0 << " m\n"
                  << "  post-fit range RMS            : " << postfitStats.rangeRmsKm * 1000.0 << " m\n"
                  << "  pre-fit range-rate RMS        : " << prefitStats.rangeRateRmsKmPerSec * 1.0e6
                  << " mm/s\n"
                  << "  post-fit range-rate RMS       : " << postfitStats.rangeRateRmsKmPerSec * 1.0e6
                  << " mm/s\n"
                  << "  prior position error          : " << priorPositionErrorKm << " km\n"
                  << "  posterior position error      : " << posteriorPositionErrorKm << " km\n"
                  << "  posterior velocity error      : "
                  << (estimatedState - truthState).segment<3>(3).norm() << " km/s\n";

        if (!filter.covariance().allFinite()) {
            std::cerr << "FAIL: posterior covariance contains non-finite values.\n";
            return EXIT_FAILURE;
        }
        if (postfitStats.weightedRms >= prefitStats.weightedRms) {
            std::cerr << "FAIL: weighted residual RMS did not decrease.\n";
            return EXIT_FAILURE;
        }
        if (postfitStats.rangeRmsKm >= prefitStats.rangeRmsKm) {
            std::cerr << "FAIL: range residual RMS did not decrease.\n";
            return EXIT_FAILURE;
        }

        kclear_c();
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Voyager position OD test failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
