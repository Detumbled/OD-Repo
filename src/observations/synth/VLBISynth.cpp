#include "observations/synth/VLBISynth.hpp"

#include "stations/StationCatalog.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

[[nodiscard]] GeometryConfig makePrimaryGeometry(const VLBIConfig& vlbi) {
    GeometryConfig geometry;
    geometry.target = vlbi.target;
    geometry.stationName = vlbi.stationOneName;
    geometry.frame = vlbi.frame;
    geometry.aberrationCorrection = vlbi.aberrationCorrection;
    return geometry;
}

void validateVLBIConfig(const VLBIConfig& vlbi) {
    if (vlbi.target.empty()) {
        throw std::invalid_argument("VLBI synthetic observation target cannot be empty.");
    }
    if (vlbi.stationOneName.empty() || vlbi.stationTwoName.empty()) {
        throw std::invalid_argument("VLBI synthetic observation stations cannot be empty.");
    }
    if (vlbi.stationOneName == vlbi.stationTwoName) {
        throw std::invalid_argument("VLBI synthetic observation requires two distinct stations.");
    }
    if (vlbi.frame.empty()) {
        throw std::invalid_argument("VLBI synthetic observation frame cannot be empty.");
    }
    if (vlbi.aberrationCorrection.empty()) {
        throw std::invalid_argument("VLBI synthetic observation aberration correction cannot be empty.");
    }
    if (!std::isfinite(vlbi.minimumElevationRad)) {
        throw std::invalid_argument("VLBI synthetic observation elevation mask must be finite.");
    }
}

[[nodiscard]] od::Station stationFromKernel(const std::string& stationName, double epochTdb) {
    const int station_naif_id = od::stationNaifIdFromName(stationName);
    return od::buildStationFromKernel(stationName, station_naif_id, epochTdb);
}

} // namespace

VLBISynth::VLBISynth(VLBIConfig vlbi, NoiseConfig noise)
    : SyntheticObservation(makePrimaryGeometry(vlbi), noise),
      vlbi_(std::move(vlbi)) {
    validateVLBIConfig(vlbi_);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(double startTdb,
                                                            double endTdb,
                                                            double stepSeconds) {
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi_);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const TargetStateProvider& targetProvider) {
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi_, targetProvider);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(double startTdb,
                                                            double endTdb,
                                                            double stepSeconds,
                                                            const std::string& stationOneName,
                                                            const std::string& stationTwoName,
                                                            const std::string& target) {
    VLBIConfig vlbi = vlbi_;
    vlbi.stationOneName = stationOneName;
    vlbi.stationTwoName = stationTwoName;
    vlbi.target = target;
    validateVLBIConfig(vlbi);
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const std::string& stationOneName,
    const std::string& stationTwoName,
    const std::string& target,
    const TargetStateProvider& targetProvider) {
    VLBIConfig vlbi = vlbi_;
    vlbi.stationOneName = stationOneName;
    vlbi.stationTwoName = stationTwoName;
    vlbi.target = target;
    validateVLBIConfig(vlbi);
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi, targetProvider);
}

const VLBIConfig& VLBISynth::vlbiConfig() const noexcept {
    return vlbi_;
}

std::vector<SyntheticObservationSample> VLBISynth::generateWithConfig(double startTdb,
                                                                      double endTdb,
                                                                      double stepSeconds,
                                                                      const VLBIConfig& vlbi) {
    const SpiceTargetStateProvider target_provider(vlbi.target, vlbi.frame, "SUN");
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi, target_provider);
}

std::vector<SyntheticObservationSample> VLBISynth::generateWithConfig(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const VLBIConfig& vlbi,
    const TargetStateProvider& targetProvider) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);
    const od::Station station_one = stationFromKernel(vlbi.stationOneName, startTdb);
    const od::Station station_two = stationFromKernel(vlbi.stationTwoName, startTdb);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const bool station_one_visible =
            targetElevationRadFor(vlbi.stationOneName, station_one, epoch, targetProvider)
            >= vlbi.minimumElevationRad;
        const bool station_two_visible =
            targetElevationRadFor(vlbi.stationTwoName, station_two, epoch, targetProvider)
            >= vlbi.minimumElevationRad;
        if (!station_one_visible || !station_two_visible) {
            continue;
        }

        const RelativeGeometry station_one_geometry =
            relativeTargetGeometryFor(vlbi.stationOneName, epoch, targetProvider);
        const RelativeGeometry station_two_geometry =
            relativeTargetGeometryFor(vlbi.stationTwoName, epoch, targetProvider);

        const double station_one_range =
            station_one_geometry.stationToTargetState.segment<3>(0).norm()
            + shapiroRangeDelayFor(vlbi.stationOneName, station_one_geometry, targetProvider);
        const double station_two_range =
            station_two_geometry.stationToTargetState.segment<3>(0).norm()
            + shapiroRangeDelayFor(vlbi.stationTwoName, station_two_geometry, targetProvider);

        const double truth = station_two_range - station_one_range;
        const double noise = drawNoise();

        samples.push_back(SyntheticObservationSample{
            epoch,
            truth,
            noise,
            truth + noise,
            noiseConfig().sigma
        });
    }

    return samples;
}

} // namespace fd::observations::synth
