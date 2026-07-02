#include "observations/synth/VLBISynth.hpp"

#include "stations/ElevationMask.hpp"

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

const VLBIConfig& VLBISynth::vlbiConfig() const noexcept {
    return vlbi_;
}

std::vector<SyntheticObservationSample> VLBISynth::generateWithConfig(double startTdb,
                                                                      double endTdb,
                                                                      double stepSeconds,
                                                                      const VLBIConfig& vlbi) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const bool station_one_visible =
            od::isVisibleAboveElevation(vlbi.target,
                                        vlbi.stationOneName,
                                        epoch,
                                        vlbi.minimumElevationRad,
                                        vlbi.aberrationCorrection);
        const bool station_two_visible =
            od::isVisibleAboveElevation(vlbi.target,
                                        vlbi.stationTwoName,
                                        epoch,
                                        vlbi.minimumElevationRad,
                                        vlbi.aberrationCorrection);
        if (!station_one_visible || !station_two_visible) {
            continue;
        }

        const RelativeGeometry station_one_geometry =
            relativeTargetGeometryFor(vlbi.target, vlbi.stationOneName, epoch);
        const RelativeGeometry station_two_geometry =
            relativeTargetGeometryFor(vlbi.target, vlbi.stationTwoName, epoch);

        const double station_one_range =
            station_one_geometry.stationToTargetState.segment<3>(0).norm()
            + shapiroRangeDelayFor(vlbi.target, vlbi.stationOneName, station_one_geometry);
        const double station_two_range =
            station_two_geometry.stationToTargetState.segment<3>(0).norm()
            + shapiroRangeDelayFor(vlbi.target, vlbi.stationTwoName, station_two_geometry);

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
