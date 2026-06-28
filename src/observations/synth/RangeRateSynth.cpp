#include "observations/synth/RangeRateSynth.hpp"

#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

constexpr double kShapiroDerivativeStepSeconds = 1.0;

} // namespace

RangeRateSynth::RangeRateSynth(GeometryConfig geometry, NoiseConfig noise)
    : SyntheticObservation(std::move(geometry), noise) {
}

std::vector<SyntheticObservationSample> RangeRateSynth::generate(double startTdb,
                                                                 double endTdb,
                                                                 double stepSeconds) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const RelativeGeometry geometry = relativeTargetGeometry(epoch);
        const Eigen::Vector3d relative_position = geometry.stationToTargetState.segment<3>(0);
        const Eigen::Vector3d relative_velocity = geometry.stationToTargetState.segment<3>(3);
        const double range = relative_position.norm();

        if (range <= 0.0) {
            throw std::runtime_error("Synthetic range-rate geometry has zero line-of-sight range.");
        }

        const double geometric_range_rate = relative_position.dot(relative_velocity) / range;
        const double shapiro_rate = (shapiroRangeDelay(epoch + kShapiroDerivativeStepSeconds)
                                     - shapiroRangeDelay(geometry))
            / kShapiroDerivativeStepSeconds;
        const double truth = geometric_range_rate + shapiro_rate;
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
