#include "observations/synth/RangeRateSynth.hpp"

#include <stdexcept>
#include <utility>

namespace fd::observations::synth {

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
        const Eigen::Matrix<double, 6, 1> relative_state = relativeTargetState(epoch);
        const Eigen::Vector3d relative_position = relative_state.segment<3>(0);
        const Eigen::Vector3d relative_velocity = relative_state.segment<3>(3);
        const double range = relative_position.norm();

        if (range <= 0.0) {
            throw std::runtime_error("Synthetic range-rate geometry has zero line-of-sight range.");
        }

        const double truth = relative_position.dot(relative_velocity) / range;
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
