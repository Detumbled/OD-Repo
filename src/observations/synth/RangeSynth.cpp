#include "observations/synth/RangeSynth.hpp"

#include <stdexcept>
#include <utility>

namespace fd::observations::synth {

RangeSynth::RangeSynth(GeometryConfig geometry, NoiseConfig noise)
    : SyntheticObservation(std::move(geometry), noise) {
}

std::vector<SyntheticObservationSample> RangeSynth::generate(double startTdb,
                                                             double endTdb,
                                                             double stepSeconds) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const Eigen::Matrix<double, 6, 1> relative_state = relativeTargetState(epoch);
        const double truth = relative_state.segment<3>(0).norm();
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
