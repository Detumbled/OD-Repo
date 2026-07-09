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
        const RelativeGeometry geometry = relativeTargetGeometry(epoch);
        const double geometric_range = geometry.stationToTargetState.segment<3>(0).norm();
        const double truth = geometric_range + shapiroRangeDelay(geometry);
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

std::vector<SyntheticObservationSample> RangeSynth::generate(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const TargetStateProvider& targetProvider) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const RelativeGeometry geometry = relativeTargetGeometry(epoch, targetProvider);
        const double geometric_range = geometry.stationToTargetState.segment<3>(0).norm();
        const double truth =
            geometric_range + shapiroRangeDelayFor(geometryConfig().stationName,
                                                   geometry,
                                                   targetProvider);
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
