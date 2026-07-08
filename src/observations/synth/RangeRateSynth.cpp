#include "observations/synth/RangeRateSynth.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

void validateCountTime(double countTimeSeconds) {
    if (!std::isfinite(countTimeSeconds) || countTimeSeconds <= 0.0) {
        throw std::invalid_argument("Synthetic range-rate count time must be positive and finite.");
    }
}

} // namespace

RangeRateSynth::RangeRateSynth(GeometryConfig geometry, NoiseConfig noise)
    : SyntheticObservation(std::move(geometry), noise) {
}

std::vector<SyntheticObservationSample> RangeRateSynth::generate(double startTdb,
                                                                 double endTdb,
                                                                 double stepSeconds) {
    return generate(startTdb, endTdb, stepSeconds, kDefaultCountTimeSeconds);
}

std::vector<SyntheticObservationSample> RangeRateSynth::generate(double startTdb,
                                                                 double endTdb,
                                                                 double stepSeconds,
                                                                 double countTimeSeconds) {
    validateCountTime(countTimeSeconds);
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);
    const double half_count_time = 0.5 * countTimeSeconds;
    const double scaled_sigma = noiseConfig().sigma / std::sqrt(countTimeSeconds);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    for (const double epoch : epochs) {
        const RelativeGeometry start_geometry = relativeTargetGeometry(epoch - half_count_time);
        const RelativeGeometry end_geometry = relativeTargetGeometry(epoch + half_count_time);
        const double start_range =
            start_geometry.stationToTargetState.segment<3>(0).norm() + shapiroRangeDelay(start_geometry);
        const double end_range =
            end_geometry.stationToTargetState.segment<3>(0).norm() + shapiroRangeDelay(end_geometry);

        if (start_range <= 0.0 || end_range <= 0.0) {
            throw std::runtime_error("Synthetic range-rate geometry has zero line-of-sight range.");
        }

        const double truth = (end_range - start_range) / countTimeSeconds;
        const double noise = drawNoise(scaled_sigma);

        samples.push_back(SyntheticObservationSample{
            epoch,
            truth,
            noise,
            truth + noise,
            scaled_sigma
        });
    }

    return samples;
}

} // namespace fd::observations::synth
