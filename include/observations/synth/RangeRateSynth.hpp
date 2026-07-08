#pragma once

#include "observations/synth/obs_synth.hpp"

namespace fd::observations::synth {

class RangeRateSynth final : public SyntheticObservation {
public:
    static constexpr double kDefaultCountTimeSeconds = 60.0;

    explicit RangeRateSynth(GeometryConfig geometry = {}, NoiseConfig noise = {});

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds) override;
    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   double countTimeSeconds);
};

} // namespace fd::observations::synth
