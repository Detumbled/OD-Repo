#pragma once

#include "observations/synth/obs_synth.hpp"

namespace fd::observations::synth {

class RangeRateSynth final : public SyntheticObservation {
public:
    explicit RangeRateSynth(GeometryConfig geometry = {}, NoiseConfig noise = {});

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds) override;
};

} // namespace fd::observations::synth
