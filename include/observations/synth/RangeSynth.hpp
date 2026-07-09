#pragma once

#include "observations/synth/obs_synth.hpp"

namespace fd::observations::synth {

class RangeSynth final : public SyntheticObservation {
public:
    explicit RangeSynth(GeometryConfig geometry = {}, NoiseConfig noise = {});

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds) override;
    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   const TargetStateProvider& targetProvider);
};

} // namespace fd::observations::synth
