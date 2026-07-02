#pragma once

#include "observations/synth/obs_synth.hpp"

#include <string>

namespace fd::observations::synth {

struct VLBIConfig {
    std::string target {"-31"};
    std::string stationOneName {"DSS-43"};
    std::string stationTwoName {"DSS-63"};
    std::string frame {"J2000"};
    std::string aberrationCorrection {"LT"};
    double minimumElevationRad {0.0};
};

class VLBISynth final : public SyntheticObservation {
public:
    explicit VLBISynth(VLBIConfig vlbi = {}, NoiseConfig noise = {});

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds) override;

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   const std::string& stationOneName,
                                                                   const std::string& stationTwoName,
                                                                   const std::string& target);

    [[nodiscard]] const VLBIConfig& vlbiConfig() const noexcept;

private:
    [[nodiscard]] std::vector<SyntheticObservationSample> generateWithConfig(double startTdb,
                                                                             double endTdb,
                                                                             double stepSeconds,
                                                                             const VLBIConfig& vlbi);

    VLBIConfig vlbi_;
};

} // namespace fd::observations::synth
