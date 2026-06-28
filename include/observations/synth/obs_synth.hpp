#pragma once

#include <Eigen/Dense>

#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace fd::observations::synth {

struct NoiseConfig {
    bool enabled {true};
    double sigma {0.0};
    std::uint32_t seed {5489U};
};

struct GeometryConfig {
    std::string target {"-31"};
    std::string stationName {"DSS-43"};
    std::string frame {"J2000"};
    std::string aberrationCorrection {"LT"};
};

struct SyntheticObservationSample {
    double epochTdb {0.0};
    double truth {0.0};
    double noise {0.0};
    double observed {0.0};
    double sigma {0.0};
};

struct RelativeGeometry {
    Eigen::Matrix<double, 6, 1> stationToTargetState {Eigen::Matrix<double, 6, 1>::Zero()};
    double lightTimeSec {0.0};
    double receiveEpochTdb {0.0};
    double emitEpochTdb {0.0};
};

class SyntheticObservation {
public:
    SyntheticObservation(GeometryConfig geometry, NoiseConfig noise);
    virtual ~SyntheticObservation();

    [[nodiscard]] virtual std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                           double endTdb,
                                                                           double stepSeconds) = 0;

    void setNoiseEnabled(bool enabled) noexcept;
    void setNoiseSigma(double sigma);
    void reseed(std::uint32_t seed);

    [[nodiscard]] const GeometryConfig& geometryConfig() const noexcept;
    [[nodiscard]] const NoiseConfig& noiseConfig() const noexcept;

    [[nodiscard]] static std::vector<double> makeEpochGrid(double startTdb,
                                                           double endTdb,
                                                           double stepSeconds);

protected:
    [[nodiscard]] Eigen::Matrix<double, 6, 1> relativeTargetState(double epochTdb) const;
    [[nodiscard]] RelativeGeometry relativeTargetGeometry(double receiveEpochTdb) const;
    [[nodiscard]] Eigen::Matrix<double, 6, 1> sunRelativeState(const std::string& body,
                                                              double epochTdb) const;
    [[nodiscard]] double shapiroRangeDelay(double receiveEpochTdb) const;
    [[nodiscard]] double shapiroRangeDelay(const RelativeGeometry& geometry) const;
    [[nodiscard]] double drawNoise();

    GeometryConfig geometry_;
    NoiseConfig noise_;
    std::mt19937 rng_;
};

} // namespace fd::observations::synth
