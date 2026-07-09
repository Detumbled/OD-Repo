#pragma once

#include "dynamics/EphemerisInterpolator.hpp"

#include <Eigen/Dense>

#include <string>

namespace fd::observations::synth {

using TargetState6 = Eigen::Matrix<double, 6, 1>;

class TargetStateProvider {
public:
    virtual ~TargetStateProvider() = default;

    // State must be expressed in geometry.frame relative to the Sun.
    [[nodiscard]] virtual TargetState6 stateAt(double epochTdb) const = 0;
};

class InterpolatedTargetStateProvider final : public TargetStateProvider {
public:
    explicit InterpolatedTargetStateProvider(const od::EphemerisInterpolator& ephemeris,
                                             std::string context = "Target ephemeris interpolation");

    [[nodiscard]] TargetState6 stateAt(double epochTdb) const override;

private:
    const od::EphemerisInterpolator* ephemeris_ {nullptr};
    std::string context_;
};

class SpiceTargetStateProvider final : public TargetStateProvider {
public:
    SpiceTargetStateProvider(std::string target,
                             std::string frame = "J2000",
                             std::string observer = "SUN");

    [[nodiscard]] TargetState6 stateAt(double epochTdb) const override;

private:
    std::string target_;
    std::string frame_;
    std::string observer_;
};

} // namespace fd::observations::synth
