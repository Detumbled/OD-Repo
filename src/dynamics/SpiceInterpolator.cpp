#include "dynamics/SpiceInterpolator.hpp"

#include "utils/CSPICE/SpiceError.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace od {
namespace {

[[nodiscard]] bool positiveFinite(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

void validateInitializeInputs(const std::string& target,
                              const std::string& center,
                              const std::string& frame,
                              double startEt,
                              double endEt,
                              double stepSec,
                              const std::string& aberrationCorrection) {
    if (target.empty()) {
        throw std::invalid_argument("SPICE interpolation target cannot be empty.");
    }
    if (center.empty()) {
        throw std::invalid_argument("SPICE interpolation center cannot be empty.");
    }
    if (frame.empty()) {
        throw std::invalid_argument("SPICE interpolation frame cannot be empty.");
    }
    if (aberrationCorrection.empty()) {
        throw std::invalid_argument("SPICE interpolation aberration correction cannot be empty.");
    }
    if (!std::isfinite(startEt) || !std::isfinite(endEt)) {
        throw std::invalid_argument("SPICE interpolation bounds must be finite.");
    }
    if (!(endEt > startEt)) {
        throw std::invalid_argument("SPICE interpolation end time must be greater than start time.");
    }
    if (!positiveFinite(stepSec)) {
        throw std::invalid_argument("SPICE interpolation step must be positive and finite.");
    }
}

} // namespace

void SpiceInterpolator::initialize(const std::string& target,
                                   const std::string& center,
                                   const std::string& frame,
                                   double startEt,
                                   double endEt,
                                   double stepSec,
                                   const std::string& aberrationCorrection) {
    validateInitializeInputs(target, center, frame, startEt, endEt, stepSec, aberrationCorrection);

    const std::size_t regularIntervals =
        static_cast<std::size_t>(std::floor((endEt - startEt) / stepSec));
    const std::size_t reserveCount = regularIntervals + 2U;
    std::vector<double> sampledTimes;
    Vector3List sampledPositions;
    Vector3List sampledVelocities;
    sampledTimes.reserve(reserveCount);
    sampledPositions.reserve(reserveCount);
    sampledVelocities.reserve(reserveCount);

    SpiceErrorModeGuard actionGuard;

    for (std::size_t i = 0; i <= regularIntervals; ++i) {
        const double et = startEt + static_cast<double>(i) * stepSec;
        if (et > endEt) {
            break;
        }

        SpiceDouble spiceState[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        SpiceDouble lightTime = 0.0;
        spkezr_c(target.c_str(),
                 et,
                 frame.c_str(),
                 aberrationCorrection.c_str(),
                 center.c_str(),
                 spiceState,
                 &lightTime);
        throwIfSpiceFailed("Failed to sample SPICE state for " + target);

        sampledTimes.push_back(et);
        sampledPositions.emplace_back(spiceState[0], spiceState[1], spiceState[2]);
        sampledVelocities.emplace_back(spiceState[3], spiceState[4], spiceState[5]);
    }

    if (sampledTimes.empty() || sampledTimes.back() < endEt) {
        SpiceDouble spiceState[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        SpiceDouble lightTime = 0.0;
        spkezr_c(target.c_str(),
                 endEt,
                 frame.c_str(),
                 aberrationCorrection.c_str(),
                 center.c_str(),
                 spiceState,
                 &lightTime);
        throwIfSpiceFailed("Failed to sample final SPICE state for " + target);

        sampledTimes.push_back(endEt);
        sampledPositions.emplace_back(spiceState[0], spiceState[1], spiceState[2]);
        sampledVelocities.emplace_back(spiceState[3], spiceState[4], spiceState[5]);
    }

    if (sampledTimes.size() < 2) {
        throw std::runtime_error("SPICE interpolation requires at least two sampled states.");
    }

    target_ = target;
    center_ = center;
    frame_ = frame;
    aberrationCorrection_ = aberrationCorrection;
    startEt_ = startEt;
    endEt_ = endEt;
    stepSec_ = stepSec;
    times_ = std::move(sampledTimes);
    positions_ = std::move(sampledPositions);
    velocities_ = std::move(sampledVelocities);
}

SpiceInterpolator::Vector3 SpiceInterpolator::getPosition(double et) const {
    if (times_.size() < 2) {
        throw std::runtime_error("Cannot interpolate an uninitialized SPICE ephemeris.");
    }
    if (et == startEt_) {
        return positions_.front();
    }
    if (et == endEt_) {
        return positions_.back();
    }

    return interpolateState(segmentIndex(et), et).head<3>();
}

SpiceInterpolator::Vector3 SpiceInterpolator::getVelocity(double et) const {
    if (times_.size() < 2) {
        throw std::runtime_error("Cannot interpolate an uninitialized SPICE ephemeris.");
    }
    if (et == startEt_) {
        return velocities_.front();
    }
    if (et == endEt_) {
        return velocities_.back();
    }

    return interpolateState(segmentIndex(et), et).tail<3>();
}

SpiceInterpolator::State6 SpiceInterpolator::getState(double et) const {
    if (times_.size() < 2) {
        throw std::runtime_error("Cannot interpolate an uninitialized SPICE ephemeris.");
    }
    if (et == startEt_) {
        State6 state;
        state.head<3>() = positions_.front();
        state.tail<3>() = velocities_.front();
        return state;
    }
    if (et == endEt_) {
        State6 state;
        state.head<3>() = positions_.back();
        state.tail<3>() = velocities_.back();
        return state;
    }

    return interpolateState(segmentIndex(et), et);
}

bool SpiceInterpolator::empty() const noexcept {
    return times_.empty();
}

std::size_t SpiceInterpolator::size() const noexcept {
    return times_.size();
}

double SpiceInterpolator::startEt() const noexcept {
    return startEt_;
}

double SpiceInterpolator::endEt() const noexcept {
    return endEt_;
}

double SpiceInterpolator::stepSec() const noexcept {
    return stepSec_;
}

std::size_t SpiceInterpolator::segmentIndex(double et) const {
    if (!std::isfinite(et)) {
        throw std::invalid_argument("SPICE interpolation time must be finite.");
    }
    if (times_.size() < 2) {
        throw std::runtime_error("Cannot interpolate an uninitialized SPICE ephemeris.");
    }
    if (et < startEt_ || et > endEt_) {
        throw std::out_of_range("SPICE interpolation time is outside the cached ephemeris span.");
    }

    std::size_t index = static_cast<std::size_t>(std::floor((et - startEt_) / stepSec_));
    if (index + 1U >= times_.size()) {
        index = times_.size() - 2U;
    }

    if (et < times_[index] || et > times_[index + 1U]) {
        const auto upper = std::upper_bound(times_.begin(), times_.end(), et);
        if (upper == times_.begin() || upper == times_.end()) {
            throw std::runtime_error("Failed to locate SPICE interpolation segment.");
        }
        index = static_cast<std::size_t>(std::distance(times_.begin(), upper) - 1);
    }

    return index;
}

SpiceInterpolator::State6 SpiceInterpolator::interpolateState(std::size_t segment, double et) const {
    const double t0 = times_[segment];
    const double t1 = times_[segment + 1U];
    const double dt = t1 - t0;
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::runtime_error("SPICE interpolation segment has invalid duration.");
    }

    const double s = (et - t0) / dt;
    const double s2 = s * s;
    const double s3 = s2 * s;

    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + s;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;

    const Vector3& p0 = positions_[segment];
    const Vector3& p1 = positions_[segment + 1U];
    const Vector3& v0 = velocities_[segment];
    const Vector3& v1 = velocities_[segment + 1U];

    State6 state;
    state.head<3>() = h00 * p0 + h10 * dt * v0 + h01 * p1 + h11 * dt * v1;

    const double dh00 = (6.0 * s2 - 6.0 * s) / dt;
    const double dh10 = 3.0 * s2 - 4.0 * s + 1.0;
    const double dh01 = (-6.0 * s2 + 6.0 * s) / dt;
    const double dh11 = 3.0 * s2 - 2.0 * s;
    state.tail<3>() = dh00 * p0 + dh10 * v0 + dh01 * p1 + dh11 * v1;

    return state;
}

} // namespace od
