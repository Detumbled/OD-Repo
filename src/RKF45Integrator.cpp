#include "RKF45Integrator.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace od {
namespace {

using State = RKF45Integrator::State;
using ConstStateRef = RKF45Integrator::ConstStateRef;
using DerivativeFunction = RKF45Integrator::DerivativeFunction;
using StateRef = RKF45Integrator::StateRef;
using StageMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

constexpr double kErrorExponent = -0.2;
constexpr double kZeroErrorScale = 5.0;

[[nodiscard]] bool isPositiveFinite(double value) noexcept {
    return std::isfinite(value) && value > 0.0;
}

void validateOptions(const RKF45Integrator::Options& options) {
    if (!isPositiveFinite(options.absoluteTolerance)) {
        throw std::invalid_argument("RKF45 absolute tolerance must be positive and finite.");
    }
    if (!isPositiveFinite(options.relativeTolerance)) {
        throw std::invalid_argument("RKF45 relative tolerance must be positive and finite.");
    }
    if (!isPositiveFinite(options.initialStep)) {
        throw std::invalid_argument("RKF45 initial step must be positive and finite.");
    }
    if (!isPositiveFinite(options.minimumStep)) {
        throw std::invalid_argument("RKF45 minimum step must be positive and finite.");
    }
    if (!isPositiveFinite(options.maximumStep)) {
        throw std::invalid_argument("RKF45 maximum step must be positive and finite.");
    }
    if (options.minimumStep > options.maximumStep) {
        throw std::invalid_argument("RKF45 minimum step cannot exceed maximum step.");
    }
    if (!isPositiveFinite(options.safetyFactor) || options.safetyFactor >= 1.0) {
        throw std::invalid_argument("RKF45 safety factor must be finite and in (0, 1).");
    }
    if (!isPositiveFinite(options.minimumScaleFactor) || !isPositiveFinite(options.maximumScaleFactor)
        || options.minimumScaleFactor > options.maximumScaleFactor) {
        throw std::invalid_argument("RKF45 scale factors must be positive, finite, and ordered.");
    }
    if (options.maximumSteps == 0) {
        throw std::invalid_argument("RKF45 maximum step count must be non-zero.");
    }
}

void validateState(const State& state, const char* label) {
    if (state.size() == 0) {
        throw std::invalid_argument(std::string(label) + " cannot be empty.");
    }

    if (!state.allFinite()) {
        throw std::invalid_argument(std::string(label) + " must contain only finite values.");
    }
}

void evaluateDynamics(const DerivativeFunction& dynamics,
                      double time,
                      ConstStateRef state,
                      StateRef derivative) {
    derivative.setZero();
    dynamics(time, state, derivative);

    if (derivative.size() != state.size()) {
        throw std::runtime_error("RKF45 dynamics returned a derivative with the wrong dimension.");
    }

    if (!derivative.allFinite()) {
        throw std::runtime_error("RKF45 dynamics returned a non-finite derivative.");
    }
}

[[nodiscard]] double normalizedError(ConstStateRef fourthOrder,
                                     ConstStateRef fifthOrder,
                                     ConstStateRef referenceState,
                                     const RKF45Integrator::Options& options) {
    return ((fifthOrder - fourthOrder).array()
            / (options.absoluteTolerance
               + options.relativeTolerance * referenceState.cwiseAbs().cwiseMax(fifthOrder.cwiseAbs()).array()))
        .matrix()
        .norm()
        / std::sqrt(static_cast<double>(fifthOrder.size()));
}

[[nodiscard]] double nextStepScale(double error, const RKF45Integrator::Options& options) {
    if (error <= std::numeric_limits<double>::epsilon()) {
        return std::min(kZeroErrorScale, options.maximumScaleFactor);
    }

    const double raw_scale = options.safetyFactor * std::pow(error, kErrorExponent);
    return std::clamp(raw_scale, options.minimumScaleFactor, options.maximumScaleFactor);
}

} // namespace

RKF45Integrator::RKF45Integrator()
    : RKF45Integrator(Options{}) {
}

RKF45Integrator::RKF45Integrator(Options options)
    : options_(options) {
    validateOptions(options_);
}

const RKF45Integrator::Options& RKF45Integrator::options() const noexcept {
    return options_;
}

RKF45Integrator::Result RKF45Integrator::integrate(double initialTime,
                                                   const State& initialState,
                                                   double finalTime,
                                                   const DerivativeFunction& dynamics) const {
    if (!std::isfinite(initialTime) || !std::isfinite(finalTime)) {
        throw std::invalid_argument("RKF45 integration times must be finite.");
    }
    if (!dynamics) {
        throw std::invalid_argument("RKF45 dynamics function is empty.");
    }

    validateState(initialState, "RKF45 initial state");

    const double direction = finalTime > initialTime ? 1.0 : -1.0;
    const double max_signed_step = direction * options_.maximumStep;
    const double min_signed_step = direction * options_.minimumStep;

    double time = initialTime;
    double step = direction * std::clamp(options_.initialStep, options_.minimumStep, options_.maximumStep);
    State state = initialState;
    State trial(initialState.size());
    State fourth_order(initialState.size());
    State fifth_order(initialState.size());
    State accepted_derivative(initialState.size());
    StageMatrix k(initialState.size(), 6);

    Result result;
    evaluateDynamics(dynamics, time, state, accepted_derivative);
    result.history.push_back({time, state, accepted_derivative});

    if (initialTime == finalTime) {
        result.state = initialState;
        result.finalTime = finalTime;
        return result;
    }

    while ((finalTime - time) * direction > 0.0) {
        if (result.acceptedSteps + result.rejectedSteps >= options_.maximumSteps) {
            throw std::runtime_error("RKF45 exceeded the maximum number of integration steps.");
        }

        const double remaining = finalTime - time;
        if (std::abs(step) > std::abs(remaining)) {
            step = remaining;
        }

        evaluateDynamics(dynamics, time, state, k.col(0));

        trial.noalias() = state + step * (1.0 / 4.0) * k.col(0);
        evaluateDynamics(dynamics, time + step * (1.0 / 4.0), trial, k.col(1));

        trial.noalias() = state + step * ((3.0 / 32.0) * k.col(0)
                                          + (9.0 / 32.0) * k.col(1));
        evaluateDynamics(dynamics, time + step * (3.0 / 8.0), trial, k.col(2));

        trial.noalias() = state + step * ((1932.0 / 2197.0) * k.col(0)
                                          + (-7200.0 / 2197.0) * k.col(1)
                                          + (7296.0 / 2197.0) * k.col(2));
        evaluateDynamics(dynamics, time + step * (12.0 / 13.0), trial, k.col(3));

        trial.noalias() = state + step * ((439.0 / 216.0) * k.col(0)
                                          - 8.0 * k.col(1)
                                          + (3680.0 / 513.0) * k.col(2)
                                          + (-845.0 / 4104.0) * k.col(3));
        evaluateDynamics(dynamics, time + step, trial, k.col(4));

        trial.noalias() = state + step * ((-8.0 / 27.0) * k.col(0)
                                          + 2.0 * k.col(1)
                                          + (-3544.0 / 2565.0) * k.col(2)
                                          + (1859.0 / 4104.0) * k.col(3)
                                          + (-11.0 / 40.0) * k.col(4));
        evaluateDynamics(dynamics, time + step * (1.0 / 2.0), trial, k.col(5));

        fourth_order.noalias() = state + step * ((25.0 / 216.0) * k.col(0)
                                                 + (1408.0 / 2565.0) * k.col(2)
                                                 + (2197.0 / 4104.0) * k.col(3)
                                                 - (1.0 / 5.0) * k.col(4));
        fifth_order.noalias() = state + step * ((16.0 / 135.0) * k.col(0)
                                                + (6656.0 / 12825.0) * k.col(2)
                                                + (28561.0 / 56430.0) * k.col(3)
                                                - (9.0 / 50.0) * k.col(4)
                                                + (2.0 / 55.0) * k.col(5));

        const double error = normalizedError(fourth_order, fifth_order, state, options_);
        const double scale = nextStepScale(error, options_);

        if (error <= 1.0) {
            state = fifth_order;
            time += step;
            ++result.acceptedSteps;
            evaluateDynamics(dynamics, time, state, accepted_derivative);
            result.history.push_back({time, state, accepted_derivative});
        } else {
            ++result.rejectedSteps;
        }

        const double candidate_step = step * scale;
        step = direction * std::clamp(std::abs(candidate_step), options_.minimumStep, options_.maximumStep);

        if (std::abs(step) <= options_.minimumStep
            && (finalTime - time) * direction > options_.minimumStep
            && error > 1.0) {
            throw std::runtime_error("RKF45 cannot satisfy tolerance before reaching the minimum step.");
        }

        if (std::abs(step) > std::abs(max_signed_step)) {
            step = max_signed_step;
        }
        if (std::abs(step) < std::abs(min_signed_step)) {
            step = min_signed_step;
        }
    }

    result.state = state;
    result.finalTime = finalTime;
    return result;
}

} // namespace od
