#ifndef DP853_INTEGRATOR_HPP
#define DP853_INTEGRATOR_HPP

#include "dynamics/EphemerisInterpolator.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace od {

class DP853Integrator {
public:
    using State6 = Eigen::Matrix<double, 6, 1>;
    using DynamicsFunction = std::function<State6(double, State6)>;

    struct Options {
        double absoluteTolerance {1.0e-10};
        double relativeTolerance {1.0e-10};
        double minimumStep {1.0e-6};
        double maximumStep {86400.0};
        double safetyFactor {0.9};
        double minimumScaleFactor {0.2};
        double maximumScaleFactor {10.0};
        double proportionalExponent {0.7 / 8.0};
        double integralExponent {0.4 / 8.0};
        std::size_t maximumSteps {1000000};
    };

    struct Result {
        using EphemerisNode = EphemerisInterpolator::Node;

        State6 state {State6::Zero()};
        double finalTime {0.0};
        std::size_t acceptedSteps {0};
        std::size_t rejectedSteps {0};
        std::vector<EphemerisNode> history;
    };

    DP853Integrator() = default;

    explicit DP853Integrator(Options options)
        : options_(options) {
        validateOptions(options_);
    }

    [[nodiscard]] const Options& options() const noexcept {
        return options_;
    }

    [[nodiscard]] State6 integrate(double t0,
                                   State6 y0,
                                   double tEnd,
                                   DynamicsFunction dynamics,
                                   double& initialStep) const {
        return integrateInternal(t0, y0, tEnd, std::move(dynamics), initialStep, false).state;
    }

    [[nodiscard]] Result integrateWithHistory(double t0,
                                              State6 y0,
                                              double tEnd,
                                              DynamicsFunction dynamics,
                                              double& initialStep) const {
        return integrateInternal(t0, y0, tEnd, std::move(dynamics), initialStep, true);
    }

private:
    struct StepResult {
        State6 state {State6::Zero()};
        State6 derivative {State6::Zero()};
        State6 err5 {State6::Zero()};
        State6 err3 {State6::Zero()};
    };

    static constexpr std::size_t kStageCount = 13;
    using StageArray = std::array<State6, kStageCount>;

    [[nodiscard]] Result integrateInternal(double t0,
                                           State6 y0,
                                           double tEnd,
                                           DynamicsFunction dynamics,
                                           double& initialStep,
                                           bool collectHistory) const {
        validateOptions(options_);
        if (!std::isfinite(t0) || !std::isfinite(tEnd)) {
            throw std::invalid_argument("DP853 integration times must be finite.");
        }
        if (!y0.allFinite()) {
            throw std::invalid_argument("DP853 initial state must be finite.");
        }
        if (!dynamics) {
            throw std::invalid_argument("DP853 dynamics function is empty.");
        }
        if (!positiveFinite(initialStep)) {
            throw std::invalid_argument("DP853 initial step must be positive and finite.");
        }

        const double direction = tEnd > t0 ? 1.0 : -1.0;
        double hAbs = std::clamp(std::abs(initialStep), options_.minimumStep, options_.maximumStep);
        double time = t0;
        State6 state = y0;
        State6 derivative = evaluateDynamics(dynamics, time, state);
        Result result;
        result.state = state;
        result.finalTime = t0;
        if (collectHistory) {
            appendHistoryNode(result, time, state, derivative);
        }

        if (t0 == tEnd) {
            return result;
        }

        double previousError = 1.0;
        bool previousStepRejected = false;

        while ((tEnd - time) * direction > 0.0) {
            if (result.acceptedSteps + result.rejectedSteps >= options_.maximumSteps) {
                throw std::runtime_error("DP853 exceeded the maximum number of integration steps.");
            }

            const double remainingAbs = std::abs(tEnd - time);
            if (hAbs > remainingAbs) {
                hAbs = remainingAbs;
            }

            const double h = direction * hAbs;
            StepResult step = computeStep(dynamics, time, state, derivative, h);
            const double error = normalizedError(state, step.state, step.err5, step.err3, h);

            if (error <= 1.0) {
                time += h;
                state = step.state;
                derivative = step.derivative;
                ++result.acceptedSteps;
                if (collectHistory) {
                    appendHistoryNode(result, time, state, derivative);
                }

                double factor = nextStepFactor(error, previousError);
                if (previousStepRejected) {
                    factor = std::min(1.0, factor);
                }
                hAbs = std::clamp(hAbs * factor, options_.minimumStep, options_.maximumStep);
                previousError = std::max(error, 1.0e-4);
                previousStepRejected = false;
            } else {
                ++result.rejectedSteps;
                previousStepRejected = true;
                hAbs = std::clamp(hAbs * rejectionStepFactor(error),
                                  options_.minimumStep,
                                  options_.maximumStep);

                if (hAbs <= options_.minimumStep
                    && std::abs(tEnd - time) > options_.minimumStep) {
                    throw std::runtime_error("DP853 cannot satisfy tolerance before reaching the minimum step.");
                }
            }
        }

        initialStep = hAbs;
        result.state = state;
        result.finalTime = tEnd;
        return result;
    }

    static constexpr std::array<double, 12> kC {{
        0.0,
        0.526001519587677318785587544488e-01,
        0.789002279381515978178381316732e-01,
        0.118350341907227396726757197510,
        0.281649658092772603273242802490,
        0.333333333333333333333333333333,
        0.25,
        0.307692307692307692307692307692,
        0.651282051282051282051282051282,
        0.6,
        0.857142857142857142857142857142,
        1.0
    }};

    static constexpr std::array<std::array<double, 12>, 12> kA {{
        {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{5.26001519587677318785587544488e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{1.97250569845378994544595329183e-2, 5.91751709536136983633785987549e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{2.95875854768068491816892993775e-2, 0.0, 8.87627564304205475450678981324e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{2.41365134159266685502369798665e-1, 0.0, -8.84549479328286085344864962717e-1, 9.24834003261792003115737966543e-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{3.7037037037037037037037037037e-2, 0.0, 0.0, 1.70828608729473871279604482173e-1, 1.25467687566822425016691814123e-1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{3.7109375e-2, 0.0, 0.0, 1.70252211019544039314978060272e-1, 6.02165389804559606850219397283e-2, -1.7578125e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{3.70920001185047927108779319836e-2, 0.0, 0.0, 1.70383925712239993810214054705e-1, 1.07262030446373284651809199168e-1, -1.53194377486244017527936158236e-2, 8.27378916381402288758473766002e-3, 0.0, 0.0, 0.0, 0.0, 0.0}},
        {{6.24110958716075717114429577812e-1, 0.0, 0.0, -3.36089262944694129406857109825, -8.68219346841726006818189891453e-1, 2.75920996994467083049415600797e1, 2.01540675504778934086186788979e1, -4.34898841810699588477366255144e1, 0.0, 0.0, 0.0, 0.0}},
        {{4.77662536438264365890433908527e-1, 0.0, 0.0, -2.48811461997166764192642586468, -5.90290826836842996371446475743e-1, 2.12300514481811942347288949897e1, 1.52792336328824235832596922938e1, -3.32882109689848629194453265587e1, -2.03312017085086261358222928593e-2, 0.0, 0.0, 0.0}},
        {{-9.3714243008598732571704021658e-1, 0.0, 0.0, 5.18637242884406370830023853209, 1.09143734899672957818500254654, -8.14978701074692612513997267357, -1.85200656599969598641566180701e1, 2.27394870993505042818970056734e1, 2.49360555267965238987089396762, -3.0467644718982195003823669022, 0.0, 0.0}},
        {{2.27331014751653820792359768449, 0.0, 0.0, -1.05344954667372501984066689879e1, -2.00087205822486249909675718444, -1.79589318631187989172765950534e1, 2.79488845294199600508499808837e1, -2.85899827713502369474065508674, -8.87285693353062954433549289258, 1.23605671757943030647266201528e1, 6.43392746015763530355970484046e-1, 0.0}}
    }};

    static constexpr std::array<double, 12> kB {{
        5.42937341165687622380535766363e-2,
        0.0,
        0.0,
        0.0,
        0.0,
        4.45031289275240888144113950566,
        1.89151789931450038304281599044,
        -5.8012039600105847814672114227,
        3.1116436695781989440891606237e-1,
        -1.52160949662516078556178806805e-1,
        2.01365400804030348374776537501e-1,
        4.47106157277725905176885569043e-2
    }};

    static constexpr std::array<double, kStageCount> kE3 {{
        kB[0] - 0.244094488188976377952755905512,
        0.0,
        0.0,
        0.0,
        0.0,
        kB[5],
        kB[6],
        kB[7],
        kB[8] - 0.733846688281611857341361741547,
        kB[9],
        kB[10],
        kB[11] - 0.220588235294117647058823529412e-1,
        0.0
    }};

    static constexpr std::array<double, kStageCount> kE5 {{
        0.1312004499419488073250102996e-1,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.1225156446376204440720569753e+1,
        -0.4957589496572501915214079952,
        0.1664377182454986536961530415e+1,
        -0.3503288487499736816886487290,
        0.3341791187130174790297318841,
        0.8192320648511571246570742613e-1,
        -0.2235530786388629525884427845e-1,
        0.0
    }};

    static bool positiveFinite(double value) noexcept {
        return std::isfinite(value) && value > 0.0;
    }

    static void validateOptions(const Options& options) {
        if (!positiveFinite(options.absoluteTolerance)) {
            throw std::invalid_argument("DP853 absolute tolerance must be positive and finite.");
        }
        if (!positiveFinite(options.relativeTolerance)) {
            throw std::invalid_argument("DP853 relative tolerance must be positive and finite.");
        }
        if (!positiveFinite(options.minimumStep)) {
            throw std::invalid_argument("DP853 minimum step must be positive and finite.");
        }
        if (!positiveFinite(options.maximumStep)) {
            throw std::invalid_argument("DP853 maximum step must be positive and finite.");
        }
        if (options.minimumStep > options.maximumStep) {
            throw std::invalid_argument("DP853 minimum step cannot exceed maximum step.");
        }
        if (!positiveFinite(options.safetyFactor) || options.safetyFactor >= 1.0) {
            throw std::invalid_argument("DP853 safety factor must be finite and in (0, 1).");
        }
        if (!positiveFinite(options.minimumScaleFactor)
            || !positiveFinite(options.maximumScaleFactor)
            || options.minimumScaleFactor > options.maximumScaleFactor) {
            throw std::invalid_argument("DP853 step scale factors must be positive, finite, and ordered.");
        }
        if (!positiveFinite(options.proportionalExponent) || !positiveFinite(options.integralExponent)) {
            throw std::invalid_argument("DP853 PI exponents must be positive and finite.");
        }
        if (options.maximumSteps == 0U) {
            throw std::invalid_argument("DP853 maximum step count must be non-zero.");
        }
    }

    [[nodiscard]] static State6 evaluateDynamics(const DynamicsFunction& dynamics,
                                                 double time,
                                                 const State6& state) {
        State6 derivative = dynamics(time, state);
        if (!derivative.allFinite()) {
            throw std::runtime_error("DP853 dynamics returned a non-finite derivative.");
        }
        return derivative;
    }

    static void appendHistoryNode(Result& result,
                                  double time,
                                  const State6& state,
                                  const State6& derivative) {
        result.history.push_back({time, Eigen::VectorXd(state), Eigen::VectorXd(derivative)});
    }

    [[nodiscard]] StepResult computeStep(const DynamicsFunction& dynamics,
                                         double time,
                                         const State6& state,
                                         const State6& derivative,
                                         double h) const {
        StageArray k {};
        k[0] = derivative;

        for (std::size_t stage = 1; stage < 12; ++stage) {
            State6 increment = State6::Zero();
            for (std::size_t previous = 0; previous < stage; ++previous) {
                increment += kA[stage][previous] * k[previous];
            }

            const State6 stageState = state + h * increment;
            k[stage] = evaluateDynamics(dynamics, time + kC[stage] * h, stageState);
        }

        State6 yNext = state;
        for (std::size_t stage = 0; stage < 12; ++stage) {
            yNext += h * kB[stage] * k[stage];
        }
        k[12] = evaluateDynamics(dynamics, time + h, yNext);

        StepResult result;
        result.state = yNext;
        result.derivative = k[12];
        for (std::size_t stage = 0; stage < kStageCount; ++stage) {
            result.err5 += kE5[stage] * k[stage];
            result.err3 += kE3[stage] * k[stage];
        }

        return result;
    }

    [[nodiscard]] double normalizedError(const State6& currentState,
                                         const State6& nextState,
                                         const State6& err5,
                                         const State6& err3,
                                         double h) const {
        State6 err5Scaled;
        State6 err3Scaled;
        for (Eigen::Index i = 0; i < State6::RowsAtCompileTime; ++i) {
            const double scale = options_.absoluteTolerance
                + options_.relativeTolerance * std::max(std::abs(currentState[i]), std::abs(nextState[i]));
            err5Scaled[i] = err5[i] / scale;
            err3Scaled[i] = err3[i] / scale;
        }

        const double err5NormSquared = err5Scaled.squaredNorm();
        const double err3NormSquared = err3Scaled.squaredNorm();
        if (err5NormSquared == 0.0 && err3NormSquared == 0.0) {
            return 0.0;
        }

        const double denominator = err5NormSquared + 0.01 * err3NormSquared;
        return std::abs(h) * err5NormSquared
            / std::sqrt(denominator * static_cast<double>(State6::RowsAtCompileTime));
    }

    [[nodiscard]] double nextStepFactor(double error, double previousError) const {
        if (error <= std::numeric_limits<double>::epsilon()) {
            return options_.maximumScaleFactor;
        }

        const double boundedPrevious = std::max(previousError, 1.0e-4);
        const double rawFactor = options_.safetyFactor
            * std::pow(error, -options_.proportionalExponent)
            * std::pow(boundedPrevious, options_.integralExponent);
        return std::clamp(rawFactor, options_.minimumScaleFactor, options_.maximumScaleFactor);
    }

    [[nodiscard]] double rejectionStepFactor(double error) const {
        if (error <= std::numeric_limits<double>::epsilon()) {
            return options_.minimumScaleFactor;
        }

        const double rawFactor = options_.safetyFactor
            * std::pow(error, -1.0 / 8.0);
        return std::clamp(rawFactor, options_.minimumScaleFactor, 1.0);
    }

    Options options_ {};
};

} // namespace od

#endif // DP853_INTEGRATOR_HPP
