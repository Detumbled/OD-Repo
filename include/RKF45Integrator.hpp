#ifndef RKF45_INTEGRATOR_HPP
#define RKF45_INTEGRATOR_HPP

#include <Eigen/Dense>

#include <cstddef>
#include <functional>

namespace od {

class RKF45Integrator {
public:
    using State = Eigen::VectorXd;
    using ConstStateRef = Eigen::Ref<const State>;
    using StateRef = Eigen::Ref<State>;
    using DerivativeFunction = std::function<void(double time, ConstStateRef state, StateRef stateDerivative)>;

    struct Options {
        double absoluteTolerance {1.0e-9};
        double relativeTolerance {1.0e-9};
        double initialStep {60.0};
        double minimumStep {1.0e-9};
        double maximumStep {3600.0};
        double safetyFactor {0.9};
        double minimumScaleFactor {0.2};
        double maximumScaleFactor {5.0};
        std::size_t maximumSteps {1000000};
    };

    struct Result {
        State state;
        double finalTime {0.0};
        std::size_t acceptedSteps {0};
        std::size_t rejectedSteps {0};
    };

    RKF45Integrator();
    explicit RKF45Integrator(Options options);

    [[nodiscard]] const Options& options() const noexcept;

    [[nodiscard]] Result integrate(double initialTime,
                                   const State& initialState,
                                   double finalTime,
                                   const DerivativeFunction& dynamics) const;

private:
    Options options_;
};

} // namespace od

#endif // RKF45_INTEGRATOR_HPP
