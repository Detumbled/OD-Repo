#include "RKF45Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

void constantDynamics(double,
                      od::RKF45Integrator::ConstStateRef,
                      od::RKF45Integrator::StateRef derivative) {
    derivative[0] = 2.0;
    derivative[1] = -0.5;
}

bool expectNear(const char* label, double actual, double expected, double tolerance) {
    const double error = std::abs(actual - expected);
    if (error <= tolerance) {
        return true;
    }

    std::cerr << "FAIL: " << label << " expected " << expected
              << " but got " << actual << " (error " << error << ")\n";
    return false;
}

} // namespace

int main() {
    try {
        od::RKF45Integrator::State initialState(2);
        initialState << 1.0, 4.0;

        od::RKF45Integrator::Options options;
        options.absoluteTolerance = 1.0e-12;
        options.relativeTolerance = 1.0e-12;
        options.initialStep = 0.5;
        options.maximumStep = 1.0;

        const od::RKF45Integrator integrator(options);
        const auto result = integrator.integrate(0.0, initialState, 10.0, constantDynamics);

        bool passed = true;
        passed &= result.history.size() == result.acceptedSteps + 1;
        if (!passed) {
            std::cerr << "FAIL: history node count does not match accepted steps.\n";
        }
        passed &= expectNear("initial history time", result.history.front().tdb, 0.0, 0.0);
        passed &= expectNear("final history time", result.history.back().tdb, 10.0, 1.0e-12);
        passed &= expectNear("final state x", result.state[0], 21.0, 1.0e-12);
        passed &= expectNear("final state y", result.state[1], -1.0, 1.0e-12);

        od::EphemerisInterpolator ephemeris;
        for (const auto& node : result.history) {
            ephemeris.addNode(node.tdb, node.state, node.derivative);
        }
        const Eigen::VectorXd interpolated = ephemeris.interpolate(3.25);
        passed &= expectNear("interpolated x", interpolated[0], 7.5, 1.0e-12);
        passed &= expectNear("interpolated y", interpolated[1], 2.375, 1.0e-12);

        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "RKF45 history test failed with exception: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
