#include "dynamics/EphemerisInterpolator.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

namespace {

[[nodiscard]] Eigen::VectorXd cubicState(double t) {
    Eigen::VectorXd state(2);
    state << t * t * t - 2.0 * t * t + t + 5.0,
             -0.5 * t * t * t + t * t + 2.0 * t - 3.0;
    return state;
}

[[nodiscard]] Eigen::VectorXd cubicDerivative(double t) {
    Eigen::VectorXd derivative(2);
    derivative << 3.0 * t * t - 4.0 * t + 1.0,
                  -1.5 * t * t + 2.0 * t + 2.0;
    return derivative;
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

bool expectThrowsOutOfRange(const od::EphemerisInterpolator& ephemeris) {
    try {
        (void) ephemeris.interpolate(-1.0);
    } catch (const std::out_of_range&) {
        return true;
    }

    std::cerr << "FAIL: interpolation outside stored history did not throw std::out_of_range.\n";
    return false;
}

} // namespace

int main() {
    try {
        od::EphemerisInterpolator ephemeris;
        ephemeris.addNode(2.0, cubicState(2.0), cubicDerivative(2.0));
        ephemeris.addNode(0.0, cubicState(0.0), cubicDerivative(0.0));

        const Eigen::VectorXd interpolated = ephemeris.interpolate(0.75);
        const Eigen::VectorXd expected = cubicState(0.75);
        const Eigen::VectorXd endpoint = ephemeris.interpolate(2.0);

        bool passed = true;
        passed &= expectNear("component 0", interpolated[0], expected[0], 1.0e-13);
        passed &= expectNear("component 1", interpolated[1], expected[1], 1.0e-13);
        passed &= expectNear("endpoint component 0", endpoint[0], cubicState(2.0)[0], 0.0);
        passed &= expectNear("endpoint component 1", endpoint[1], cubicState(2.0)[1], 0.0);
        passed &= expectThrowsOutOfRange(ephemeris);

        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "Ephemeris interpolator test failed with exception: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
