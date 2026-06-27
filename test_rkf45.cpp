#include "RKF45Integrator.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

constexpr double kEarthMuKm3PerSec2 = 398600.4418;
constexpr double kPi = 3.141592653589793238462643383279502884;

class Rkf45KeplerTest {
public:
    bool run() const {
        const double radius_km = 7000.0;
        const double speed_km_per_sec = std::sqrt(kEarthMuKm3PerSec2 / radius_km);
        // Propagate for 60 periods:
        const double period_sec = 60.0*2.0 * kPi * std::sqrt((radius_km * radius_km * radius_km) / kEarthMuKm3PerSec2);

        od::RKF45Integrator::State initial_state(6);
        initial_state << radius_km,
                         0.0,
                         0.0,
                         0.0,
                         speed_km_per_sec,
                         0.0;

        od::RKF45Integrator::Options options;
        options.absoluteTolerance = 1.0e-11;
        options.relativeTolerance = 1.0e-11;
        options.initialStep = 30.0;
        options.minimumStep = 1.0e-6;
        options.maximumStep = 120.0;

        const od::RKF45Integrator integrator(options);
        const auto result = integrator.integrate(0.0, initial_state, period_sec, twoBodyDynamics);

        const double position_error_km = normDifference(result.state, initial_state, 0);
        const double velocity_error_km_per_sec = normDifference(result.state, initial_state, 3);
        const double initial_energy = specificEnergy(initial_state);
        const double initial_angular_momentum_z = specificAngularMomentumZ(initial_state);
        const double relative_energy_error = std::abs((specificEnergy(result.state) - initial_energy)
                                                      / initial_energy);
        const double relative_angular_momentum_error = std::abs(
            (specificAngularMomentumZ(result.state) - initial_angular_momentum_z) / initial_angular_momentum_z);

        std::cout << std::scientific << std::setprecision(6)
                  << "RKF45 Kepler one-period test\n"
                  << "  accepted steps          : " << result.acceptedSteps << '\n'
                  << "  rejected steps          : " << result.rejectedSteps << '\n'
                  << "  position closure error  : " << position_error_km << " km\n"
                  << "  velocity closure error  : " << velocity_error_km_per_sec << " km/s\n"
                  << "  relative energy error   : " << relative_energy_error << '\n'
                  << "  relative h_z error      : " << relative_angular_momentum_error << '\n';

        bool passed = true;
        passed &= expectLessThan("position closure", position_error_km, 1.0e-5);
        passed &= expectLessThan("velocity closure", velocity_error_km_per_sec, 1.0e-8);
        passed &= expectLessThan("relative specific energy", relative_energy_error, 5.0e-10);
        passed &= expectLessThan("relative angular momentum", relative_angular_momentum_error, 5.0e-10);

        return passed;
    }

private:
    static void twoBodyDynamics(double,
                                od::RKF45Integrator::ConstStateRef state,
                                od::RKF45Integrator::StateRef state_derivative) {
        if (state.size() != 6) {
            throw std::invalid_argument("Kepler test state must have six elements.");
        }

        const double x = state[0];
        const double y = state[1];
        const double z = state[2];
        const double radius_squared = x * x + y * y + z * z;
        const double radius = std::sqrt(radius_squared);
        const double radius_cubed = radius_squared * radius;

        state_derivative[0] = state[3];
        state_derivative[1] = state[4];
        state_derivative[2] = state[5];
        state_derivative[3] = -kEarthMuKm3PerSec2 * x / radius_cubed;
        state_derivative[4] = -kEarthMuKm3PerSec2 * y / radius_cubed;
        state_derivative[5] = -kEarthMuKm3PerSec2 * z / radius_cubed;
    }

    [[nodiscard]] static double normDifference(const od::RKF45Integrator::State& lhs,
                                               const od::RKF45Integrator::State& rhs,
                                               std::size_t offset) {
        const double dx = lhs[offset] - rhs[offset];
        const double dy = lhs[offset + 1] - rhs[offset + 1];
        const double dz = lhs[offset + 2] - rhs[offset + 2];
        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    [[nodiscard]] static double specificEnergy(const od::RKF45Integrator::State& state) {
        const double radius = std::sqrt(state[0] * state[0] + state[1] * state[1] + state[2] * state[2]);
        const double speed_squared = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];
        return 0.5 * speed_squared - kEarthMuKm3PerSec2 / radius;
    }

    [[nodiscard]] static double specificAngularMomentumZ(const od::RKF45Integrator::State& state) {
        return state[0] * state[4] - state[1] * state[3];
    }

    static bool expectLessThan(const std::string& label, double actual, double limit) {
        if (actual < limit) {
            return true;
        }

        std::cerr << "FAIL: " << label << " error " << actual
                  << " exceeded limit " << limit << '\n';
        return false;
    }
};

} // namespace

int main() {
    try {
        const Rkf45KeplerTest test;
        const auto start_time = std::chrono::steady_clock::now();
        const bool passed = test.run();
        const auto end_time = std::chrono::steady_clock::now();
        const std::chrono::duration<double, std::milli> elapsed_time = end_time - start_time;

        std::cout << "  elapsed time            : " << elapsed_time.count() << " ms\n";
        return passed ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (const std::exception& error) {
        std::cerr << "RKF45 Kepler test failed with exception: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
