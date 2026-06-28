#include "perturbations/Gravitational.hpp"
#include "perturbations/SRP.hpp"

#include <SpiceUsr.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

constexpr const char* kMetaKernel = "../kernels.tm";
constexpr const char* kStartUtc = "1979-03-05T00:00:00";
constexpr const char* kTarget = "-31";
constexpr const char* kCentralBody = "SUN";
constexpr const char* kFrame = "J2000";
constexpr const char* kJupiterBarycenter = "5";
constexpr double kJupiterSystemMuKm3PerSec2 = 126712764.8;
constexpr double kReflectivity = 1.3;
constexpr double kAreaM2 = 10.0;
constexpr double kMassKg = 721.9;
constexpr double kRelativeTolerance = 1.0e-12;

void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar short_message[1841] = {0};
    SpiceChar long_message[1841] = {0};
    getmsg_c("SHORT", sizeof(short_message), short_message);
    getmsg_c("LONG", sizeof(long_message), long_message);
    reset_c();

    throw std::runtime_error(context + ": " + short_message + " | " + long_message);
}

[[nodiscard]] Eigen::Matrix<double, 6, 1> stateRelativeTo(const std::string& target,
                                                          double tdb,
                                                          const std::string& observer) {
    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkezr_c(target.c_str(), tdb, kFrame, "NONE", observer.c_str(), spice_state, &light_time);
    throwIfSpiceFailed("Failed to get SPICE state for " + target);

    Eigen::Matrix<double, 6, 1> state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spice_state[i];
    }
    return state;
}

[[nodiscard]] Eigen::Vector3d positionRelativeTo(const std::string& target,
                                                 double tdb,
                                                 const std::string& observer) {
    SpiceDouble spice_position[3] = {0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkpos_c(target.c_str(), tdb, kFrame, "NONE", observer.c_str(), spice_position, &light_time);
    throwIfSpiceFailed("Failed to get SPICE position for " + target);
    return {spice_position[0], spice_position[1], spice_position[2]};
}

[[nodiscard]] Eigen::Vector3d referenceThirdBodyAcceleration(const Eigen::Vector3d& sc_position_km,
                                                             const Eigen::Vector3d& body_position_km,
                                                             double body_mu) {
    const Eigen::Vector3d spacecraft_to_body = body_position_km - sc_position_km;
    const double sc_body_radius = spacecraft_to_body.norm();
    const double body_radius = body_position_km.norm();

    return body_mu
        * (spacecraft_to_body / (sc_body_radius * sc_body_radius * sc_body_radius)
           - body_position_km / (body_radius * body_radius * body_radius));
}

[[nodiscard]] Eigen::Vector3d referenceSrpAcceleration(const Eigen::Vector3d& sc_position_km) {
    const double distance_km = sc_position_km.norm();
    const Eigen::Vector3d anti_sunward_unit = sc_position_km / distance_km;
    const double acceleration_m_per_s2 =
        fd::perturbations::SolarRadiationPressure::kSolarPressure1AU_N_m2
        * kReflectivity
        * (kAreaM2 / kMassKg)
        * std::pow(fd::perturbations::SolarRadiationPressure::kAu_km / distance_km, 2.0);

    return (acceleration_m_per_s2 / 1000.0) * anti_sunward_unit;
}

[[nodiscard]] bool relativeClose(const Eigen::Vector3d& actual,
                                 const Eigen::Vector3d& expected,
                                 double tolerance) {
    const double scale = std::max(1.0, expected.norm());
    return (actual - expected).norm() / scale < tolerance;
}

void printVector(const std::string& label, const Eigen::Vector3d& value) {
    std::cout << "  " << label << " = [" << value.transpose()
              << "] km/s^2, norm = " << value.norm() << '\n';
}

} // namespace

int main() {
    try {
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
        errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));

        furnsh_c(kMetaKernel);
        throwIfSpiceFailed("Failed to load perturbation test meta-kernel");

        SpiceDouble start_et = 0.0;
        str2et_c(kStartUtc, &start_et);
        throwIfSpiceFailed("Failed to convert perturbation test start epoch");

        const Eigen::Matrix<double, 6, 1> voyager_state = stateRelativeTo(kTarget, start_et, kCentralBody);
        const Eigen::Vector3d sc_position = voyager_state.segment<3>(0);
        const Eigen::Vector3d jupiter_position = positionRelativeTo(kJupiterBarycenter, start_et, kCentralBody);

        fd::perturbations::ThirdBodyGravity third_body(kCentralBody, kFrame);
        third_body.addBody({kJupiterBarycenter, kJupiterSystemMuKm3PerSec2});
        const Eigen::Vector3d third_body_acceleration =
            third_body.computeAcceleration(start_et, sc_position);
        const Eigen::Vector3d third_body_reference =
            referenceThirdBodyAcceleration(sc_position, jupiter_position, kJupiterSystemMuKm3PerSec2);

        fd::perturbations::SolarRadiationPressure srp(kCentralBody,
                                                      kFrame,
                                                      kReflectivity,
                                                      kAreaM2,
                                                      kMassKg);
        const Eigen::Vector3d srp_acceleration = srp.computeAcceleration(start_et, sc_position);
        const Eigen::Vector3d srp_reference = referenceSrpAcceleration(sc_position);

        std::cout << std::scientific << std::setprecision(12)
                  << "Perturbation verification at " << kStartUtc << '\n'
                  << "  Voyager Sun-centered position norm = " << sc_position.norm() << " km\n";
        printVector("Jupiter third-body acceleration", third_body_acceleration);
        printVector("SRP acceleration", srp_acceleration);
        std::cout << "  third-body relative error = "
                  << (third_body_acceleration - third_body_reference).norm()
                     / std::max(1.0, third_body_reference.norm()) << '\n'
                  << "  SRP relative error        = "
                  << (srp_acceleration - srp_reference).norm()
                     / std::max(1.0, srp_reference.norm()) << '\n';

        if (!third_body_acceleration.allFinite() || third_body_acceleration.norm() <= 0.0) {
            std::cerr << "FAIL: third-body acceleration is invalid.\n";
            return EXIT_FAILURE;
        }
        if (!srp_acceleration.allFinite() || srp_acceleration.norm() <= 0.0) {
            std::cerr << "FAIL: SRP acceleration is invalid.\n";
            return EXIT_FAILURE;
        }
        if (!relativeClose(third_body_acceleration, third_body_reference, kRelativeTolerance)) {
            std::cerr << "FAIL: third-body acceleration does not match reference computation.\n";
            return EXIT_FAILURE;
        }
        if (!relativeClose(srp_acceleration, srp_reference, kRelativeTolerance)) {
            std::cerr << "FAIL: SRP acceleration does not match reference computation.\n";
            return EXIT_FAILURE;
        }

        kclear_c();
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Perturbation verification failed: " << error.what() << '\n';
        kclear_c();
        return EXIT_FAILURE;
    }
}
