#include "perturbations/Gravitational.hpp"

#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::perturbations {
namespace {

using od::throwIfSpiceFailed;

void validateBody(const MassiveBody& body) {
    if (body.name.empty()) {
        throw std::invalid_argument("Third-body perturbing body name cannot be empty.");
    }
    if (!std::isfinite(body.mu) || body.mu <= 0.0) {
        throw std::invalid_argument("Third-body gravitational parameter must be positive and finite.");
    }
}

[[nodiscard]] Eigen::Vector3d bodyPositionRelativeToCentral(const std::string& body,
                                                            double tdb,
                                                            const std::string& frame,
                                                            const std::string& centralBody) {
    od::SpiceErrorModeGuard action_guard;

    SpiceDouble position[3] = {0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkpos_c(body.c_str(),
             tdb,
             frame.c_str(),
             "NONE",
             centralBody.c_str(),
             position,
             &light_time);

    throwIfSpiceFailed("Failed to compute third-body position for " + body);

    return {position[0], position[1], position[2]};
}

} // namespace

ThirdBodyGravity::ThirdBodyGravity(std::string centralBody, std::string frame)
    : centralBody_(std::move(centralBody)),
      frame_(std::move(frame)) {
    if (centralBody_.empty()) {
        throw std::invalid_argument("Third-body central body cannot be empty.");
    }
    if (frame_.empty()) {
        throw std::invalid_argument("Third-body frame cannot be empty.");
    }
}

void ThirdBodyGravity::addBody(const MassiveBody& body) {
    validateBody(body);
    perturbingBodies_.push_back(body);
}

void ThirdBodyGravity::clearBodies() {
    perturbingBodies_.clear();
}

Eigen::Vector3d ThirdBodyGravity::computeAcceleration(double tdb,
                                                      const Eigen::Vector3d& scPosition) const {
    if (!std::isfinite(tdb) || !scPosition.allFinite()) {
        throw std::invalid_argument("Third-body acceleration inputs must be finite.");
    }

    Eigen::Vector3d acceleration = Eigen::Vector3d::Zero();

    for (const MassiveBody& body : perturbingBodies_) {
        validateBody(body);

        const Eigen::Vector3d central_to_body =
            bodyPositionRelativeToCentral(body.name, tdb, frame_, centralBody_);
        const Eigen::Vector3d spacecraft_to_body = central_to_body - scPosition;

        const double body_radius = central_to_body.norm();
        const double sc_body_radius = spacecraft_to_body.norm();
        if (body_radius <= 0.0 || sc_body_radius <= 0.0) {
            throw std::runtime_error("Third-body geometry has a zero separation.");
        }

        acceleration.noalias() += body.mu
            * (spacecraft_to_body / (sc_body_radius * sc_body_radius * sc_body_radius)
               - central_to_body / (body_radius * body_radius * body_radius));
    }

    return acceleration;
}

} // namespace fd::perturbations
