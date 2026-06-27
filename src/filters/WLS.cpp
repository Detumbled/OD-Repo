#include "filters/WLS.hpp"

#include "StationCatalog.hpp"

#include <SpiceUsr.h>

#include <array>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace fd::filters {
namespace {

using SpiceErrorActionBuffer = std::array<SpiceChar, 32>;

class SpiceErrorActionGuard {
public:
    SpiceErrorActionGuard() {
        erract_c("GET", static_cast<SpiceInt>(previous_action_.size()), previous_action_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
    }

    SpiceErrorActionGuard(const SpiceErrorActionGuard&) = delete;
    SpiceErrorActionGuard& operator=(const SpiceErrorActionGuard&) = delete;

    ~SpiceErrorActionGuard() {
        erract_c("SET", 0, previous_action_.data());
    }

private:
    SpiceErrorActionBuffer previous_action_ {};
};

[[nodiscard]] std::string spiceTargetName(const std::string& stationName) {
    try {
        return std::to_string(od::stationNaifIdFromName(stationName));
    } catch (const std::invalid_argument&) {
        return stationName;
    }
}

void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar short_message[1841] = {0};
    SpiceChar long_message[1841] = {0};
    getmsg_c("SHORT", sizeof(short_message), short_message);
    getmsg_c("LONG", sizeof(long_message), long_message);
    reset_c();

    std::ostringstream message;
    message << context << ": " << short_message << " | " << long_message;
    throw std::runtime_error(message.str());
}

void requireFiniteMatrix(const Eigen::MatrixXd& matrix, const char* label) {
    if (!matrix.allFinite()) {
        throw std::invalid_argument(std::string(label) + " must contain only finite values.");
    }
}

} // namespace

WLS::WLS() = default;

WLS::WLS(std::string integrationCenter)
    : integration_center_(std::move(integrationCenter)) {
    if (integration_center_.empty()) {
        throw std::invalid_argument("WLS integration center cannot be empty.");
    }
}

void WLS::setInitialState(const Eigen::VectorXd& x0,
                          const Eigen::MatrixXd& P0,
                          double epoch) {
    Filter::setInitialState(x0, P0, epoch);

    prior_state_ = x0;
    prior_covariance_ = P0;
    prior_information_ = invertSymmetricPositiveDefinite(P0, "WLS prior covariance");

    last_state_correction_.setZero(x0.size());
    last_postfit_residuals_.resize(0);
    last_information_matrix_.setZero(x0.size(), x0.size());
}

void WLS::processBatch(const Eigen::VectorXd& residuals,
                       const Eigen::MatrixXd& designMatrix,
                       const Eigen::MatrixXd& measurementCovariance) {
    if (!initialized_) {
        throw std::logic_error("WLS filter must be initialized before processing a batch.");
    }

    validateBatchInputs(residuals, designMatrix, measurementCovariance, x_hat_.size());

    // Apply the measurement weighting as solves instead of explicitly forming R^-1.
    const Eigen::MatrixXd weighted_design = solveSymmetricPositiveDefinite(
        measurementCovariance,
        designMatrix,
        "WLS measurement covariance");
    const Eigen::VectorXd weighted_residuals = solveSymmetricPositiveDefinite(
        measurementCovariance,
        residuals,
        "WLS measurement covariance");

    // Classical a-priori batch normal equations:
    // Lambda = H^T R^-1 H + P0^-1
    // N      = H^T R^-1 r + P0^-1 (x_prior - x_nominal)
    const Eigen::VectorXd prior_deviation = prior_state_ - x_hat_;
    last_information_matrix_.noalias() = designMatrix.transpose() * weighted_design;
    last_information_matrix_ += prior_information_;

    Eigen::VectorXd normal_rhs(x_hat_.size());
    normal_rhs.noalias() = designMatrix.transpose() * weighted_residuals;
    normal_rhs += prior_information_ * prior_deviation;

    last_state_correction_ = solveSymmetricPositiveDefinite(
        last_information_matrix_,
        normal_rhs,
        "WLS normal information matrix");

    x_hat_ += last_state_correction_;
    P_ = invertSymmetricPositiveDefinite(last_information_matrix_, "WLS normal information matrix");
    last_postfit_residuals_ = residuals - designMatrix * last_state_correction_;
}

WLS::StationState WLS::getStationState(const std::string& stationName,
                                       double tdb) const {
    if (stationName.empty()) {
        throw std::invalid_argument("SPICE station name cannot be empty.");
    }
    if (!std::isfinite(tdb)) {
        throw std::invalid_argument("SPICE station epoch must be finite.");
    }

    // Kernel loading remains the caller's responsibility; this method only queries SPICE.
    SpiceErrorActionGuard action_guard;

    const std::string target = spiceTargetName(stationName);
    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;

    spkezr_c(target.c_str(),
             tdb,
             spice_frame_.c_str(),
             spice_aberration_correction_.c_str(),
             integration_center_.c_str(),
             spice_state,
             &light_time);

    throwIfSpiceFailed("Failed to fetch station state for " + stationName);

    StationState state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spice_state[i];
    }

    return state;
}

void WLS::setSpiceFrame(std::string frame) {
    if (frame.empty()) {
        throw std::invalid_argument("SPICE frame cannot be empty.");
    }

    spice_frame_ = std::move(frame);
}

void WLS::setSpiceAberrationCorrection(std::string aberrationCorrection) {
    if (aberrationCorrection.empty()) {
        throw std::invalid_argument("SPICE aberration correction cannot be empty.");
    }

    spice_aberration_correction_ = std::move(aberrationCorrection);
}

void WLS::setIntegrationCenter(std::string integrationCenter) {
    if (integrationCenter.empty()) {
        throw std::invalid_argument("WLS integration center cannot be empty.");
    }

    integration_center_ = std::move(integrationCenter);
}

const Eigen::VectorXd& WLS::priorState() const noexcept {
    return prior_state_;
}

const Eigen::MatrixXd& WLS::priorCovariance() const noexcept {
    return prior_covariance_;
}

const Eigen::VectorXd& WLS::lastStateCorrection() const noexcept {
    return last_state_correction_;
}

const Eigen::VectorXd& WLS::lastPostfitResiduals() const noexcept {
    return last_postfit_residuals_;
}

const Eigen::MatrixXd& WLS::lastInformationMatrix() const noexcept {
    return last_information_matrix_;
}

void WLS::validateBatchInputs(const Eigen::VectorXd& residuals,
                              const Eigen::MatrixXd& designMatrix,
                              const Eigen::MatrixXd& measurementCovariance,
                              int stateDimension) {
    if (residuals.size() == 0) {
        throw std::invalid_argument("WLS residual vector cannot be empty.");
    }
    if (designMatrix.rows() != residuals.size() || designMatrix.cols() != stateDimension) {
        throw std::invalid_argument("WLS design matrix dimensions are inconsistent with residual/state dimensions.");
    }
    if (measurementCovariance.rows() != residuals.size()
        || measurementCovariance.cols() != residuals.size()) {
        throw std::invalid_argument("WLS measurement covariance must be square with measurement dimension.");
    }
    if (!residuals.allFinite()) {
        throw std::invalid_argument("WLS residuals must contain only finite values.");
    }

    requireFiniteMatrix(designMatrix, "WLS design matrix");
    requireFiniteMatrix(measurementCovariance, "WLS measurement covariance");

    if (!measurementCovariance.isApprox(measurementCovariance.transpose(), 1.0e-12)) {
        throw std::invalid_argument("WLS measurement covariance must be symmetric.");
    }
}

Eigen::MatrixXd WLS::solveSymmetricPositiveDefinite(const Eigen::MatrixXd& matrix,
                                                    const Eigen::MatrixXd& rhs,
                                                    const char* context) {
    if (matrix.rows() != matrix.cols() || matrix.rows() != rhs.rows()) {
        throw std::invalid_argument(std::string(context) + " solve has incompatible dimensions.");
    }

    const Eigen::LDLT<Eigen::MatrixXd> ldlt(matrix);
    if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
        const Eigen::MatrixXd solution = ldlt.solve(rhs);
        if (solution.allFinite()) {
            return solution;
        }
    }

    const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(matrix);
    if (qr.rank() < matrix.cols()) {
        throw std::runtime_error(std::string(context) + " is rank deficient.");
    }

    const Eigen::MatrixXd solution = qr.solve(rhs);
    if (!solution.allFinite()) {
        throw std::runtime_error(std::string(context) + " solve produced non-finite values.");
    }

    return solution;
}

Eigen::VectorXd WLS::solveSymmetricPositiveDefinite(const Eigen::MatrixXd& matrix,
                                                    const Eigen::VectorXd& rhs,
                                                    const char* context) {
    Eigen::MatrixXd rhs_matrix(rhs.size(), 1);
    rhs_matrix.col(0) = rhs;

    const Eigen::MatrixXd solution = solveSymmetricPositiveDefinite(matrix, rhs_matrix, context);
    return solution.col(0);
}

Eigen::MatrixXd WLS::invertSymmetricPositiveDefinite(const Eigen::MatrixXd& matrix,
                                                     const char* context) {
    if (matrix.rows() != matrix.cols()) {
        throw std::invalid_argument(std::string(context) + " must be square.");
    }

    const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols());
    Eigen::MatrixXd inverse = solveSymmetricPositiveDefinite(matrix, identity, context);
    const Eigen::MatrixXd inverse_transpose = inverse.transpose();
    inverse += inverse_transpose;
    inverse *= 0.5;
    return inverse;
}

} // namespace fd::filters
