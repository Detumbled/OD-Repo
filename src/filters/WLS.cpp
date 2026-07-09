#include "filters/WLS.hpp"

#include "stations/StationCatalog.hpp"
#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::filters {
namespace {

using od::throwIfSpiceFailed;

[[nodiscard]] std::string spiceTargetName(const std::string& stationName) {
    try {
        return std::to_string(od::stationNaifIdFromName(stationName));
    } catch (const std::invalid_argument&) {
        return stationName;
    }
}

void requireFiniteMatrix(const Eigen::MatrixXd& matrix, const char* label) {
    if (!matrix.allFinite()) {
        throw std::invalid_argument(std::string(label) + " must contain only finite values.");
    }
}

[[nodiscard]] bool isEffectivelyDiagonal(const Eigen::MatrixXd& matrix) {
    constexpr double kAbsoluteTolerance = 1.0e-30;

    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index column = 0; column < matrix.cols(); ++column) {
            if (row != column && std::abs(matrix(row, column)) > kAbsoluteTolerance) {
                return false;
            }
        }
    }

    return true;
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
    setApriori(P0);
}

void WLS::processBatch(const Eigen::VectorXd& residuals,
                       const Eigen::MatrixXd& designMatrix,
                       const Eigen::MatrixXd& measurementCovariance) {
    if (!initialized_) {
        throw std::logic_error("WLS filter must be initialized before processing a batch.");
    }

    validateBatchInputs(residuals, designMatrix, measurementCovariance, x_hat_.size());

    if (isEffectivelyDiagonal(measurementCovariance)) {
        processBatchDiagonal(residuals, designMatrix, measurementCovariance.diagonal());
        return;
    }

    Eigen::MatrixXd weighted_design;
    Eigen::VectorXd weighted_residuals;
    // Apply the measurement weighting as solves instead of explicitly forming R^-1.
    weighted_design = solveSymmetricPositiveDefinite(
        measurementCovariance,
        designMatrix,
        "WLS measurement covariance");
    weighted_residuals = solveSymmetricPositiveDefinite(
        measurementCovariance,
        residuals,
        "WLS measurement covariance");

    applyWeightedNormalEquations(residuals, designMatrix, weighted_design, weighted_residuals);
}

void WLS::setApriori(const Eigen::MatrixXd& covariance) {
    if (covariance.rows() == 0 || covariance.rows() != covariance.cols()) {
        throw std::invalid_argument("WLS a-priori covariance must be non-empty and square.");
    }
    requireFiniteMatrix(covariance, "WLS a-priori covariance");
    if (!covariance.isApprox(covariance.transpose(), 1.0e-12)) {
        throw std::invalid_argument("WLS a-priori covariance must be symmetric.");
    }

    if (!initialized_) {
        x_hat_.setZero(covariance.rows());
        P_ = covariance;
        prior_state_ = x_hat_;
        epoch_tk_ = 0.0;
        initialized_ = true;
    } else {
        if (covariance.rows() != x_hat_.size()) {
            throw std::invalid_argument("WLS a-priori covariance dimensions must match the state dimension.");
        }
        P_ = covariance;
        if (prior_state_.size() != x_hat_.size()) {
            prior_state_ = x_hat_;
        }
    }

    prior_covariance_ = covariance;
    prior_information_ = invertSymmetricPositiveDefinite(covariance, "WLS prior covariance");

    last_state_correction_.setZero(x_hat_.size());
    last_postfit_residuals_.resize(0);
    resetAccumulator();
}

void WLS::resetAccumulator() {
    if (!initialized_) {
        throw std::logic_error("WLS filter must be initialized before resetting the accumulator.");
    }
    if (prior_information_.rows() != x_hat_.size() || prior_information_.cols() != x_hat_.size()) {
        throw std::logic_error("WLS a-priori information matrix is not initialized.");
    }

    information_matrix_ = prior_information_;
    information_vector_ = prior_information_ * (prior_state_ - x_hat_);
    accumulated_measurements_ = 0;
    last_information_matrix_ = information_matrix_;
}

Eigen::VectorXd WLS::solve() const {
    if (!initialized_) {
        throw std::logic_error("WLS filter must be initialized before solving.");
    }
    if (information_matrix_.rows() != x_hat_.size()
        || information_matrix_.cols() != x_hat_.size()
        || information_vector_.size() != x_hat_.size()) {
        throw std::logic_error("WLS information accumulator is not initialized.");
    }

    return solveSymmetricPositiveDefinite(
        information_matrix_,
        information_vector_,
        "WLS normal information matrix");
}

Eigen::VectorXd WLS::solveAndUpdate() {
    const Eigen::VectorXd correction = solve();

    last_state_correction_ = correction;
    last_information_matrix_ = information_matrix_;
    x_hat_ += correction;
    P_ = invertSymmetricPositiveDefinite(information_matrix_, "WLS normal information matrix");
    last_postfit_residuals_.resize(0);
    return correction;
}

std::size_t WLS::accumulatedMeasurementCount() const noexcept {
    return accumulated_measurements_;
}

const Eigen::MatrixXd& WLS::informationMatrix() const noexcept {
    return information_matrix_;
}

const Eigen::VectorXd& WLS::informationVector() const noexcept {
    return information_vector_;
}

void WLS::processBatchDiagonal(const Eigen::VectorXd& residuals,
                               const Eigen::MatrixXd& designMatrix,
                               const Eigen::VectorXd& measurementVariances) {
    if (!initialized_) {
        throw std::logic_error("WLS filter must be initialized before processing a batch.");
    }
    if (residuals.size() == 0) {
        throw std::invalid_argument("WLS residual vector cannot be empty.");
    }
    if (designMatrix.rows() != residuals.size() || designMatrix.cols() != x_hat_.size()) {
        throw std::invalid_argument("WLS design matrix dimensions are inconsistent with residual/state dimensions.");
    }
    if (measurementVariances.size() != residuals.size()) {
        throw std::invalid_argument("WLS measurement variance vector must match the residual dimension.");
    }
    if (!residuals.allFinite() || !measurementVariances.allFinite()) {
        throw std::invalid_argument("WLS residuals and measurement variances must be finite.");
    }

    requireFiniteMatrix(designMatrix, "WLS design matrix");

    for (Eigen::Index row = 0; row < measurementVariances.size(); ++row) {
        const double variance = measurementVariances[row];
        if (!(variance > 0.0) || !std::isfinite(variance)) {
            throw std::invalid_argument("WLS measurement variances must be positive and finite.");
        }
    }

    resetAccumulator();
    for (Eigen::Index row = 0; row < residuals.size(); ++row) {
        addMeasurementVariance(residuals[row], designMatrix.row(row), measurementVariances[row]);
    }

    const Eigen::VectorXd correction = solveAndUpdate();
    last_postfit_residuals_ = residuals - designMatrix * correction;
}

void WLS::applyWeightedNormalEquations(const Eigen::VectorXd& residuals,
                                       const Eigen::MatrixXd& designMatrix,
                                       const Eigen::MatrixXd& weightedDesign,
                                       const Eigen::VectorXd& weightedResiduals) {
    // Classical a-priori batch normal equations:
    // Lambda = H^T R^-1 H + P0^-1
    // N      = H^T R^-1 r + P0^-1 (x_prior - x_nominal)
    const Eigen::VectorXd prior_deviation = prior_state_ - x_hat_;
    information_matrix_.noalias() = designMatrix.transpose() * weightedDesign;
    information_matrix_ += prior_information_;

    Eigen::VectorXd normal_rhs(x_hat_.size());
    normal_rhs.noalias() = designMatrix.transpose() * weightedResiduals;
    normal_rhs += prior_information_ * prior_deviation;
    information_vector_ = normal_rhs;
    accumulated_measurements_ = static_cast<std::size_t>(residuals.size());
    last_information_matrix_ = information_matrix_;

    last_state_correction_ = solveAndUpdate();

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
    od::SpiceErrorModeGuard action_guard;

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
