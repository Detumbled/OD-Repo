#include "filters/BatchLeastSquaresDriver.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::filters {
namespace {

[[nodiscard]] bool hasConfiguredConvergenceCheck(
    const BatchLeastSquaresDriver::ConvergenceCriteria& criteria) {
    return criteria.correctionNormTolerance > 0.0
        || criteria.weightedRmsTolerance > 0.0
        || criteria.linearizedPostfitWeightedRmsTolerance > 0.0
        || criteria.correctionComponentTolerances.size() > 0
        || !criteria.correctionBlockTolerances.empty();
}

} // namespace

bool BatchLeastSquaresDriver::Result::converged() const noexcept {
    return status == Status::Converged;
}

std::size_t BatchLeastSquaresDriver::Result::iterationCount() const noexcept {
    return iterations.size();
}

const char* BatchLeastSquaresDriver::Result::statusName() const noexcept {
    return BatchLeastSquaresDriver::statusName(status);
}

BatchLeastSquaresDriver::BatchLeastSquaresDriver(WLS& filter,
                                                 LinearizationFunction linearize)
    : BatchLeastSquaresDriver(filter, std::move(linearize), ConvergenceCriteria {}) {}

BatchLeastSquaresDriver::BatchLeastSquaresDriver(WLS& filter,
                                                 LinearizationFunction linearize,
                                                 ConvergenceCriteria criteria)
    : filter_(filter),
      linearize_(std::move(linearize)),
      criteria_(std::move(criteria)) {
    if (!linearize_) {
        throw std::invalid_argument("BatchLeastSquaresDriver requires a linearization function.");
    }
}

void BatchLeastSquaresDriver::setConvergenceCriteria(ConvergenceCriteria criteria) {
    criteria_ = std::move(criteria);
}

void BatchLeastSquaresDriver::setConvergenceFunction(ConvergenceFunction convergenceFunction) {
    convergenceFunction_ = std::move(convergenceFunction);
}

const BatchLeastSquaresDriver::ConvergenceCriteria&
BatchLeastSquaresDriver::convergenceCriteria() const noexcept {
    return criteria_;
}

BatchLeastSquaresDriver::Result BatchLeastSquaresDriver::solve(
    const IterationCallback& iterationCallback) {
    if (!filter_.initialized()) {
        throw std::logic_error("BatchLeastSquaresDriver requires an initialized WLS filter.");
    }
    validateCriteria(criteria_, filter_.state().size());

    Result result;
    result.iterations.reserve(criteria_.maxIterations);

    for (std::size_t iteration = 1; iteration <= criteria_.maxIterations; ++iteration) {
        const Eigen::VectorXd state_before = filter_.state();
        LinearizedProblem problem = linearize_(state_before);
        validateProblem(problem, state_before.size());

        const double weighted_rms = weightedRms(problem.residuals, problem.measurementVariances);
        filter_.processBatchDiagonal(problem.residuals,
                                     problem.designMatrix,
                                     problem.measurementVariances);

        IterationSummary summary;
        summary.iteration = iteration;
        summary.stateBefore = state_before;
        summary.stateAfter = filter_.state();
        summary.residuals = std::move(problem.residuals);
        summary.linearizedPostfitResiduals = filter_.lastPostfitResiduals();
        summary.correction = filter_.lastStateCorrection();
        summary.correctionNorm = summary.correction.norm();
        summary.weightedRms = weighted_rms;
        summary.linearizedPostfitWeightedRms =
            weightedRms(summary.linearizedPostfitResiduals, problem.measurementVariances);
        summary.measurementCount = filter_.accumulatedMeasurementCount();

        if (iterationCallback) {
            iterationCallback(summary);
        }

        const bool converged = convergenceFunction_
            ? convergenceFunction_(summary, criteria_)
            : defaultConverged(summary, criteria_);

        result.iterations.push_back(std::move(summary));
        if (converged) {
            result.status = Status::Converged;
            return result;
        }
    }

    result.status = Status::MaxIterations;
    return result;
}

const char* BatchLeastSquaresDriver::statusName(Status status) noexcept {
    switch (status) {
    case Status::Converged:
        return "CONVERGED";
    case Status::MaxIterations:
        return "MAX_ITERATIONS";
    }
    return "UNKNOWN";
}

double BatchLeastSquaresDriver::weightedRms(const Eigen::VectorXd& residuals,
                                            const Eigen::VectorXd& variances) {
    if (residuals.size() == 0 || residuals.size() != variances.size()) {
        throw std::invalid_argument("Weighted RMS inputs have inconsistent dimensions.");
    }

    double weighted_sum = 0.0;
    for (Eigen::Index row = 0; row < residuals.size(); ++row) {
        const double variance = variances[row];
        if (!(variance > 0.0) || !std::isfinite(variance)) {
            throw std::invalid_argument("Weighted RMS variances must be positive and finite.");
        }
        weighted_sum += residuals[row] * residuals[row] / variance;
    }
    return std::sqrt(weighted_sum / static_cast<double>(residuals.size()));
}

void BatchLeastSquaresDriver::validateCriteria(const ConvergenceCriteria& criteria,
                                               Eigen::Index stateDimension) {
    if (criteria.maxIterations == 0U) {
        throw std::invalid_argument("BatchLeastSquaresDriver max iterations must be positive.");
    }
    if (criteria.minimumIterations == 0U || criteria.minimumIterations > criteria.maxIterations) {
        throw std::invalid_argument("Minimum iterations must be positive and no larger than max iterations.");
    }
    if (!std::isfinite(criteria.correctionNormTolerance)
        || criteria.correctionNormTolerance < 0.0) {
        throw std::invalid_argument("Correction norm tolerance must be finite and non-negative.");
    }
    if (!std::isfinite(criteria.weightedRmsTolerance) || criteria.weightedRmsTolerance < 0.0) {
        throw std::invalid_argument("Weighted RMS tolerance must be finite and non-negative.");
    }
    if (!std::isfinite(criteria.linearizedPostfitWeightedRmsTolerance)
        || criteria.linearizedPostfitWeightedRmsTolerance < 0.0) {
        throw std::invalid_argument("Linearized postfit weighted RMS tolerance must be finite and non-negative.");
    }
    if (criteria.correctionComponentTolerances.size() > 0) {
        if (criteria.correctionComponentTolerances.size() != stateDimension) {
            throw std::invalid_argument("Component correction tolerances must match the state dimension.");
        }
        if (!criteria.correctionComponentTolerances.allFinite()
            || (criteria.correctionComponentTolerances.array() <= 0.0).any()) {
            throw std::invalid_argument("Component correction tolerances must be positive and finite.");
        }
    }

    for (const CorrectionBlockTolerance& block : criteria.correctionBlockTolerances) {
        if (block.offset < 0 || block.size <= 0 || block.offset + block.size > stateDimension) {
            throw std::invalid_argument("Correction block tolerance is outside the state dimension.");
        }
        if (!(block.normTolerance > 0.0) || !std::isfinite(block.normTolerance)) {
            throw std::invalid_argument("Correction block tolerance must be positive and finite.");
        }
    }
}

void BatchLeastSquaresDriver::validateProblem(const LinearizedProblem& problem,
                                              Eigen::Index stateDimension) {
    if (problem.residuals.size() == 0) {
        throw std::invalid_argument("BatchLeastSquaresDriver residual vector cannot be empty.");
    }
    if (problem.designMatrix.rows() != problem.residuals.size()
        || problem.designMatrix.cols() != stateDimension) {
        throw std::invalid_argument("BatchLeastSquaresDriver design matrix dimensions are inconsistent.");
    }
    if (problem.measurementVariances.size() != problem.residuals.size()) {
        throw std::invalid_argument("BatchLeastSquaresDriver variances must match residual dimension.");
    }
    if (!problem.residuals.allFinite()
        || !problem.designMatrix.allFinite()
        || !problem.measurementVariances.allFinite()
        || (problem.measurementVariances.array() <= 0.0).any()) {
        throw std::invalid_argument("BatchLeastSquaresDriver problem data must be finite with positive variances.");
    }
}

bool BatchLeastSquaresDriver::defaultConverged(const IterationSummary& iteration,
                                               const ConvergenceCriteria& criteria) {
    if (!hasConfiguredConvergenceCheck(criteria)) {
        return false;
    }
    if (iteration.iteration < criteria.minimumIterations) {
        return false;
    }
    if (criteria.weightedRmsTolerance > 0.0
        && iteration.weightedRms > criteria.weightedRmsTolerance) {
        return false;
    }
    if (criteria.linearizedPostfitWeightedRmsTolerance > 0.0
        && iteration.linearizedPostfitWeightedRms > criteria.linearizedPostfitWeightedRmsTolerance) {
        return false;
    }
    if (criteria.correctionNormTolerance > 0.0
        && iteration.correctionNorm > criteria.correctionNormTolerance) {
        return false;
    }
    if (criteria.correctionComponentTolerances.size() > 0
        && (iteration.correction.cwiseAbs().array()
            > criteria.correctionComponentTolerances.array()).any()) {
        return false;
    }
    for (const CorrectionBlockTolerance& block : criteria.correctionBlockTolerances) {
        if (iteration.correction.segment(block.offset, block.size).norm() > block.normTolerance) {
            return false;
        }
    }
    return true;
}

} // namespace fd::filters
