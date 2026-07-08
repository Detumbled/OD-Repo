#pragma once

#include "filters/WLS.hpp"

#include <Eigen/Dense>

#include <cstddef>
#include <functional>
#include <vector>

namespace fd::filters {

class BatchLeastSquaresDriver final {
public:
    struct LinearizedProblem {
        Eigen::VectorXd residuals;
        Eigen::MatrixXd designMatrix;
        Eigen::VectorXd measurementVariances;
    };

    struct CorrectionBlockTolerance {
        Eigen::Index offset {0};
        Eigen::Index size {0};
        double normTolerance {0.0};
    };

    struct ConvergenceCriteria {
        std::size_t maxIterations {10};
        std::size_t minimumIterations {1};
        double correctionNormTolerance {0.0};
        double weightedRmsTolerance {0.0};
        double linearizedPostfitWeightedRmsTolerance {0.0};
        Eigen::VectorXd correctionComponentTolerances;
        std::vector<CorrectionBlockTolerance> correctionBlockTolerances;
    };

    struct IterationSummary {
        std::size_t iteration {0};
        Eigen::VectorXd stateBefore;
        Eigen::VectorXd stateAfter;
        Eigen::VectorXd residuals;
        Eigen::VectorXd linearizedPostfitResiduals;
        Eigen::VectorXd correction;
        double correctionNorm {0.0};
        double weightedRms {0.0};
        double linearizedPostfitWeightedRms {0.0};
        std::size_t measurementCount {0};
    };

    enum class Status {
        Converged,
        MaxIterations
    };

    struct Result {
        Status status {Status::MaxIterations};
        std::vector<IterationSummary> iterations;

        [[nodiscard]] bool converged() const noexcept;
        [[nodiscard]] std::size_t iterationCount() const noexcept;
        [[nodiscard]] const char* statusName() const noexcept;
    };

    using LinearizationFunction = std::function<LinearizedProblem(const Eigen::VectorXd& state)>;
    using ConvergenceFunction = std::function<bool(const IterationSummary& iteration,
                                                   const ConvergenceCriteria& criteria)>;
    using IterationCallback = std::function<void(const IterationSummary& iteration)>;

    BatchLeastSquaresDriver(WLS& filter,
                            LinearizationFunction linearize);

    BatchLeastSquaresDriver(WLS& filter,
                            LinearizationFunction linearize,
                            ConvergenceCriteria criteria);

    void setConvergenceCriteria(ConvergenceCriteria criteria);
    void setConvergenceFunction(ConvergenceFunction convergenceFunction);

    [[nodiscard]] const ConvergenceCriteria& convergenceCriteria() const noexcept;
    [[nodiscard]] Result solve(const IterationCallback& iterationCallback = {});

    [[nodiscard]] static const char* statusName(Status status) noexcept;

private:
    [[nodiscard]] static double weightedRms(const Eigen::VectorXd& residuals,
                                            const Eigen::VectorXd& variances);
    static void validateCriteria(const ConvergenceCriteria& criteria,
                                 Eigen::Index stateDimension);
    static void validateProblem(const LinearizedProblem& problem,
                                Eigen::Index stateDimension);
    [[nodiscard]] static bool defaultConverged(const IterationSummary& iteration,
                                               const ConvergenceCriteria& criteria);

    WLS& filter_;
    LinearizationFunction linearize_;
    ConvergenceCriteria criteria_;
    ConvergenceFunction convergenceFunction_;
};

} // namespace fd::filters
