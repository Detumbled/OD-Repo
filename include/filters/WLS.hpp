#pragma once

#include "filters/filter.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace fd::filters {

class WLS final : public Filter {
public:
    using StationState = Eigen::Matrix<double, 6, 1>;
    using MeasurementRow = Eigen::RowVectorXd;

    WLS();
    explicit WLS(std::string integrationCenter);

    void setInitialState(const Eigen::VectorXd& x0,
                         const Eigen::MatrixXd& P0,
                         double epoch) override;

    void processBatch(const Eigen::VectorXd& residuals,
                      const Eigen::MatrixXd& designMatrix,
                      const Eigen::MatrixXd& measurementCovariance) override;

    void processBatchDiagonal(const Eigen::VectorXd& residuals,
                              const Eigen::MatrixXd& designMatrix,
                              const Eigen::VectorXd& measurementVariances);

    void setApriori(const Eigen::MatrixXd& covariance);
    void resetAccumulator();

    template <typename Derived>
    void addMeasurement(double residual,
                        const Eigen::MatrixBase<Derived>& jacobianRow,
                        double sigma) {
        if (!(sigma > 0.0) || !std::isfinite(sigma)) {
            throw std::invalid_argument("WLS measurement sigma must be positive and finite.");
        }

        addMeasurementVariance(residual, jacobianRow, sigma * sigma);
    }

    template <typename Derived>
    void addMeasurementVariance(double residual,
                                const Eigen::MatrixBase<Derived>& jacobianRow,
                                double variance) {
        if (!initialized_) {
            throw std::logic_error("WLS filter must be initialized before adding measurements.");
        }
        if (jacobianRow.rows() != 1 || jacobianRow.cols() != x_hat_.size()) {
            throw std::invalid_argument("WLS measurement Jacobian row dimensions are inconsistent with the state.");
        }
        if (!std::isfinite(residual) || !jacobianRow.allFinite()) {
            throw std::invalid_argument("WLS measurement residual and Jacobian row must be finite.");
        }
        if (!(variance > 0.0) || !std::isfinite(variance)) {
            throw std::invalid_argument("WLS measurement variance must be positive and finite.");
        }

        const double weight = 1.0 / variance;
        information_matrix_.noalias() += weight * jacobianRow.transpose() * jacobianRow;
        information_vector_.noalias() += (weight * residual) * jacobianRow.transpose();
        last_information_matrix_ = information_matrix_;
        ++accumulated_measurements_;
    }

    [[nodiscard]] Eigen::VectorXd solve() const;
    [[nodiscard]] Eigen::VectorXd solveAndUpdate();
    [[nodiscard]] std::size_t accumulatedMeasurementCount() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& informationMatrix() const noexcept;
    [[nodiscard]] const Eigen::VectorXd& informationVector() const noexcept;

    [[nodiscard]] StationState getStationState(const std::string& stationName,
                                               double tdb) const;

    void setSpiceFrame(std::string frame);
    void setSpiceAberrationCorrection(std::string aberrationCorrection);
    void setIntegrationCenter(std::string integrationCenter);

    [[nodiscard]] const Eigen::VectorXd& priorState() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& priorCovariance() const noexcept;
    [[nodiscard]] const Eigen::VectorXd& lastStateCorrection() const noexcept;
    [[nodiscard]] const Eigen::VectorXd& lastPostfitResiduals() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& lastInformationMatrix() const noexcept;

private:
    static void validateBatchInputs(const Eigen::VectorXd& residuals,
                                    const Eigen::MatrixXd& designMatrix,
                                    const Eigen::MatrixXd& measurementCovariance,
                                    int stateDimension);

    [[nodiscard]] static Eigen::MatrixXd solveSymmetricPositiveDefinite(
        const Eigen::MatrixXd& matrix,
        const Eigen::MatrixXd& rhs,
        const char* context);

    [[nodiscard]] static Eigen::VectorXd solveSymmetricPositiveDefinite(
        const Eigen::MatrixXd& matrix,
        const Eigen::VectorXd& rhs,
        const char* context);

    [[nodiscard]] static Eigen::MatrixXd invertSymmetricPositiveDefinite(
        const Eigen::MatrixXd& matrix,
        const char* context);

    void applyWeightedNormalEquations(const Eigen::VectorXd& residuals,
                                      const Eigen::MatrixXd& designMatrix,
                                      const Eigen::MatrixXd& weightedDesign,
                                      const Eigen::VectorXd& weightedResiduals);

    Eigen::VectorXd prior_state_;
    Eigen::MatrixXd prior_covariance_;
    Eigen::MatrixXd prior_information_;

    Eigen::MatrixXd information_matrix_;
    Eigen::VectorXd information_vector_;
    std::size_t accumulated_measurements_ {0};

    Eigen::VectorXd last_state_correction_;
    Eigen::VectorXd last_postfit_residuals_;
    Eigen::MatrixXd last_information_matrix_;

    std::string integration_center_ {"EARTH"};
    std::string spice_frame_ {"J2000"};
    std::string spice_aberration_correction_ {"LT+S"};
};

} // namespace fd::filters
