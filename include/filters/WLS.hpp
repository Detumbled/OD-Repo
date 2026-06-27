#pragma once

#include "filters/filter.hpp"

#include <Eigen/Dense>

#include <string>

namespace fd::filters {

class WLS final : public Filter {
public:
    using StationState = Eigen::Matrix<double, 6, 1>;

    WLS();
    explicit WLS(std::string integrationCenter);

    void setInitialState(const Eigen::VectorXd& x0,
                         const Eigen::MatrixXd& P0,
                         double epoch) override;

    void processBatch(const Eigen::VectorXd& residuals,
                      const Eigen::MatrixXd& designMatrix,
                      const Eigen::MatrixXd& measurementCovariance) override;

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

    Eigen::VectorXd prior_state_;
    Eigen::MatrixXd prior_covariance_;
    Eigen::MatrixXd prior_information_;

    Eigen::VectorXd last_state_correction_;
    Eigen::VectorXd last_postfit_residuals_;
    Eigen::MatrixXd last_information_matrix_;

    std::string integration_center_ {"EARTH"};
    std::string spice_frame_ {"J2000"};
    std::string spice_aberration_correction_ {"LT+S"};
};

} // namespace fd::filters
