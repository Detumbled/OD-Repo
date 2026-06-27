#pragma once

#include <Eigen/Dense>

namespace fd::filters {

class Filter {
public:
    virtual ~Filter();

    virtual void processBatch(const Eigen::VectorXd& residuals,
                              const Eigen::MatrixXd& designMatrix,
                              const Eigen::MatrixXd& measurementCovariance) = 0;

    virtual void setInitialState(const Eigen::VectorXd& x0,
                                 const Eigen::MatrixXd& P0,
                                 double epoch);

    [[nodiscard]] const Eigen::VectorXd& state() const noexcept;
    [[nodiscard]] const Eigen::MatrixXd& covariance() const noexcept;
    [[nodiscard]] double epoch() const noexcept;
    [[nodiscard]] bool initialized() const noexcept;

protected:
    static void validateStateAndCovariance(const Eigen::VectorXd& state,
                                           const Eigen::MatrixXd& covariance);

    Eigen::VectorXd x_hat_;
    Eigen::MatrixXd P_;
    double epoch_tk_ {0.0};
    bool initialized_ {false};
};

} // namespace fd::filters
