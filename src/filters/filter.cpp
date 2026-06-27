#include "filters/filter.hpp"

#include <cmath>
#include <stdexcept>

namespace fd::filters {

Filter::~Filter() = default;

void Filter::setInitialState(const Eigen::VectorXd& x0,
                             const Eigen::MatrixXd& P0,
                             double epoch) {
    if (!std::isfinite(epoch)) {
        throw std::invalid_argument("Filter epoch must be finite.");
    }

    validateStateAndCovariance(x0, P0);

    x_hat_ = x0;
    P_ = P0;
    epoch_tk_ = epoch;
    initialized_ = true;
}

const Eigen::VectorXd& Filter::state() const noexcept {
    return x_hat_;
}

const Eigen::MatrixXd& Filter::covariance() const noexcept {
    return P_;
}

double Filter::epoch() const noexcept {
    return epoch_tk_;
}

bool Filter::initialized() const noexcept {
    return initialized_;
}

void Filter::validateStateAndCovariance(const Eigen::VectorXd& state,
                                        const Eigen::MatrixXd& covariance) {
    if (state.size() == 0) {
        throw std::invalid_argument("Filter state cannot be empty.");
    }
    if (covariance.rows() != state.size() || covariance.cols() != state.size()) {
        throw std::invalid_argument("Filter covariance dimensions must match the state dimension.");
    }
    if (!state.allFinite() || !covariance.allFinite()) {
        throw std::invalid_argument("Filter state and covariance must contain only finite values.");
    }
    if (!covariance.isApprox(covariance.transpose(), 1.0e-12)) {
        throw std::invalid_argument("Filter covariance must be symmetric.");
    }
}

} // namespace fd::filters
