#include "dynamics/EphemerisInterpolator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace od {
namespace {

void validateNode(double tdb, const Eigen::VectorXd& state, const Eigen::VectorXd& derivative) {
    if (!std::isfinite(tdb)) {
        throw std::invalid_argument("Ephemeris node time must be finite.");
    }
    if (state.size() == 0) {
        throw std::invalid_argument("Ephemeris node state cannot be empty.");
    }
    if (derivative.size() != state.size()) {
        throw std::invalid_argument("Ephemeris node derivative dimension must match state dimension.");
    }
    if (!state.allFinite() || !derivative.allFinite()) {
        throw std::invalid_argument("Ephemeris node state and derivative must be finite.");
    }
}

void validateCompatibleNode(const EphemerisInterpolator::Node& node, Eigen::Index dimension) {
    if (node.state.size() != dimension || node.derivative.size() != dimension) {
        throw std::runtime_error("Ephemeris history contains inconsistent node dimensions.");
    }
}

} // namespace

void EphemerisInterpolator::addNode(double tdb,
                                    const Eigen::VectorXd& state,
                                    const Eigen::VectorXd& derivative) {
    validateNode(tdb, state, derivative);

    if (!history_.empty() && history_.front().state.size() != state.size()) {
        throw std::invalid_argument("Ephemeris node dimension does not match existing history.");
    }

    Node node{tdb, state, derivative};
    const auto insertion = std::lower_bound(history_.begin(),
                                            history_.end(),
                                            tdb,
                                            [](const Node& existing, double value) {
                                                return existing.tdb < value;
                                            });

    if (insertion != history_.end() && insertion->tdb == tdb) {
        *insertion = std::move(node);
        return;
    }

    history_.insert(insertion, std::move(node));
}

Eigen::VectorXd EphemerisInterpolator::interpolate(double tdb) const {
    if (!std::isfinite(tdb)) {
        throw std::invalid_argument("Ephemeris interpolation time must be finite.");
    }
    if (history_.empty()) {
        throw std::runtime_error("Cannot interpolate an empty ephemeris history.");
    }

    const Node& first = history_.front();
    const Node& last = history_.back();
    validateCompatibleNode(first, first.state.size());
    validateCompatibleNode(last, first.state.size());

    if (tdb == first.tdb) {
        return first.state;
    }
    if (tdb == last.tdb) {
        return last.state;
    }
    if (tdb < first.tdb || tdb > last.tdb) {
        throw std::out_of_range("Ephemeris interpolation time is outside the stored history.");
    }
    if (history_.size() < 2) {
        throw std::runtime_error("At least two ephemeris nodes are required for interpolation.");
    }

    const auto upper = std::upper_bound(history_.begin(),
                                        history_.end(),
                                        tdb,
                                        [](double value, const Node& node) {
                                            return value < node.tdb;
                                        });
    if (upper == history_.begin() || upper == history_.end()) {
        throw std::runtime_error("Failed to locate ephemeris interpolation interval.");
    }

    const Node& node1 = *upper;
    const Node& node0 = *(upper - 1);
    validateCompatibleNode(node0, first.state.size());
    validateCompatibleNode(node1, first.state.size());

    const double delta_t = node1.tdb - node0.tdb;
    if (!(delta_t > 0.0) || !std::isfinite(delta_t)) {
        throw std::runtime_error("Ephemeris interpolation interval must have positive finite duration.");
    }

    const double s = (tdb - node0.tdb) / delta_t;
    const double s2 = s * s;
    const double s3 = s2 * s;
    const double h00 = 2.0 * s3 - 3.0 * s2 + 1.0;
    const double h10 = s3 - 2.0 * s2 + s;
    const double h01 = -2.0 * s3 + 3.0 * s2;
    const double h11 = s3 - s2;

    return h00 * node0.state
        + h10 * delta_t * node0.derivative
        + h01 * node1.state
        + h11 * delta_t * node1.derivative;
}

const std::vector<EphemerisInterpolator::Node>& EphemerisInterpolator::history() const noexcept {
    return history_;
}

} // namespace od
