#pragma once

#include <Eigen/Dense>

#include <vector>

namespace od {

class EphemerisInterpolator {
public:
    struct Node {
        double tdb {0.0};
        Eigen::VectorXd state;
        Eigen::VectorXd derivative;
    };

    void addNode(double tdb, const Eigen::VectorXd& state, const Eigen::VectorXd& derivative);

    [[nodiscard]] Eigen::VectorXd interpolate(double tdb) const;

    [[nodiscard]] const std::vector<Node>& history() const noexcept;

private:
    std::vector<Node> history_;
};

} // namespace od
