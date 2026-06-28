#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>

namespace fd::perturbations {

struct MassiveBody {
    std::string name;
    double mu {0.0};
};

class ThirdBodyGravity {
public:
    ThirdBodyGravity(std::string centralBody, std::string frame);

    void addBody(const MassiveBody& body);
    void clearBodies();

    [[nodiscard]] Eigen::Vector3d computeAcceleration(double tdb,
                                                      const Eigen::Vector3d& scPosition) const;

private:
    std::string centralBody_;
    std::string frame_;
    std::vector<MassiveBody> perturbingBodies_;
};

} // namespace fd::perturbations
