#ifndef SPICE_INTERPOLATOR_HPP
#define SPICE_INTERPOLATOR_HPP

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <cstddef>
#include <string>
#include <vector>

namespace od {

class SpiceInterpolator {
public:
    using State6 = Eigen::Matrix<double, 6, 1>;
    using Vector3 = Eigen::Vector3d;
    using Vector3List = std::vector<Vector3, Eigen::aligned_allocator<Vector3>>;

    void initialize(const std::string& target,
                    const std::string& center,
                    const std::string& frame,
                    double startEt,
                    double endEt,
                    double stepSec,
                    const std::string& aberrationCorrection = "NONE");

    [[nodiscard]] Vector3 getPosition(double et) const;
    [[nodiscard]] Vector3 getVelocity(double et) const;
    [[nodiscard]] State6 getState(double et) const;

    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] double startEt() const noexcept;
    [[nodiscard]] double endEt() const noexcept;
    [[nodiscard]] double stepSec() const noexcept;

private:
    [[nodiscard]] std::size_t segmentIndex(double et) const;
    [[nodiscard]] State6 interpolateState(std::size_t segment, double et) const;

    std::string target_;
    std::string center_;
    std::string frame_;
    std::string aberrationCorrection_ {"NONE"};
    double startEt_ {0.0};
    double endEt_ {0.0};
    double stepSec_ {0.0};
    std::vector<double> times_;
    Vector3List positions_;
    Vector3List velocities_;
};

} // namespace od

#endif // SPICE_INTERPOLATOR_HPP
