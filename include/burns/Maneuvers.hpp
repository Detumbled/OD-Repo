#pragma once

#include <Eigen/Dense>
#include <string>



class Maneuvers {
    public:
        virtual void applyManeuver(const double t, Eigen::VectorXd& state) const = 0;
        virtual ~Maneuvers() = default;

        Maneuvers() = default;
        Maneuvers(const )

    private: 
        



};