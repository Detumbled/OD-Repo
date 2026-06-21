 #include <iostream>
#include <ceres/jet.h>

// 1. Define the 7-State Dynamics (Example 4.2.1)
template <typename T>
void ComputeDynamics(const T state[7], T state_dot[7]) {
    // Extract variables from the augmented state vector
    T x  = state[0];
    T y  = state[1];
    T u  = state[2];
    T v  = state[3];
    T mu = state[4]; // mu is now a state variable, not a constant
    T xs = state[5]; // station X
    T ys = state[6]; // station Y

    // Calculate r = sqrt(X^2 + Y^2) 
    // (Note: Example 4.2.1 is a planar 2D problem)
    T r2 = x*x + y*y;
    T r  = ceres::sqrt(r2);
    T r3 = r2 * r;

    // F1, F2: Velocity dynamics
    state_dot[0] = u;
    state_dot[1] = v;

    // F3, F4: Acceleration dynamics
    state_dot[2] = -mu * x / r3;
    state_dot[3] = -mu * y / r3;

    // F5, F6, F7: Static parameters have zero derivative
    state_dot[4] = T(0.0);
    state_dot[5] = T(0.0);
    state_dot[6] = T(0.0);
}

int main() {
    // Nominal reference state X* [X, Y, U, V, mu, Xs, Ys]
    double X_star[7] = {7000.0, 0.0, 0.0, 7.5, 398600.4415, 6378.0, 0.0};

    // Setup Ceres Jets for a 7-dimensional problem
    ceres::Jet<double, 7> X_jet[7];
    ceres::Jet<double, 7> X_dot_jet[7];

    // Initialize values and seed the Identity matrix for AutoDiff
    for (int i = 0; i < 7; ++i) {
        X_jet[i].a = X_star[i];       
        X_jet[i].v.setZero();         
        X_jet[i].v[i] = 1.0;          
    }

    // Evaluate the augmented dynamics
    ComputeDynamics(X_jet, X_dot_jet);

    // Extract and print the 7x7 A(t) Matrix
    std::cout << "The Corrected 7x7 A(t) Matrix for Example 4.2.1:\n\n";
    
    for (int i = 0; i < 7; ++i) { // Rows: F1 through F7
        for (int j = 0; j < 7; ++j) { // Columns: X1 through X7
            std::cout << X_dot_jet[i].v[j] << "\t";
        }
        std::cout << "\n";
    }

    return 0;
}