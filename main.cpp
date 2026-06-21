#include <iostream>
#include <vector>
#include <random>
#include <cmath>

// External Libraries
#include <SpiceUsr.h>
#include <ceres/ceres.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

// ============================================================================
// 1. DYNAMICS & INTEGRATION (Templated for Ceres AutoDiff)
// ============================================================================

template <typename T>
void SunGravity(const T state[6], T state_dot[6]) {
    // Sun Gravitational Parameter (km^3 / s^2)
    T mu_sun = T(1.327124400419393e11); 

    T x = state[0], y = state[1], z = state[2];
    T r2 = x*x + y*y + z*z;
    T r3 = r2 * ceres::sqrt(r2);

    state_dot[0] = state[3];
    state_dot[1] = state[4];
    state_dot[2] = state[5];
    state_dot[3] = -mu_sun * x / r3;
    state_dot[4] = -mu_sun * y / r3;
    state_dot[5] = -mu_sun * z / r3;
}

template <typename T>
void RK4_Step(T state[6], double dt) {
    T k1[6], k2[6], k3[6], k4[6], temp[6];
    T dt_T = T(dt);

    SunGravity(state, k1);
    for(int i=0; i<6; ++i) temp[i] = state[i] + T(0.5) * dt_T * k1[i];
    
    SunGravity(temp, k2);
    for(int i=0; i<6; ++i) temp[i] = state[i] + T(0.5) * dt_T * k2[i];
    
    SunGravity(temp, k3);
    for(int i=0; i<6; ++i) temp[i] = state[i] + dt_T * k3[i];
    
    SunGravity(temp, k4);
    for(int i=0; i<6; ++i) {
        state[i] += (dt_T / T(6.0)) * (k1[i] + T(2.0)*k2[i] + T(2.0)*k3[i] + k4[i]);
    }
}

// ============================================================================
// 2. CERES COST FUNCTOR (Range Observation Model)
// ============================================================================

struct RangeResidual {
    RangeResidual(double t_obs, double t_epoch, double measured_range, const double station_pos[3])
        : t_obs_(t_obs), t_epoch_(t_epoch), range_(measured_range) {
        station_[0] = station_pos[0]; station_[1] = station_pos[1]; station_[2] = station_pos[2];
    }

    template <typename T>
    bool operator()(const T* const initial_state, T* residual) const {
        T current_state[6];
        for(int i=0; i<6; ++i) current_state[i] = initial_state[i];

        // Propagate state to observation time
        double dt = 3600.0; 
        double t_current = t_epoch_;
        while (t_current < t_obs_) {
            double step = std::min(dt, t_obs_ - t_current);
            RK4_Step(current_state, step);
            t_current += step;
        }

        // Compute Range
        T dx = current_state[0] - T(station_[0]);
        T dy = current_state[1] - T(station_[1]);
        T dz = current_state[2] - T(station_[2]);
        T computed_range = ceres::sqrt(dx*dx + dy*dy + dz*dz);

        // Residual (Y - G(X))
        residual[0] = T(range_) - computed_range;
        return true;
    }

private:
    double t_obs_, t_epoch_, range_;
    double station_[3];
};

// ============================================================================
// 3. UTILITIES
// ============================================================================
void ensureNoSpiceError(const char* msg) {
    if (failed_c()) {
        char err[1024];
        getmsg_c("LONG", 1024, err);
        std::cerr << "SPICE Error [" << msg << "]: " << err << std::endl;
        reset_c();
        exit(1);
    }
}

// ============================================================================
// 4. MAIN APPLICATION
// ============================================================================

int main() {
    // ---------------------------------------------------------
    // A. SPICE Initialization
    // ---------------------------------------------------------
    erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
    errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
    
    // Make sure your kernels.tm includes leapseconds, planets, and voyager
    furnsh_c("../kernels.tm");
    ensureNoSpiceError("Loading kernels");

    SpiceDouble et_start, et_end;
    str2et_c("1978-02-01T00:00:00", &et_start); 
    str2et_c("1990-04-01T00:00:00", &et_end);
    ensureNoSpiceError("Time conversion");

    // ---------------------------------------------------------
    // B. Generate True Trajectories for ImPlot (Grand Tour)
    // ---------------------------------------------------------
    const double SCALE = 1000000.0; // Scale to millions of km for UI limits
    int steps = 2000;
    double step_size = (et_end - et_start) / steps;

    std::vector<double> e_x, e_y, j_x, j_y, s_x, s_y, v_x, v_y;
    SpiceDouble lt;

    std::cout << "Extracting Grand Tour ephemeris..." << std::endl;
    for (int i = 0; i <= steps; ++i) {
        SpiceDouble et = et_start + (i * step_size);
        SpiceDouble st_e[6], st_j[6], st_s[6], st_v[6];
        
        spkezr_c("EARTH", et, "J2000", "NONE", "SUN", st_e, &lt);
        spkezr_c("5",     et, "J2000", "NONE", "SUN", st_j, &lt);
        spkezr_c("6",     et, "J2000", "NONE", "SUN", st_s, &lt);
        spkezr_c("-31",   et, "J2000", "NONE", "SUN", st_v, &lt);

        e_x.push_back(st_e[0]/SCALE); e_y.push_back(st_e[1]/SCALE);
        j_x.push_back(st_j[0]/SCALE); j_y.push_back(st_j[1]/SCALE);
        s_x.push_back(st_s[0]/SCALE); s_y.push_back(st_s[1]/SCALE);
        v_x.push_back(st_v[0]/SCALE); v_y.push_back(st_v[1]/SCALE);
    }

    // ---------------------------------------------------------
    // C. Perform Orbit Determination (10 Days early in mission)
    // ---------------------------------------------------------
    std::cout << "Running WLS Orbit Determination..." << std::endl;
    SpiceDouble od_start = et_start;
    SpiceDouble od_end = od_start + (10 * 86400.0); // 10 days of tracking

    // Get exact initial state
    double true_initial_state[6];
    spkezr_c("-31", od_start, "J2000", "NONE", "SUN", true_initial_state, &lt);

    // Create Synthetic Observations
    std::vector<double> obs_times, obs_ranges;
    std::default_random_engine gen;
    std::normal_distribution<double> noise(0.0, 0.010); // 10m noise

    for (double et = od_start; et <= od_end; et += 86400.0) {
        double st_v[6], st_e[6];
        spkezr_c("-31", et, "J2000", "NONE", "SUN", st_v, &lt);
        spkezr_c("EARTH", et, "J2000", "NONE", "SUN", st_e, &lt); // Assuming station at Earth center for simplicity
        
        double dx = st_v[0] - st_e[0], dy = st_v[1] - st_e[1], dz = st_v[2] - st_e[2];
        obs_times.push_back(et);
        obs_ranges.push_back(sqrt(dx*dx + dy*dy + dz*dz) + noise(gen));
    }

    // Perturb initial state (Bad Guess)
    double estimated_state[6];
    for(int i=0; i<6; ++i) estimated_state[i] = true_initial_state[i];
    estimated_state[0] += 50000.0; // 50,000 km error
    estimated_state[3] += 0.5;     // 0.5 km/s velocity error

    // Setup and Solve Ceres Problem
    ceres::Problem problem;
    for (size_t i = 0; i < obs_times.size(); ++i) {
        double st_e[6];
        spkezr_c("EARTH", obs_times[i], "J2000", "NONE", "SUN", st_e, &lt);
        
        problem.AddResidualBlock(
            new ceres::AutoDiffCostFunction<RangeResidual, 1, 6>(
                new RangeResidual(obs_times[i], od_start, obs_ranges[i], st_e)
            ),
            new ceres::HuberLoss(1.0),
            estimated_state
        );
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.max_num_iterations = 100; // Let it run a bit longer!
    
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // ---------------------------------------------------------
    // D. GUI / Render Loop Setup
    // ---------------------------------------------------------
 // ---------------------------------------------------------
    // D. GUI / Render Loop Setup
    // ---------------------------------------------------------
        if (!glfwInit()) return -1;

        // --- MAC SPECIFIC OPENGL HINTS ---
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); 
    #endif
        // ---------------------------------

        GLFWwindow* window = glfwCreateWindow(1280, 800, "Deep Space OD Simulator", NULL, NULL);
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1); // Enable vsync

        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        
        // --- MAC SPECIFIC SHADER VERSION ---
        ImGui_ImplOpenGL3_Init("#version 330 core");

    // Main render loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowSize(ImVec2(1200, 750), ImGuiCond_FirstUseEver);
        ImGui::Begin("Statistical Orbit Determination");

        ImGui::Text("Voyager 1 Grand Tour (1978 - 1990)");
        ImGui::Text("Ceres Initial Error: X: %.2f km | Vx: %.4f km/s", 
            std::abs(true_initial_state[0] - estimated_state[0]),
            std::abs(true_initial_state[3] - estimated_state[3]));

        if (ImPlot::BeginPlot("Solar System J2000 (Millions of km)", ImVec2(-1, -1))) {
            ImPlot::SetupAxes("X (10^6 km)", "Y (10^6 km)", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
            ImPlot::SetupAxesLimits(-2000, 2000, -2000, 2000); // Frame Jupiter/Saturn scale

            // Plot Planetary Orbits
            ImPlot::PlotLine("Earth", e_x.data(), e_y.data(), steps);
            
            ImPlot::PlotLine("Jupiter", j_x.data(), j_y.data(), steps);
            
  
            ImPlot::PlotLine("Saturn", s_x.data(), s_y.data(), steps);

            // Plot Voyager True Trajectory

            ImPlot::PlotLine("Voyager 1 (True)", v_x.data(), v_y.data(), steps);

            ImPlot::EndPlot();
        }
        ImGui::End();

        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.05f, 0.05f, 0.05f, 1.00f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}