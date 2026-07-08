#include "filters/BatchLeastSquaresDriver.hpp"
#include "filters/WLS.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {

constexpr double kRangeSigmaKm = 0.010;
constexpr double kArcDurationSec = 8.0 * 3600.0;
constexpr double kCadenceSec = 300.0;

struct SyntheticBatch {
    Eigen::VectorXd residuals;
    Eigen::VectorXd observed_ranges;
    Eigen::VectorXd epochs_sec;
    Eigen::MatrixXd design_matrix;
    Eigen::MatrixXd measurement_covariance;
    Eigen::VectorXd nonlinear_prefit_residuals;
    Eigen::Vector3d station_position_km;
    Eigen::VectorXd true_state;
    Eigen::VectorXd prior_state;
};

[[nodiscard]] Eigen::Vector3d positionAt(const Eigen::VectorXd& epochState,
                                         double dt) {
    return epochState.segment<3>(0) + dt * epochState.segment<3>(3);
}

[[nodiscard]] double rangeFromStation(const Eigen::VectorXd& epochState,
                                      const Eigen::Vector3d& stationPositionKm,
                                      double dt) {
    return (positionAt(epochState, dt) - stationPositionKm).norm();
}

[[nodiscard]] SyntheticBatch generateSyntheticRangeBatch() {
    // In an operational test this is where the caller would load kernels:
    // furnsh_c("kernels.tm");
    // This unit-style test keeps station geometry deterministic and kernel-free.
    const Eigen::Vector3d station_position_km(6378.1363, -120.0, 35.0);

    Eigen::VectorXd true_state(6);
    true_state << 7078.0, -240.0, 1260.0,
                  0.12, 7.35, 1.05;

    Eigen::VectorXd prior_error(6);
    prior_error << 5.0, -3.0, 2.0,
                   0.003, -0.004, 0.002;

    const Eigen::VectorXd prior_state = true_state + prior_error;
    const int observation_count = static_cast<int>(kArcDurationSec / kCadenceSec) + 1;

    Eigen::VectorXd residuals(observation_count);
    Eigen::VectorXd observed_ranges(observation_count);
    Eigen::VectorXd epochs_sec(observation_count);
    Eigen::VectorXd nonlinear_prefit_residuals(observation_count);
    Eigen::MatrixXd design_matrix(observation_count, 6);

    std::mt19937 generator(421337);
    std::normal_distribution<double> noise_distribution(0.0, kRangeSigmaKm);

    for (int row = 0; row < observation_count; ++row) {
        const double dt = static_cast<double>(row) * kCadenceSec;

        const Eigen::Vector3d prior_position = positionAt(prior_state, dt);
        const Eigen::Vector3d line_of_sight = prior_position - station_position_km;
        const double prior_range = line_of_sight.norm();
        const double true_range = rangeFromStation(true_state, station_position_km, dt);
        const double observed_range = true_range + noise_distribution(generator);

        observed_ranges[row] = observed_range;
        epochs_sec[row] = dt;
        nonlinear_prefit_residuals[row] = observed_range - prior_range;
        residuals[row] = nonlinear_prefit_residuals[row];

        const Eigen::Vector3d range_unit = line_of_sight / prior_range;
        design_matrix.block<1, 3>(row, 0) = range_unit.transpose();
        design_matrix.block<1, 3>(row, 3) = dt * range_unit.transpose();
    }

    SyntheticBatch batch;
    batch.residuals = std::move(residuals);
    batch.observed_ranges = std::move(observed_ranges);
    batch.epochs_sec = std::move(epochs_sec);
    batch.design_matrix = std::move(design_matrix);
    batch.measurement_covariance = Eigen::MatrixXd::Identity(observation_count, observation_count)
        * (kRangeSigmaKm * kRangeSigmaKm);
    batch.nonlinear_prefit_residuals = std::move(nonlinear_prefit_residuals);
    batch.station_position_km = station_position_km;
    batch.true_state = std::move(true_state);
    batch.prior_state = std::move(prior_state);
    return batch;
}

[[nodiscard]] Eigen::VectorXd computeRangeResiduals(const Eigen::VectorXd& observedRanges,
                                                    const Eigen::VectorXd& epochsSec,
                                                    const Eigen::Vector3d& stationPositionKm,
                                                    const Eigen::VectorXd& epochState) {
    Eigen::VectorXd residuals(observedRanges.size());

    for (Eigen::Index row = 0; row < observedRanges.size(); ++row) {
        residuals[row] = observedRanges[row]
            - rangeFromStation(epochState, stationPositionKm, epochsSec[row]);
    }

    return residuals;
}

[[nodiscard]] fd::filters::BatchLeastSquaresDriver::LinearizedProblem linearizeRangeBatch(
    const SyntheticBatch& batch,
    const Eigen::VectorXd& epochState) {
    fd::filters::BatchLeastSquaresDriver::LinearizedProblem problem;
    problem.residuals.resize(batch.observed_ranges.size());
    problem.designMatrix.resize(batch.observed_ranges.size(), epochState.size());
    problem.measurementVariances =
        Eigen::VectorXd::Constant(batch.observed_ranges.size(), kRangeSigmaKm * kRangeSigmaKm);

    for (Eigen::Index row = 0; row < batch.observed_ranges.size(); ++row) {
        const double dt = batch.epochs_sec[row];
        const Eigen::Vector3d position = positionAt(epochState, dt);
        const Eigen::Vector3d line_of_sight = position - batch.station_position_km;
        const double range = line_of_sight.norm();
        const Eigen::Vector3d range_unit = line_of_sight / range;

        problem.residuals[row] = batch.observed_ranges[row] - range;
        problem.designMatrix.row(row).setZero();
        problem.designMatrix.block<1, 3>(row, 0) = range_unit.transpose();
        problem.designMatrix.block<1, 3>(row, 3) = dt * range_unit.transpose();
    }

    return problem;
}

[[nodiscard]] double rms(const Eigen::VectorXd& values) {
    return std::sqrt(values.squaredNorm() / static_cast<double>(values.size()));
}

void printStateDelta(const std::string& label,
                     const Eigen::VectorXd& estimate,
                     const Eigen::VectorXd& truth) {
    const Eigen::VectorXd error = estimate - truth;
    std::cout << label << '\n'
              << "  position error norm : " << error.segment<3>(0).norm() << " km\n"
              << "  velocity error norm : " << error.segment<3>(3).norm() << " km/s\n"
              << "  full error          : " << error.transpose() << '\n';
}

} // namespace

int main() {
    try {
        const SyntheticBatch batch = generateSyntheticRangeBatch();

        Eigen::MatrixXd prior_covariance = Eigen::MatrixXd::Zero(6, 6);
        prior_covariance.diagonal() << 100.0, 100.0, 100.0,
                                       1.0e-4, 1.0e-4, 1.0e-4;

        fd::filters::WLS filter;
        filter.setInitialState(batch.prior_state, prior_covariance, 0.0);
        filter.processBatch(batch.residuals, batch.design_matrix, batch.measurement_covariance);

        fd::filters::WLS streaming_filter;
        streaming_filter.setInitialState(batch.prior_state, prior_covariance, 0.0);
        for (Eigen::Index row = 0; row < batch.residuals.size(); ++row) {
            streaming_filter.addMeasurement(batch.residuals[row],
                                            batch.design_matrix.row(row),
                                            kRangeSigmaKm);
        }
        const Eigen::VectorXd streaming_correction = streaming_filter.solveAndUpdate();
        const Eigen::VectorXd streaming_postfit_residuals =
            batch.residuals - batch.design_matrix * streaming_correction;

        const Eigen::VectorXd nonlinear_postfit_residuals = computeRangeResiduals(batch.observed_ranges,
                                                                                  batch.epochs_sec,
                                                                                  batch.station_position_km,
                                                                                  filter.state());

        fd::filters::WLS driver_filter;
        driver_filter.setInitialState(batch.prior_state, prior_covariance, 0.0);
        fd::filters::BatchLeastSquaresDriver::ConvergenceCriteria driver_criteria;
        driver_criteria.maxIterations = 8;
        driver_criteria.correctionBlockTolerances = {
            {0, 3, 1.0e-4},
            {3, 3, 1.0e-8}
        };
        fd::filters::BatchLeastSquaresDriver driver(
            driver_filter,
            [&batch](const Eigen::VectorXd& state) {
                return linearizeRangeBatch(batch, state);
            },
            driver_criteria);
        const fd::filters::BatchLeastSquaresDriver::Result driver_result = driver.solve();
        const Eigen::VectorXd driver_postfit_residuals =
            computeRangeResiduals(batch.observed_ranges,
                                  batch.epochs_sec,
                                  batch.station_position_km,
                                  driver_filter.state());

        const double prefit_rms_km = rms(batch.nonlinear_prefit_residuals);
        const double postfit_rms_km = rms(filter.lastPostfitResiduals());
        const double nonlinear_postfit_rms_km = rms(nonlinear_postfit_residuals);
        const double driver_postfit_rms_km = rms(driver_postfit_residuals);
        const double prior_position_error_km = (batch.prior_state.segment<3>(0)
                                                - batch.true_state.segment<3>(0)).norm();
        const double posterior_position_error_km = (filter.state().segment<3>(0)
                                                    - batch.true_state.segment<3>(0)).norm();
        const double driver_position_error_km = (driver_filter.state().segment<3>(0)
                                                - batch.true_state.segment<3>(0)).norm();

        std::cout << std::fixed << std::setprecision(9)
                  << "Synthetic WLS batch test\n"
                  << "  observations          : " << batch.residuals.size() << '\n'
                  << "  streamed measurements : " << streaming_filter.accumulatedMeasurementCount() << '\n'
                  << "  range sigma           : " << kRangeSigmaKm * 1000.0 << " m\n"
                  << "  pre-fit O-C RMS       : " << prefit_rms_km * 1000.0 << " m\n"
                  << "  post-fit O-C RMS      : " << postfit_rms_km * 1000.0 << " m\n"
                  << "  nonlinear post RMS    : " << nonlinear_postfit_rms_km * 1000.0 << " m\n"
                  << "  driver status         : " << driver_result.statusName() << '\n'
                  << "  driver iterations     : " << driver_result.iterationCount() << '\n'
                  << "  driver post RMS       : " << driver_postfit_rms_km * 1000.0 << " m\n"
                  << "  state correction      : " << filter.lastStateCorrection().transpose() << '\n';

        printStateDelta("Prior state error", batch.prior_state, batch.true_state);
        printStateDelta("Posterior state error", filter.state(), batch.true_state);

        if (postfit_rms_km >= 0.25 * prefit_rms_km) {
            std::cerr << "FAIL: WLS did not sufficiently reduce linearized range residual RMS.\n";
            return EXIT_FAILURE;
        }

        if (nonlinear_postfit_rms_km >= 0.25 * prefit_rms_km) {
            std::cerr << "FAIL: WLS did not sufficiently reduce nonlinear range residual RMS.\n";
            return EXIT_FAILURE;
        }

        if (posterior_position_error_km >= prior_position_error_km) {
            std::cerr << "FAIL: WLS did not reduce the epoch position error.\n";
            return EXIT_FAILURE;
        }

        if (!driver_result.converged()) {
            std::cerr << "FAIL: batch least-squares driver did not report convergence.\n";
            return EXIT_FAILURE;
        }

        if (driver_postfit_rms_km >= nonlinear_postfit_rms_km) {
            std::cerr << "FAIL: batch least-squares driver did not improve nonlinear postfit RMS.\n";
            return EXIT_FAILURE;
        }

        if (driver_position_error_km >= prior_position_error_km) {
            std::cerr << "FAIL: batch least-squares driver did not reduce the epoch position error.\n";
            return EXIT_FAILURE;
        }

        if (!filter.state().isApprox(streaming_filter.state(), 1.0e-10)
            || !filter.lastStateCorrection().isApprox(streaming_correction, 1.0e-10)
            || !filter.lastPostfitResiduals().isApprox(streaming_postfit_residuals, 1.0e-10)
            || !filter.lastInformationMatrix().isApprox(streaming_filter.informationMatrix(), 1.0e-7)) {
            std::cerr << "FAIL: sequential WLS accumulation does not match the batch wrapper.\n";
            return EXIT_FAILURE;
        }

        if (!filter.covariance().allFinite() || !filter.lastInformationMatrix().allFinite()) {
            std::cerr << "FAIL: WLS produced non-finite covariance or information matrix.\n";
            return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "WLS synthetic test failed with exception: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
