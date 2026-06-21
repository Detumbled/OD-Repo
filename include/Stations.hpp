#ifndef STATIONS_HPP
#define STATIONS_HPP

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>

namespace od {

class Station {
public:
    using Vector3 = std::array<double, 3>;

    Station() = default;

    Station(std::string name,
            const Vector3& ecef_position_km,
            const Vector3& ecef_bias_km = {0.0, 0.0, 0.0},
            const Vector3& ecef_noise_sigma_km = {0.0, 0.0, 0.0},
            double latitude_rad = 0.0,
            double longitude_rad = 0.0,
            double altitude_km = 0.0,
            double elevation_mask_rad = 0.0,
            bool enabled = true)
        : name_(std::move(name)),
          ecef_position_km_(ecef_position_km),
          ecef_bias_km_(ecef_bias_km),
          ecef_noise_sigma_km_(ecef_noise_sigma_km),
          latitude_rad_(latitude_rad),
          longitude_rad_(longitude_rad),
          altitude_km_(altitude_km),
          elevation_mask_rad_(elevation_mask_rad),
          enabled_(enabled) {
        validateNoiseSigma();
    }

    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] const Vector3& ecefPositionKm() const noexcept { return ecef_position_km_; }
    [[nodiscard]] const Vector3& ecefBiasKm() const noexcept { return ecef_bias_km_; }
    [[nodiscard]] const Vector3& ecefNoiseSigmaKm() const noexcept { return ecef_noise_sigma_km_; }
    [[nodiscard]] double latitudeRad() const noexcept { return latitude_rad_; }
    [[nodiscard]] double longitudeRad() const noexcept { return longitude_rad_; }
    [[nodiscard]] double altitudeKm() const noexcept { return altitude_km_; }
    [[nodiscard]] double elevationMaskRad() const noexcept { return elevation_mask_rad_; }
    [[nodiscard]] bool enabled() const noexcept { return enabled_; }

    void setName(std::string name) { name_ = std::move(name); }

    void setEcefPositionKm(const Vector3& ecef_position_km) noexcept {
        ecef_position_km_ = ecef_position_km;
    }

    void setEcefBiasKm(const Vector3& ecef_bias_km) noexcept {
        ecef_bias_km_ = ecef_bias_km;
    }

    void setEcefNoiseSigmaKm(const Vector3& ecef_noise_sigma_km) {
        ecef_noise_sigma_km_ = ecef_noise_sigma_km;
        validateNoiseSigma();
    }

    void setGeodeticCoordinates(double latitude_rad,
                                double longitude_rad,
                                double altitude_km) noexcept {
        latitude_rad_ = latitude_rad;
        longitude_rad_ = longitude_rad;
        altitude_km_ = altitude_km;
    }

    void setElevationMaskRad(double elevation_mask_rad) noexcept {
        elevation_mask_rad_ = elevation_mask_rad;
    }

    void setEnabled(bool enabled) noexcept { enabled_ = enabled; }

    [[nodiscard]] Vector3 biasedPositionKm() const noexcept {
        return {
            ecef_position_km_[0] + ecef_bias_km_[0],
            ecef_position_km_[1] + ecef_bias_km_[1],
            ecef_position_km_[2] + ecef_bias_km_[2]
        };
    }


    [[nodiscard]] bool hasBias() const noexcept {
        return !isNearZero(ecef_bias_km_);
    }

    [[nodiscard]] bool hasNoise() const noexcept {
        return !isNearZero(ecef_noise_sigma_km_);
    }

private:
    static constexpr double kZeroTolerance = 1.0e-15;

    [[nodiscard]] static bool isNearZero(const Vector3& vector) noexcept {
        return std::abs(vector[0]) < kZeroTolerance
            && std::abs(vector[1]) < kZeroTolerance
            && std::abs(vector[2]) < kZeroTolerance;
    }

    void validateNoiseSigma() const {
        for (double sigma : ecef_noise_sigma_km_) {
            if (sigma < 0.0) {
                throw std::invalid_argument("Station noise sigma must be non-negative.");
            }
        }
    }

    std::string name_;
    Vector3 ecef_position_km_ {0.0, 0.0, 0.0};
    Vector3 ecef_bias_km_ {0.0, 0.0, 0.0};
    Vector3 ecef_noise_sigma_km_ {0.0, 0.0, 0.0};
    double latitude_rad_ {0.0};
    double longitude_rad_ {0.0};
    double altitude_km_ {0.0};
    double elevation_mask_rad_ {0.0};
    bool enabled_ {true};
};

} // namespace od

#endif // STATIONS_HPP
