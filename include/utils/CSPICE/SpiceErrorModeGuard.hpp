#ifndef UTILS_CSPICE_SPICE_ERROR_MODE_GUARD_HPP
#define UTILS_CSPICE_SPICE_ERROR_MODE_GUARD_HPP

#include <SpiceUsr.h>

#include <array>

namespace od {

class SpiceErrorModeGuard {
public:
    SpiceErrorModeGuard()
        : SpiceErrorModeGuard("RETURN", nullptr) {
    }

    explicit SpiceErrorModeGuard(const char* action)
        : SpiceErrorModeGuard(action, nullptr) {
    }

    SpiceErrorModeGuard(const char* action, const char* report) {
        erract_c("GET", static_cast<SpiceInt>(previousAction_.size()), previousAction_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>(action));

        if (report != nullptr) {
            restoreReport_ = true;
            errprt_c("GET", static_cast<SpiceInt>(previousReport_.size()), previousReport_.data());
            errprt_c("SET", 0, const_cast<SpiceChar*>(report));
        }
    }

    SpiceErrorModeGuard(const SpiceErrorModeGuard&) = delete;
    SpiceErrorModeGuard& operator=(const SpiceErrorModeGuard&) = delete;

    ~SpiceErrorModeGuard() {
        erract_c("SET", 0, previousAction_.data());
        if (restoreReport_) {
            errprt_c("SET", 0, previousReport_.data());
        }
    }

private:
    std::array<SpiceChar, 32> previousAction_ {};
    std::array<SpiceChar, 256> previousReport_ {};
    bool restoreReport_ {false};
};

} // namespace od

#endif // UTILS_CSPICE_SPICE_ERROR_MODE_GUARD_HPP
