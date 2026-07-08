#ifndef CSPICE_WRAPPER_HPP
#define CSPICE_WRAPPER_HPP

#include <array>
#include <stdexcept>
#include <string>

// Tell C++ not to mangle the C function names
#ifdef __cplusplus
extern "C" {
#endif

#include "SpiceUsr.h"

#ifdef __cplusplus
}
#endif

// ---------------------------------------------------------
// Global CSPICE Error Handling Utilities
// Marked 'inline' so they can be included in multiple .cpp 
// files without violating the One Definition Rule (ODR).
// ---------------------------------------------------------

inline void throwSpiceError(const std::string& stage) {
    SpiceChar shortMsg[1841] = {0};
    SpiceChar longMsg[1841] = {0};

    getmsg_c("SHORT", sizeof(shortMsg), shortMsg);
    getmsg_c("LONG", sizeof(longMsg), longMsg);
    reset_c();

    throw std::runtime_error(stage + ": " + shortMsg + " | " + longMsg);
}

inline void ensureNoSpiceError(const std::string& stage) {
    if (failed_c()) {
        throwSpiceError(stage);
    }
}

namespace od {

class SpiceErrorActionGuard {
public:
    SpiceErrorActionGuard() {
        erract_c("GET", static_cast<SpiceInt>(previousAction_.size()), previousAction_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
    }

    SpiceErrorActionGuard(const SpiceErrorActionGuard&) = delete;
    SpiceErrorActionGuard& operator=(const SpiceErrorActionGuard&) = delete;

    ~SpiceErrorActionGuard() {
        erract_c("SET", 0, previousAction_.data());
    }

private:
    std::array<SpiceChar, 32> previousAction_ {};
};

class SpiceErrorModeGuard {
public:
    SpiceErrorModeGuard() {
        erract_c("GET", static_cast<SpiceInt>(previousAction_.size()), previousAction_.data());
        errprt_c("GET", static_cast<SpiceInt>(previousReport_.size()), previousReport_.data());
        erract_c("SET", 0, const_cast<SpiceChar*>("RETURN"));
        errprt_c("SET", 0, const_cast<SpiceChar*>("NONE"));
    }

    SpiceErrorModeGuard(const SpiceErrorModeGuard&) = delete;
    SpiceErrorModeGuard& operator=(const SpiceErrorModeGuard&) = delete;

    ~SpiceErrorModeGuard() {
        erract_c("SET", 0, previousAction_.data());
        errprt_c("SET", 0, previousReport_.data());
    }

private:
    std::array<SpiceChar, 32> previousAction_ {};
    std::array<SpiceChar, 32> previousReport_ {};
};

inline void throwIfSpiceFailed(const std::string& context) {
    if (!failed_c()) {
        return;
    }

    SpiceChar shortMessage[1841] = {0};
    SpiceChar longMessage[1841] = {0};
    getmsg_c("SHORT", sizeof(shortMessage), shortMessage);
    getmsg_c("LONG", sizeof(longMessage), longMessage);
    reset_c();

    throw std::runtime_error(context + ": " + shortMessage + " | " + longMessage);
}

} // namespace od

#endif // CSPICE_WRAPPER_HPP
