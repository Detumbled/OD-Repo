#ifndef CSPICE_WRAPPER_HPP
#define CSPICE_WRAPPER_HPP

#include <string>
#include <stdexcept>

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

#endif // CSPICE_WRAPPER_HPP