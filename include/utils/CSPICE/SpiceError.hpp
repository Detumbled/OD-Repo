#ifndef UTILS_CSPICE_SPICE_ERROR_HPP
#define UTILS_CSPICE_SPICE_ERROR_HPP

#include "utils/CSPICE/SpiceErrorModeGuard.hpp"

#include <SpiceUsr.h>

#include <stdexcept>
#include <string>

namespace od {

[[noreturn]] inline void throwSpiceError(const std::string& context) {
    SpiceChar shortMessage[1841] = {0};
    SpiceChar longMessage[1841] = {0};

    getmsg_c("SHORT", sizeof(shortMessage), shortMessage);
    getmsg_c("LONG", sizeof(longMessage), longMessage);
    reset_c();

    throw std::runtime_error(context + ": " + shortMessage + " | " + longMessage);
}

inline void throwIfSpiceFailed(const std::string& context) {
    if (failed_c()) {
        throwSpiceError(context);
    }
}

inline void ensureNoSpiceError(const std::string& context) {
    throwIfSpiceFailed(context);
}

} // namespace od

#endif // UTILS_CSPICE_SPICE_ERROR_HPP
