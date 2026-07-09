#ifndef CSPICE_WRAPPER_HPP
#define CSPICE_WRAPPER_HPP

#include "utils/CSPICE/SpiceError.hpp"

#include <string>

inline void throwSpiceError(const std::string& stage) {
    od::throwSpiceError(stage);
}

inline void ensureNoSpiceError(const std::string& stage) {
    od::ensureNoSpiceError(stage);
}

#endif // CSPICE_WRAPPER_HPP
