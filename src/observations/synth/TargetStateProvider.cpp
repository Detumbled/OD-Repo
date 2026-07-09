#include "observations/synth/TargetStateProvider.hpp"

#include "utils/CSPICE/SpiceError.hpp"

#include <SpiceUsr.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

using od::throwIfSpiceFailed;

} // namespace

InterpolatedTargetStateProvider::InterpolatedTargetStateProvider(
    const od::EphemerisInterpolator& ephemeris,
    std::string context)
    : ephemeris_(&ephemeris),
      context_(std::move(context)) {
}

TargetState6 InterpolatedTargetStateProvider::stateAt(double epochTdb) const {
    if (ephemeris_ == nullptr) {
        throw std::runtime_error("Interpolated target state provider has no ephemeris.");
    }
    if (!std::isfinite(epochTdb)) {
        throw std::invalid_argument("Target provider epoch must be finite.");
    }

    const Eigen::VectorXd interpolated = ephemeris_->interpolate(epochTdb);
    if (interpolated.size() != 6 || !interpolated.allFinite()) {
        throw std::runtime_error(context_ + " produced an invalid target state.");
    }

    TargetState6 state;
    state = interpolated;
    return state;
}

SpiceTargetStateProvider::SpiceTargetStateProvider(std::string target,
                                                   std::string frame,
                                                   std::string observer)
    : target_(std::move(target)),
      frame_(std::move(frame)),
      observer_(std::move(observer)) {
    if (target_.empty()) {
        throw std::invalid_argument("SPICE target state provider target cannot be empty.");
    }
    if (frame_.empty()) {
        throw std::invalid_argument("SPICE target state provider frame cannot be empty.");
    }
    if (observer_.empty()) {
        throw std::invalid_argument("SPICE target state provider observer cannot be empty.");
    }
}

TargetState6 SpiceTargetStateProvider::stateAt(double epochTdb) const {
    if (!std::isfinite(epochTdb)) {
        throw std::invalid_argument("SPICE target provider epoch must be finite.");
    }

    od::SpiceErrorModeGuard action_guard;

    SpiceDouble spice_state[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    SpiceDouble light_time = 0.0;
    spkezr_c(target_.c_str(),
             epochTdb,
             frame_.c_str(),
             "NONE",
             observer_.c_str(),
             spice_state,
             &light_time);
    throwIfSpiceFailed("Failed to compute SPICE target provider state for " + target_);

    TargetState6 state;
    for (Eigen::Index i = 0; i < state.size(); ++i) {
        state[i] = spice_state[i];
    }
    return state;
}

} // namespace fd::observations::synth
