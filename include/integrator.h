#pragma once

#include "trajectory.h"
#include <functional>
#include <vector>

using DerivativeFunc = std::function<StateVec(const StateVec&, double, const TrajectoryParams&)>;

StateVec rk4Step(const StateVec& state, double t, double dt,
                 const TrajectoryParams& params, DerivativeFunc f);
std::vector<TrajectoryRecord> integrateTrajectory(
    const StateVec& initial_state, double t_start,
    const TrajectoryParams& params, DerivativeFunc f);
