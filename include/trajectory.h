#pragma once

#include <array>

constexpr int STATE_DIM = 5;
using StateVec = std::array<double, STATE_DIM>;

enum StateIdx {
    IDX_V = 0,
    IDX_THETA = 1,
    IDX_X = 2,
    IDX_Y = 3,
    IDX_M = 4
};

struct TrajectoryRecord {
    double t;
    double V;
    double theta;
    double x;
    double y;
    double m;
    double L;
    double D;
    double T;
    double Ma;
    double alpha;
    double phase;
};

struct TrajectoryParams {
    double cruise_alt;
    double cruise_mach;
    double initial_mass;
    double fuel_mass;
    double cruise_alpha;
};

StateVec trajectoryDerivatives(const StateVec& state, double t, const TrajectoryParams& params);
bool checkTermination(const StateVec& state, double t, const TrajectoryParams& params,
                      double initial_fuel_mass);
TrajectoryParams defaultTrajectoryParams();
TrajectoryRecord evaluateTrajectory(const StateVec& state, double t, const TrajectoryParams& params);
