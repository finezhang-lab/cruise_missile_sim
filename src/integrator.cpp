#include "integrator.h"
#include "config.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

void requireFinite(const StateVec& state) {
    for (double value : state) {
        if (!std::isfinite(value)) throw std::runtime_error("Non-finite integration state");
    }
}

StateVec shifted(const StateVec& state, const StateVec& slope, double scale) {
    StateVec result{};
    for (int i = 0; i < STATE_DIM; ++i) result[i] = state[i] + scale * slope[i];
    return result;
}

}

StateVec rk4Step(const StateVec& state, double t, double dt,
                 const TrajectoryParams& params, DerivativeFunc f) {
    if (!f || !std::isfinite(t) || !std::isfinite(dt) || dt <= 0.0 || t + dt == t) {
        throw std::invalid_argument("RK4 requires a derivative and a positive finite step");
    }
    requireFinite(state);
    const StateVec k1 = f(state, t, params);
    requireFinite(k1);
    const StateVec k2 = f(shifted(state, k1, dt / 2.0), t + dt / 2.0, params);
    requireFinite(k2);
    const StateVec k3 = f(shifted(state, k2, dt / 2.0), t + dt / 2.0, params);
    requireFinite(k3);
    const StateVec k4 = f(shifted(state, k3, dt), t + dt, params);
    requireFinite(k4);
    StateVec result{};
    for (int i = 0; i < STATE_DIM; ++i) {
        result[i] = state[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    requireFinite(result);
    return result;
}

std::vector<TrajectoryRecord> integrateTrajectory(
    const StateVec& initial_state, double t_start,
    const TrajectoryParams& params, DerivativeFunc f) {
    requireFinite(initial_state);
    const double values[] = {params.cruise_alt, params.cruise_mach, params.initial_mass,
                             params.fuel_mass, params.cruise_alpha, t_start};
    for (double value : values) {
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite trajectory parameter");
    }
    if (!f || t_start < 0.0 || params.cruise_alt < 0.0 || params.cruise_mach <= 0.0
        || params.fuel_mass < 0.0 || params.initial_mass <= config::MASS_BOOSTER + params.fuel_mass
        || initial_state[IDX_M] <= 0.0 || initial_state[IDX_V] <= 0.0) {
        throw std::invalid_argument("Invalid trajectory parameters or initial state");
    }
    StateVec state = initial_state;
    double t = t_start;
    const double post_boost_mass = params.initial_mass - config::MASS_BOOSTER;
    const double dry_mass = post_boost_mass - params.fuel_mass;
    bool separated = state[IDX_M] <= post_boost_mass + 1.0;
    if (!separated && t >= config::BOOSTER_BURN_TIME) {
        state[IDX_M] -= config::MASS_BOOSTER - config::MASS_BOOSTER_PROPELLANT;
        separated = true;
    }
    std::vector<TrajectoryRecord> records;
    records.push_back(evaluateTrajectory(state, t, params));
    int cruise_steps = 0;
    while (!checkTermination(state, t, params, params.fuel_mass)) {
        const bool boost = !separated && t < config::BOOSTER_BURN_TIME;
        double end = std::min(t + (boost ? config::DT_BOOST : config::DT_CRUISE),
                              config::MAX_SIM_TIME);
        if (boost) {
            const double boundaries[] = {0.5, 1.0, config::BOOSTER_BURN_TIME - 0.5,
                                          config::BOOSTER_BURN_TIME};
            for (double boundary : boundaries) {
                if (boundary > t && end >= boundary - 1e-12) end = boundary;
            }
        }
        const double dt = end - t;
        DerivativeFunc step_f = f;
        if (boost && end == config::BOOSTER_BURN_TIME) {
            step_f = [f, end](const StateVec& s, double time, const TrajectoryParams& p) {
                return f(s, time >= end ? std::nextafter(end, 0.0) : time, p);
            };
        }
        StateVec next = rk4Step(state, t, dt, params, step_f);
        const auto crossed = [boost, dry_mass](const StateVec& s) {
            return s[IDX_Y] < 0.0 || s[IDX_V] <= 0.0 || (!boost && s[IDX_M] <= dry_mass);
        };
        bool terminal = crossed(next);
        if (terminal) {
            double low = 0.0;
            double high = dt;
            for (int i = 0; i < 45; ++i) {
                const double middle = 0.5 * (low + high);
                const auto candidate = rk4Step(state, t, middle, params, step_f);
                if (crossed(candidate)) high = middle;
                else low = middle;
            }
            end = t + high;
            next = rk4Step(state, t, high, params, step_f);
            if (next[IDX_Y] < 0.0) next[IDX_Y] = 0.0;
            if (next[IDX_V] < 0.0) next[IDX_V] = 0.0;
            if (!boost && next[IDX_M] < dry_mass) next[IDX_M] = dry_mass;
        }
        state = next;
        t = end;
        if (boost || ++cruise_steps % 10 == 0 || terminal || t >= config::MAX_SIM_TIME) {
            records.push_back(evaluateTrajectory(state, t, params));
        }
        if (!separated && t == config::BOOSTER_BURN_TIME) {
            state[IDX_M] -= config::MASS_BOOSTER - config::MASS_BOOSTER_PROPELLANT;
            separated = true;
            records.push_back(evaluateTrajectory(state, t, params));
            cruise_steps = 0;
        }
        if (terminal) break;
    }
    if (records.back().t != t || records.back().m != state[IDX_M]) {
        records.push_back(evaluateTrajectory(state, t, params));
    }
    return records;
}
