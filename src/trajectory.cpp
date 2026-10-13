#include "trajectory.h"
#include "config.h"
#include "atmosphere.h"
#include "aerodynamics.h"
#include "propulsion.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

TrajectoryParams defaultTrajectoryParams() {
    return {config::CRUISE_ALT_SEA, config::CRUISE_MACH,
            config::MASS_INITIAL, config::MASS_FUEL, 0.0};
}

TrajectoryRecord evaluateTrajectory(const StateVec& state, double t,
                                    const TrajectoryParams& params) {
    const double V = std::max(state[IDX_V], 1e-6);
    const double m = std::max(state[IDX_M], 1e-6);
    const double h = std::max(state[IDX_Y], 0.0);
    const auto atm = calcAtmosphere(h);
    const double Ma = V / atm.a;
    const double q = 0.5 * atm.rho * V * V;
    const double phase_mass = state[IDX_M] + config::MASS_INITIAL - params.initial_mass;
    const auto prop = calcPropulsion(t, h, V, Ma, phase_mass);
    const bool boost = prop.phase == FlightPhase::BOOST;
    double alpha = params.cruise_alpha;
    if (boost) {
        const double ratio = std::clamp((t - 1.0) / (config::BOOSTER_BURN_TIME - 1.0), 0.0, 1.0);
        alpha = (5.0 - 3.0 * ratio) * config::DEG2RAD;
    } else if (alpha == 0.0) {
        const double target_theta = std::clamp((params.cruise_alt - state[IDX_Y]) / (10.0 * V),
                                               -10.0 * config::DEG2RAD, 10.0 * config::DEG2RAD);
        const double theta_rate = (target_theta - state[IDX_THETA]) / 3.0;
        const double normal_force = m * (config::g0 * std::cos(state[IDX_THETA]) + V * theta_rate);
        double low = -2.0 * config::DEG2RAD;
        double high = 12.0 * config::DEG2RAD;
        alpha = std::clamp(calcTrimAlpha(Ma, q, normal_force), low, high);
        for (int i = 0; i < 35; ++i) {
            const auto aero = calcAerodynamics(alpha, Ma, q);
            if (aero.L + prop.thrust * std::sin(alpha) < normal_force) low = alpha;
            else high = alpha;
            alpha = 0.5 * (low + high);
        }
    }
    const auto aero = calcAerodynamics(alpha, Ma, q);
    return {t, state[IDX_V], state[IDX_THETA] * config::RAD2DEG,
            state[IDX_X], state[IDX_Y], state[IDX_M], aero.L, aero.D,
            prop.thrust, Ma, alpha * config::RAD2DEG, boost ? 0.0 : 1.0};
}

StateVec trajectoryDerivatives(const StateVec& state, double t,
                               const TrajectoryParams& params) {
    for (double value : state) {
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite trajectory state");
    }
    if (state[IDX_M] <= 0.0) throw std::invalid_argument("Mass must be positive");
    const auto rec = evaluateTrajectory(state, t, params);
    const double alpha = rec.alpha * config::DEG2RAD;
    const double V = std::max(state[IDX_V], 1e-6);
    const double theta = state[IDX_THETA];
    const double m = state[IDX_M];
    const double phase_mass = m + config::MASS_INITIAL - params.initial_mass;
    const auto prop = calcPropulsion(t, std::max(state[IDX_Y], 0.0), V, rec.Ma, phase_mass);
    return {(rec.T * std::cos(alpha) - rec.D) / m - config::g0 * std::sin(theta),
            (rec.T * std::sin(alpha) + rec.L - m * config::g0 * std::cos(theta)) / (m * V),
            state[IDX_V] * std::cos(theta), state[IDX_V] * std::sin(theta),
            -prop.massFlowRate};
}

bool checkTermination(const StateVec& state, double t, const TrajectoryParams& params,
                      double initial_fuel_mass) {
    if (!std::isfinite(t)) return true;
    for (double value : state) {
        if (!std::isfinite(value)) return true;
    }
    if (t >= config::MAX_SIM_TIME || state[IDX_Y] < 0.0 || state[IDX_V] <= 0.0
        || state[IDX_M] <= 0.0) return true;
    if (state[IDX_Y] == 0.0 && state[IDX_THETA] < 0.0) return true;
    const double phase_mass = state[IDX_M] + config::MASS_INITIAL - params.initial_mass;
    const bool cruise = t >= config::BOOSTER_BURN_TIME
        || getFlightPhase(t, phase_mass) != FlightPhase::BOOST;
    const double dry_mass = params.initial_mass - config::MASS_BOOSTER - initial_fuel_mass;
    return cruise && state[IDX_M] <= dry_mass + 1e-8;
}
