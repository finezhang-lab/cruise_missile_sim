// ============================================================================
// integrator.cpp -- RK4 数值积分器实现
// 经典四阶龙格-库塔法:
//   k1 = f(xn, tn)
//   k2 = f(xn + dt/2 * k1, tn + dt/2)
//   k3 = f(xn + dt/2 * k2, tn + dt/2)
//   k4 = f(xn + dt * k3, tn + dt)
//   xn+1 = xn + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
// ============================================================================

#include "integrator.h"
#include "config.h"
#include "atmosphere.h"
#include "aerodynamics.h"
#include "propulsion.h"
#include <cmath>
#include <algorithm>
#include <iostream>

// ---------------------------------------------------------------------------
// 状态向量加法
// ---------------------------------------------------------------------------
static StateVec stateAdd(const StateVec& a, const StateVec& b) {
    StateVec r;
    for (int i = 0; i < STATE_DIM; ++i)
        r[i] = a[i] + b[i];
    return r;
}

// ---------------------------------------------------------------------------
// 状态向量标量乘法
// ---------------------------------------------------------------------------
static StateVec stateScale(const StateVec& a, double s) {
    StateVec r;
    for (int i = 0; i < STATE_DIM; ++i)
        r[i] = a[i] * s;
    return r;
}

// ---------------------------------------------------------------------------
// 单步 RK4 积分
// ---------------------------------------------------------------------------
StateVec rk4Step(const StateVec& state, double t, double dt,
                 const TrajectoryParams& params,
                 DerivativeFunc f) {
    // k1 = f(xn, tn)
    StateVec k1 = f(state, t, params);

    // k2 = f(xn + dt/2 * k1, tn + dt/2)
    StateVec x2 = stateAdd(state, stateScale(k1, dt / 2.0));
    StateVec k2 = f(x2, t + dt / 2.0, params);

    // k3 = f(xn + dt/2 * k2, tn + dt/2)
    StateVec x3 = stateAdd(state, stateScale(k2, dt / 2.0));
    StateVec k3 = f(x3, t + dt / 2.0, params);

    // k4 = f(xn + dt * k3, tn + dt)
    StateVec x4 = stateAdd(state, stateScale(k3, dt));
    StateVec k4 = f(x4, t + dt, params);

    // xn+1 = xn + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    StateVec result;
    for (int i = 0; i < STATE_DIM; ++i) {
        result[i] = state[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return result;
}

// ---------------------------------------------------------------------------
// 构建单条 TrajectoryRecord
// ---------------------------------------------------------------------------
static TrajectoryRecord buildRecord(const StateVec& state, double t,
                                    const TrajectoryParams& params) {
    TrajectoryRecord rec;
    rec.t     = t;
    rec.V     = state[IDX_V];
    rec.theta = state[IDX_THETA] * config::RAD2DEG;
    rec.x     = state[IDX_X];
    rec.y     = state[IDX_Y];
    rec.m     = state[IDX_M];

    // 计算附加量
    AtmosphereState atm = calcAtmosphere(std::max(state[IDX_Y], 0.0));
    rec.Ma = state[IDX_V] / atm.a;
    double q = 0.5 * atm.rho * state[IDX_V] * state[IDX_V];

    double mass_without_booster = config::MASS_LAUNCH - config::MASS_BOOSTER;
    FlightPhase phase = (state[IDX_V] > 150.0 || state[IDX_M] <= mass_without_booster + 1.0)
                      ? FlightPhase::CRUISE : getFlightPhase(t, state[IDX_M]);

    if (params.enable_typical_profile) {
        if (phase == FlightPhase::BOOST) {
            rec.phase = 0.0; // 助推段
        } else if (state[IDX_X] < params.high_cruise_dist) {
            if (state[IDX_Y] < params.climb_target_alt - 100.0) {
                rec.phase = 1.0; // 爬升段
            } else {
                rec.phase = 2.0; // 5000m 高空平飞段
            }
        } else {
            if (state[IDX_Y] > params.sea_skim_alt + 30.0) {
                rec.phase = 3.0; // 俯冲下滑段
            } else {
                rec.phase = 4.0; // 15m 超低空掠海平飞段
            }
        }
    } else {
        rec.phase = (phase == FlightPhase::BOOST) ? 0.0 : 1.0;
    }

    // 使用统一的飞控与推力计算逻辑，确保记录输出与积分步进严格一致
    FlightControlResult fc = calcFlightControl(state, t, params);
    rec.alpha = fc.alpha * config::RAD2DEG;
    rec.T     = fc.thrust;

    AeroResult aero = calcAerodynamics(fc.alpha, rec.Ma, q);
    rec.L = aero.L;
    rec.D = aero.D;

    return rec;
}

// ---------------------------------------------------------------------------
// 全弹道积分
// 自动根据飞行阶段切换步长
// ---------------------------------------------------------------------------
std::vector<TrajectoryRecord> integrateTrajectory(
    const StateVec& initial_state,
    double t_start,
    const TrajectoryParams& params,
    DerivativeFunc f)
{
    std::vector<TrajectoryRecord> records;

    StateVec state = initial_state;
    double t = t_start;
    double fuel_at_cruise_start = params.fuel_mass;
    // 若初始质量已无助推器，或初速已处于巡航高亚声速，视为已分离
    double mass_without_booster = config::MASS_LAUNCH - config::MASS_BOOSTER;
    bool booster_separated = (state[IDX_M] <= mass_without_booster + 1.0) || (state[IDX_V] > 150.0);

    // 记录初始状态
    records.push_back(buildRecord(state, t, params));

    int cruise_step_count = 0; // 巡航段记录计数器

    while (!checkTermination(state, t, params, fuel_at_cruise_start)) {
        // 根据飞行阶段选择步长 (一旦处于巡航段，不可逆回助推段)
        FlightPhase phase = booster_separated ? FlightPhase::CRUISE : getFlightPhase(t, state[IDX_M]);
        double dt = (phase == FlightPhase::BOOST) ? config::DT_BOOST : config::DT_CRUISE;

        // 助推器分离: 在助推段结束瞬间抛掉壳体质量
        if (phase == FlightPhase::CRUISE && !booster_separated) {
            booster_separated = true;
            // 抛掉助推器壳体 (助推器总质量 - 推进剂质量 = 壳体质量)
            double shell_mass = config::MASS_BOOSTER - config::MASS_BOOSTER_PROPELLANT;
            state[IDX_M] -= shell_mass;
            std::cout << "  [助推器分离] t=" << t << " s, 抛掉壳体 "
                      << shell_mass << " kg, 当前质量 "
                      << state[IDX_M] << " kg\n";
        }

        // RK4 积分一步
        StateVec new_state = rk4Step(state, t, dt, params, f);

        // 质量下限保护
        new_state[IDX_M] = std::max(new_state[IDX_M], config::MASS_EMPTY);

        state = new_state;
        t += dt;

        // 每隔一定时间步记录一次 (避免数据量过大)
        // 助推段: 每步记录; 巡航段: 每 10 步记录一次
        if (phase == FlightPhase::BOOST) {
            records.push_back(buildRecord(state, t, params));
            cruise_step_count = 0;
        } else {
            cruise_step_count++;
            if (cruise_step_count % 10 == 0) {
                records.push_back(buildRecord(state, t, params));
            }
        }
    }

    // 记录终止状态
    records.push_back(buildRecord(state, t, params));

    return records;
}
