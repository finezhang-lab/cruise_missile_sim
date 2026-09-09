// ============================================================================
// trajectory.cpp -- 导弹质点弹道运动方程实现
// 二维纵向质点运动方程:
//   m * dV/dt     = T*cos(alpha) - D - m*g*sin(theta)
//   m*V * dtheta/dt = T*sin(alpha) + L - m*g*cos(theta)
//   dx/dt         = V * cos(theta)
//   dy/dt         = V * sin(theta)
//   dm/dt         = -massFlowRate
// ============================================================================

#include "trajectory.h"
#include "config.h"
#include "atmosphere.h"
#include "aerodynamics.h"
#include "propulsion.h"
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// 默认工况参数 (基准工况 C-1: 海上巡航)
// ---------------------------------------------------------------------------
TrajectoryParams defaultTrajectoryParams() {
    TrajectoryParams p;
    p.cruise_alt   = config::CRUISE_ALT_SEA; // 15 m
    p.cruise_mach  = config::CRUISE_MACH;     // 0.72
    p.initial_mass = config::MASS_LAUNCH;     // 1450 kg
    p.fuel_mass    = config::MASS_FUEL;       // 450 kg
    p.cruise_alpha = 0.0;                     // 自动计算配平攻角
    return p;
}

// ---------------------------------------------------------------------------
// 右端函数: 计算状态向量的时间导数
// ---------------------------------------------------------------------------
StateVec trajectoryDerivatives(const StateVec& state, double t, const TrajectoryParams& params) {
    double V     = state[IDX_V];
    double theta = state[IDX_THETA];
    double y     = state[IDX_Y];
    double m     = state[IDX_M];

    // 防止速度和质量出现非物理值
    V = std::max(V, 1.0);
    m = std::max(m, config::MASS_EMPTY);

    // 大气参数
    AtmosphereState atm = calcAtmosphere(std::max(y, 0.0));
    double Ma = V / atm.a;

    // 动压
    double q = 0.5 * atm.rho * V * V;

    // 飞行阶段
    FlightPhase phase = getFlightPhase(t, m);

    // 攻角确定
    double alpha;
    if (phase == FlightPhase::BOOST) {
        // 助推段: 采用程序攻角 (爬升转弯)
        // 简单策略: 前 1 s 保持 theta0, 之后线性转到巡航攻角
        if (t < 1.0) {
            alpha = 5.0 * config::DEG2RAD; // 初始拉起攻角
        } else {
            // 线性过渡到较小攻角
            double ratio = (t - 1.0) / (config::BOOSTER_BURN_TIME - 1.0);
            ratio = std::min(ratio, 1.0);
            alpha = (5.0 * (1.0 - ratio) + 2.0 * ratio) * config::DEG2RAD;
        }
    } else {
        // 巡航段
        if (params.cruise_alpha > 0.0) {
            alpha = params.cruise_alpha;
        } else {
            // 自动计算配平攻角: L = m*g
            double weight = m * config::g0;
            alpha = calcTrimAlpha(Ma, q, weight);
            // 限制攻角范围
            alpha = std::clamp(alpha, -2.0 * config::DEG2RAD, 12.0 * config::DEG2RAD);
        }
    }

    // 气动力
    AeroResult aero = calcAerodynamics(alpha, Ma, q);
    double L = aero.L;
    double D = aero.D;

    // 推力
    PropulsionResult prop = calcPropulsion(t, y, V, Ma, m);
    double T = prop.thrust;
    double m_dot = prop.massFlowRate;

    // 助推器壳体分离在 integrator.cpp 中处理 (质量跳变)

    // ---- 运动方程右端函数 ----
    StateVec dstate;

    // dV/dt = [T*cos(alpha) - D - m*g*sin(theta)] / m
    dstate[IDX_V] = (T * std::cos(alpha) - D - m * config::g0 * std::sin(theta)) / m;

    // dtheta/dt = [T*sin(alpha) + L - m*g*cos(theta)] / (m*V)
    dstate[IDX_THETA] = (T * std::sin(alpha) + L - m * config::g0 * std::cos(theta)) / (m * V);

    // dx/dt = V * cos(theta)
    dstate[IDX_X] = V * std::cos(theta);

    // dy/dt = V * sin(theta)
    dstate[IDX_Y] = V * std::sin(theta);

    // dm/dt = -m_dot
    dstate[IDX_M] = -m_dot;

    return dstate;
}

// ---------------------------------------------------------------------------
// 终止条件检查
// ---------------------------------------------------------------------------
bool checkTermination(const StateVec& state, double t, const TrajectoryParams& params,
                      double initial_fuel_mass) {
    double y = state[IDX_Y];
    double m = state[IDX_M];

    // 条件 1: 超时
    if (t >= config::MAX_SIM_TIME) return true;

    // 条件 2: 触地 (高度 < 0)
    if (y < 0.0 && t > config::BOOSTER_BURN_TIME + 5.0) return true;

    // 条件 3: 巡航段燃油耗尽
    // 助推段结束后，导弹质量 = 发射质量 - 助推器(含壳体) - 已消耗燃油
    double mass_without_booster = config::MASS_LAUNCH - config::MASS_BOOSTER;
    double fuel_consumed_cruise = mass_without_booster - m;

    if (t > config::BOOSTER_BURN_TIME && fuel_consumed_cruise >= params.fuel_mass) {
        return true;
    }

    return false;
}
