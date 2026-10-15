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
// 综合飞控计算: 包含助推段程序转弯、巡航段高度保持闭环及自动油门推力调节
// ---------------------------------------------------------------------------
FlightControlResult calcFlightControl(const StateVec& state, double t, const TrajectoryParams& params) {
    FlightControlResult res{};
    double V     = std::max(state[IDX_V], 1.0);
    double theta = state[IDX_THETA];
    double y     = state[IDX_Y];
    double m     = std::max(state[IDX_M], config::MASS_EMPTY);

    AtmosphereState atm = calcAtmosphere(std::max(y, 0.0));
    double Ma = V / atm.a;
    double q  = 0.5 * atm.rho * V * V;

    // 飞行阶段判断 (初速处于巡航高亚声速或质量已扣除助推器，均为巡航段)
    FlightPhase phase = (V > 150.0 || m < params.initial_mass - 1.0)
                      ? FlightPhase::CRUISE : getFlightPhase(t, m);

    // 1. 攻角控制 (姿态与高度导引律)
    double alpha = 0.0;
    if (phase == FlightPhase::BOOST) {
        // 助推段: 采用程序转弯导引
        // 发射仰角约 30 deg，助推时间内平滑引导至接近水平 (0 deg)
        const double theta_initial = 30.0 * config::DEG2RAD;
        double ratio = std::clamp(t / config::BOOSTER_BURN_TIME, 0.0, 1.0);
        double theta_target = theta_initial * (1.0 - ratio);

        // 姿态闭环: 攻角由目标倾角偏差控制，保证法向力克服低速重力下沉
        double alpha_cmd = 1.5 * (theta_target - theta) + 3.0 * config::DEG2RAD;
        alpha = std::clamp(alpha_cmd, -1.0 * config::DEG2RAD, 10.0 * config::DEG2RAD);
    } else {
        // 巡航段 / 过渡段
        if (params.cruise_alpha > 0.0) {
            alpha = params.cruise_alpha;
        } else {
            // 自动高度捕获与高度保持闭环
            // 1) 静态重力配平攻角 (L = m*g)
            double weight = m * config::g0;
            double alpha_trim = calcTrimAlpha(Ma, q, weight);

            // 2) 高度捕获导引: 目标弹道倾角 theta_cmd
            double h_target = params.cruise_alt;
            double h_err = h_target - y;
            double h_dot = V * std::sin(theta);

            // 限制俯冲角/爬升角在安全范围内 [-2.5 deg, +4.0 deg]
            double theta_cmd = std::clamp(0.005 * h_err - 0.02 * (h_dot / std::max(V, 50.0)),
                                          -2.5 * config::DEG2RAD, 4.0 * config::DEG2RAD);

            // 3) 攻角控制指令: 跟踪 theta_cmd
            // 攻角上限放宽至 11.0 deg (失速角 14 deg)，保证低速 (90~120 m/s) 过渡段升力充足
            double alpha_cmd = alpha_trim + 1.5 * (theta_cmd - theta);
            alpha = std::clamp(alpha_cmd, -1.0 * config::DEG2RAD, 11.0 * config::DEG2RAD);
        }
    }
    res.alpha = alpha;

    // 2. 气动力初算 (用于推力平衡)
    AeroResult aero = calcAerodynamics(alpha, Ma, q);

    // 3. 推力与油门调节 (Auto-throttle)
    PropulsionResult prop = calcPropulsion(t, y, V, Ma, m);
    double T     = prop.thrust;
    double m_dot = prop.massFlowRate;

    if (phase == FlightPhase::CRUISE) {
        double V_target = params.cruise_mach * atm.a;
        if (V >= V_target - 0.5) {
            // 已达巡航速度: 自动油门闭环，推力平衡气动阻力与重力分量，维持巡航马赫数
            double T_req = aero.D + m * config::g0 * std::sin(theta) + 15.0 * (V_target - V);
            double T_max = prop.thrust;
            double T_min = 0.2 * T_max;
            T = std::clamp(T_req, T_min, T_max);

            // 节流状态下的燃油消耗率等比例折算: m_dot = sfc * T
            if (T_max > 1.0) {
                m_dot = prop.massFlowRate * (T / T_max);
            }
        }
        // 若 V < V_target - 0.5，发动机维持全推力加速
    }

    res.thrust = T;
    res.massFlowRate = m_dot;
    return res;
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
    double q  = 0.5 * atm.rho * V * V;

    // 飞控与推力计算
    FlightControlResult fc = calcFlightControl(state, t, params);
    double alpha = fc.alpha;
    double T     = fc.thrust;
    double m_dot = fc.massFlowRate;

    // 气动力
    AeroResult aero = calcAerodynamics(alpha, Ma, q);
    double L = aero.L;
    double D = aero.D;

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

    // 条件 1: 超时保护 (放宽至 30000s 约 8.3h，避免低速长航时或大载油量工况被数值提前截断)
    double max_sim_time = 30000.0;
    if (t >= max_sim_time) return true;

    // 条件 2: 触地 (高度 < 0)，给予发射初期 0.1s 离架容差
    if (y < 0.0 && t > 0.1) return true;

    // 条件 3: 巡航段燃油耗尽
    // 动态解算真实的巡航起点总质量，杜绝超额质量扫描时的死锁
    double mass_without_booster = config::MASS_LAUNCH - config::MASS_BOOSTER;
    double initial_cruise_mass  = (params.initial_mass > mass_without_booster)
                                ? (params.initial_mass - config::MASS_BOOSTER)
                                : params.initial_mass;
    double fuel_consumed_cruise = initial_cruise_mass - m;

    if (t > config::BOOSTER_BURN_TIME && fuel_consumed_cruise >= (params.fuel_mass - 1e-4)) {
        return true;
    }

    return false;
}
