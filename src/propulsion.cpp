// ============================================================================
// propulsion.cpp -- 动力系统模型实现（姚飞负责）
//
// 对应课题阶段二 B-1/B-2/B-3 与阶段三/四的动力接入。
// 模块只通过 propulsion.h 接口通信；本文件不得修改公共头文件签名（签名由张统一管理）。
//
// 模型层级（教学/课题抽象，不用于现实型号性能估算）：
//   1) 固体火箭助推器：梯形推力曲线（上升→恒推→下降），恒定比冲，
//      质量流率 m_dot = T / (Isp * g0)，保证 ∫m_dot dt = 推进剂质量。
//   2) 涡扇发动机：基准推力乘密度比修正与马赫数修正；SFC 乘高度/速度修正。
//   3) 阶段切换：助推段 t<=BOOSTER_BURN_TIME 且质量仍含助推器；
//      助推器分离为瞬时质量跳变，由 integrator.cpp 处理，本模块只返回
//      助推段 / 巡航段两种持续状态。
//
// 单位：全部 SI（N、kg、s、m、m/s；角度内部用 rad）。
// 参数来源：include/config.h（张管理）与 docs/planning 工况表；
//          助推器/涡扇/质量等参数为待核查参考值。
// ============================================================================

#include "propulsion.h"
#include "config.h"
#include "atmosphere.h"

#include <cmath>
#include <algorithm>

namespace {

// ---- 助推器梯形推力曲线几何参数 ----
constexpr double BOOST_RAMP_UP_S   = 0.5;   // 上升段时间 [s]
constexpr double BOOST_RAMP_DOWN_S = 0.5;   // 下降段时间 [s]

// ---- 涡扇推力修正经验系数 ----
constexpr double ENGINE_ALT_EXP    = 0.7;   // 密度比修正指数 σ^0.7
constexpr double ENGINE_SPEED_COEF = 0.05;  // 冲压效应系数 (1+0.05·Ma)

// ---- 涡扇耗油率修正经验系数 ----
constexpr double SFC_ALT_COEF    = 0.02;    // 高度修正 (1 - 0.02·h/1000)
constexpr double SFC_ALT_LIMIT_M = 11000.0; // 高度修正上限 [m]（对流层顶）
constexpr double SFC_SPEED_COEF  = 0.05;    // 速度修正 (1+0.05·Ma)

// 助推器梯形曲线总冲 [N·s]：梯形面积 = T_max * (t_b - (t_ramp_up + t_ramp_down)/2)
double boosterTotalImpulseNs() {
    const double t_b   = config::BOOSTER_BURN_TIME;
    const double T_max = config::BOOSTER_THRUST;
    return T_max * (t_b - (BOOST_RAMP_UP_S + BOOST_RAMP_DOWN_S) / 2.0);
}

// 助推器比冲 [s]：Isp = I_total / (m_propellant * g0)
double boosterIspS() {
    return boosterTotalImpulseNs()
         / (config::MASS_BOOSTER_PROPELLANT * config::g0);
}

} // namespace

// ---------------------------------------------------------------------------
// 飞行阶段判断
// 助推段：t <= BOOSTER_BURN_TIME 且质量仍包含助推器（> 分离后质量 + 1kg 容差）
// 巡航段：其余情况。SEPARATION 为瞬时事件，由 integrator.cpp 处理质量跳变。
// ---------------------------------------------------------------------------
FlightPhase getFlightPhase(double t, double mass) {
    const double mass_after_boost = config::MASS_LAUNCH - config::MASS_BOOSTER;
    if (t <= config::BOOSTER_BURN_TIME && mass > mass_after_boost + 1.0) {
        return FlightPhase::BOOST;
    }
    return FlightPhase::CRUISE;
}

// ---------------------------------------------------------------------------
// 固体火箭助推器推力（梯形曲线）
//   t ∈ [0, t_ramp_up]                 : T = T_max * t / t_ramp_up
//   t ∈ (t_ramp_up, t_b - t_ramp_down] : T = T_max
//   t ∈ (t_b - t_ramp_down, t_b]       : T = T_max * (t_b - t) / t_ramp_down
//   t < 0 或 t > t_b                   : T = 0
// ---------------------------------------------------------------------------
double calcBoosterThrust(double t) {
    const double T_max = config::BOOSTER_THRUST;
    const double t_b   = config::BOOSTER_BURN_TIME;

    if (t < 0.0 || t > t_b) {
        return 0.0;
    }
    if (t <= BOOST_RAMP_UP_S) {
        return T_max * (t / BOOST_RAMP_UP_S);
    }
    if (t <= t_b - BOOST_RAMP_DOWN_S) {
        return T_max;
    }
    return T_max * (t_b - t) / BOOST_RAMP_DOWN_S;
}

// ---------------------------------------------------------------------------
// 涡扇发动机推力模型
// 假设：config::ENGINE_THRUST 按“海平面静态基准推力”使用（与 config.h 中
//       “巡航推力”标注的口径差异为待确认项）。
//   T = T0 * σ^0.7 * (1 + 0.05·Ma),  σ = ρ(h)/ρ(0)
// 高度修正：进气密度比；速度修正：亚声速冲压效应（弱正贡献）。
// 参数 V 与 Ma 冗余（V = Ma·a），保留以匹配头文件签名。
// ---------------------------------------------------------------------------
double calcEngineThrust(double h, double V, double Ma) {
    (void)V; // V 与 Ma 冗余，签名兼容保留

    const double T0 = config::ENGINE_THRUST;

    const AtmosphereState atm    = calcAtmosphere(h);
    const AtmosphereState atm_sl = calcAtmosphere(0.0);

    const double sigma = atm.rho / atm_sl.rho;   // 密度比
    const double T     = T0 * std::pow(sigma, ENGINE_ALT_EXP);
    const double speed_factor = 1.0 + ENGINE_SPEED_COEF * Ma;

    return T * speed_factor;
}

// ---------------------------------------------------------------------------
// 涡扇发动机耗油率 [kg/(N·s)]
//   sfc = sfc0 * (1 - 0.02·min(h,11km)/1000) * (1 + 0.05·Ma)
// sfc0 = config::ENGINE_SFC（已换算为 kg/(N·s)）
// ---------------------------------------------------------------------------
double calcEngineSFC(double h, double V, double Ma) {
    (void)V; // V 与 Ma 冗余，签名兼容保留

    const double sfc0 = config::ENGINE_SFC;

    const double alt_factor = 1.0 - SFC_ALT_COEF
        * std::min(h, SFC_ALT_LIMIT_M) / 1000.0;
    const double spd_factor = 1.0 + SFC_SPEED_COEF * Ma;

    return sfc0 * alt_factor * spd_factor;
}

// ---------------------------------------------------------------------------
// 综合推力计算：按飞行阶段自动选择动力源
// 助推段：T = T_booster + T_engine；m_dot = m_dot_booster + m_dot_engine
//         其中 m_dot_booster = T_booster / (Isp·g0)，与推力曲线一致
// 巡航段：T = T_engine；m_dot = sfc * T_engine
// ---------------------------------------------------------------------------
PropulsionResult calcPropulsion(double t, double h, double V, double Ma, double mass) {
    PropulsionResult res;
    res.phase = getFlightPhase(t, mass);

    // 涡扇推力与耗油率（两个阶段都需要）
    const double T_engine     = calcEngineThrust(h, V, Ma);
    const double sfc          = calcEngineSFC(h, V, Ma);
    const double m_dot_engine = sfc * T_engine;

    if (res.phase == FlightPhase::BOOST) {
        const double T_booster = calcBoosterThrust(t);
        res.thrust = T_booster + T_engine;

        // 助推器质量流率与推力曲线同形，保证 ∫m_dot dt = 推进剂质量
        const double m_dot_booster = T_booster / (boosterIspS() * config::g0);

        res.massFlowRate = m_dot_booster + m_dot_engine;
    } else {
        res.thrust = T_engine;
        res.massFlowRate = m_dot_engine;
    }

    return res;
}
