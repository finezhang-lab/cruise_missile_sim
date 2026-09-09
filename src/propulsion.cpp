// ============================================================================
// propulsion.cpp -- 动力系统模型实现
// 包含固体火箭助推器和涡扇发动机的推力/耗油特性
// ============================================================================

#include "propulsion.h"
#include "config.h"
#include "atmosphere.h"
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// 飞行阶段判断
// 助推段: t <= BOOSTER_BURN_TIME 且质量仍为发射质量级（含助推器）
// 巡航段: 助推器已分离，或时间已超过助推段
// ---------------------------------------------------------------------------
FlightPhase getFlightPhase(double t, double mass) {
    double mass_after_boost = config::MASS_LAUNCH - config::MASS_BOOSTER;
    if (t <= config::BOOSTER_BURN_TIME && mass > mass_after_boost + 1.0) {
        return FlightPhase::BOOST;
    }
    return FlightPhase::CRUISE;
}

// ---------------------------------------------------------------------------
// 固体火箭助推器推力 (恒推力近似)
// 实际推力曲线可近似为梯形:
//   - 0 ~ 0.5 s: 推力上升段 (线性增至额定)
//   - 0.5 ~ (tb-0.5) s: 恒定推力
//   - (tb-0.5) ~ tb s: 推力下降段 (线性降至 0)
// ---------------------------------------------------------------------------
double calcBoosterThrust(double t) {
    if (t < 0.0 || t > config::BOOSTER_BURN_TIME) {
        return 0.0;
    }

    double T_max = config::BOOSTER_THRUST;
    double t_ramp = 0.5;  // 上升段时间 [s]
    double t_decay = 0.5; // 下降段时间 [s]
    double t_b = config::BOOSTER_BURN_TIME;

    if (t <= t_ramp) {
        // 上升段
        return T_max * (t / t_ramp);
    } else if (t <= t_b - t_decay) {
        // 恒定推力段
        return T_max;
    } else {
        // 下降段
        return T_max * (t_b - t) / t_decay;
    }
}

// ---------------------------------------------------------------------------
// 涡扇发动机推力模型
// 基准推力: 2670 N (海平面静态)
// 推力随高度和速度的修正:
//   - 高度修正: 推力随空气密度降低而下降
//   - 速度修正: 冲压效应在亚声速段有微弱正贡献
// ---------------------------------------------------------------------------
double calcEngineThrust(double h, double V, double Ma) {
    double T0 = config::ENGINE_THRUST; // 海平面静态推力

    // 获取当前高度大气参数
    AtmosphereState atm = calcAtmosphere(h);
    AtmosphereState atm_sl = calcAtmosphere(0.0);

    // 密度比修正 (涡扇推力近似正比于进气密度)
    double sigma = atm.rho / atm_sl.rho;
    double T = T0 * std::pow(sigma, 0.7); // 经验指数 0.7

    // 速度修正 (冲压效应, 亚声速段微弱正贡献)
    // 参考: T = T_static * (1 + 0.1*Ma) for low subsonic
    double speed_factor = 1.0 + 0.05 * Ma;

    return T * speed_factor;
}

// ---------------------------------------------------------------------------
// 涡扇发动机耗油率 [kg/(N*s)]
// 基准: SFC ≈ 0.07 kg/(N*h) = 0.07/3600 kg/(N*s)
// 高度和速度对 SFC 有轻微影响
// ---------------------------------------------------------------------------
double calcEngineSFC(double h, double V, double Ma) {
    double sfc0 = config::ENGINE_SFC; // 基准耗油率 [kg/(N*s)]

    // 高度修正: 高空 SFC 略有下降 (热效率提高)
    // 经验公式: SFC = SFC0 * (1 - 0.02 * h/1000) for h < 11km
    double alt_factor = 1.0 - 0.02 * std::min(h, 11000.0) / 1000.0;

    // 速度修正: 高速 SFC 略有上升
    double spd_factor = 1.0 + 0.05 * Ma;

    return sfc0 * alt_factor * spd_factor;
}

// ---------------------------------------------------------------------------
// 综合推力计算
// 根据飞行阶段自动选择助推器或涡扇发动机
// ---------------------------------------------------------------------------
PropulsionResult calcPropulsion(double t, double h, double V, double Ma, double mass) {
    PropulsionResult res;
    res.phase = getFlightPhase(t, mass);

    if (res.phase == FlightPhase::BOOST) {
        // ---- 助推段 ----
        double T_booster = calcBoosterThrust(t);
        // 涡扇发动机在助推段也工作（助推器提供额外推力）
        double T_engine = calcEngineThrust(h, V, Ma);
        res.thrust = T_booster + T_engine;

        // 助推器质量流率: 推进剂总质量 / 工作时间
        double m_dot_booster = config::MASS_BOOSTER_PROPELLANT / config::BOOSTER_BURN_TIME;

        // 涡扇质量流率
        double sfc = calcEngineSFC(h, V, Ma);
        double m_dot_engine = sfc * T_engine;

        res.massFlowRate = m_dot_booster + m_dot_engine;
    } else {
        // ---- 巡航段 ----
        double T_engine = calcEngineThrust(h, V, Ma);
        res.thrust = T_engine;

        double sfc = calcEngineSFC(h, V, Ma);
        res.massFlowRate = sfc * T_engine;
    }

    return res;
}
