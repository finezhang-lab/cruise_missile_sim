// ============================================================================
// atmosphere.cpp -- ISA 1976 国际标准大气模型实现
// ============================================================================

#include "atmosphere.h"
#include "config.h"
#include <cmath>
#include <algorithm>

AtmosphereState calcAtmosphere(double h) {
    // 高度下限
    h = std::max(h, 0.0);

    AtmosphereState atm;

    // 海平面基准值
    constexpr double T0   = 288.15;    // 海平面温度 [K]
    constexpr double p0   = 101325.0;  // 海平面气压 [Pa]
    constexpr double L_tropo = 0.0065; // 对流层温度递减率 [K/m]
    constexpr double h_tropo = 11000.0;// 对流层顶高度 [m]
    constexpr double T_tropo = T0 - L_tropo * h_tropo; // 对流层顶温度 [K]

    if (h <= h_tropo) {
        // ---- 对流层 (0 ~ 11 km) ----
        // T = T0 - L * h
        atm.T = T0 - L_tropo * h;
        // p = p0 * (T/T0)^(g / (R*L))
        double exponent = config::g0 / (config::R_air * L_tropo);
        atm.p = p0 * std::pow(atm.T / T0, exponent);
    } else {
        // ---- 平流层底层 (11 ~ 20 km)，等温层 ----
        // 先算对流层顶气压
        double exponent_tropo = config::g0 / (config::R_air * L_tropo);
        double p_tropo = p0 * std::pow(T_tropo / T0, exponent_tropo);

        atm.T = T_tropo; // 等温
        // p = p_tropo * exp(-g*(h - h_tropo) / (R*T_tropo))
        atm.p = p_tropo * std::exp(-config::g0 * (h - h_tropo) / (config::R_air * T_tropo));
    }

    // 密度: rho = p / (R * T)
    atm.rho = atm.p / (config::R_air * atm.T);

    // 声速: a = sqrt(gamma * R * T)
    atm.a = std::sqrt(config::gamma_air * config::R_air * atm.T);

    return atm;
}
