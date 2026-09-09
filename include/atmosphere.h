#pragma once

// ============================================================================
// atmosphere.h -- ISA 1976 国际标准大气模型
// 输入: 几何高度 h [m]
// 输出: 温度 T [K], 气压 p [Pa], 密度 rho [kg/m^3], 声速 a [m/s]
// 覆盖范围: 0 ~ 20 km (对流层 + 平流层底层)
// ============================================================================

struct AtmosphereState {
    double T;      // 温度 [K]
    double p;      // 气压 [Pa]
    double rho;    // 密度 [kg/m^3]
    double a;      // 声速 [m/s]
};

// 根据几何高度计算标准大气参数
AtmosphereState calcAtmosphere(double h);
