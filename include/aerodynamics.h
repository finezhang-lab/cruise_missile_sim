#pragma once

// ============================================================================
// aerodynamics.h -- 气动力计算模块接口声明
//
// 本文件只声明函数接口(函数名、参数、返回值)，具体实现见 src/aerodynamics.cpp。
// 其他文件(main.cpp、trajectory.cpp)只需 #include "aerodynamics.h" 即可调用。
// ============================================================================

#include <string>
#include "atmosphere.h"   // AtmosphereState (writeAeroPolarCSV 的参数用到)

// 气动力计算结果
struct AeroResult {
    double CL;   // 升力系数 [-]
    double CD;   // 阻力系数 [-]
    double L;    // 升力 [N]
    double D;    // 阻力 [N]
    double LD;   // 升阻比 [-]
};

// 全弹升力线斜率 CL_alpha [1/rad] (基于机翼面积)
double calcCLalpha(double Ma);

// 诱导阻力因子 K (CDi = K·CL²)
double calcInducedDragFactor(double Ma);

// 零升阻力系数 CD0 [-]
double calcCD0(double Ma);

// 计算气动力: alpha [rad], Ma 马赫数, q 动压 [Pa]
AeroResult calcAerodynamics(double alpha, double Ma, double q);

// 压力中心系数 x_cp/L [-]
double calcPressureCenterCoeff(double Ma);

// 质心系数 x_cg/L [-]
double calcCGLocation();

// 俯仰力矩系数 Cm [-]
double calcPitchMomentCoeff(double alpha, double Ma);

// 配平攻角 [rad]
double calcTrimAlpha(double Ma, double q, double weight);

// 生成 CL/CD 随攻角变化的扫描数据 CSV
void writeAeroPolarCSV(const std::string& filename, double Ma,
                       const AtmosphereState& atm,
                       double alphaMinDeg, double alphaMaxDeg, double alphaStepDeg);
