#pragma once

// ============================================================================
// aerodynamics.h -- 气动力计算模块
// 包含: 升力系数、阻力系数、升阻比、静稳定性分析
// 参考面积: 机翼面积 S_WING
// ============================================================================

struct AeroResult {
    double CL;     // 升力系数
    double CD;     // 阻力系数
    double L;      // 升力 [N]
    double D;      // 阻力 [N]
    double LD;     // 升阻比 L/D
};

// 计算给定攻角和马赫数下的气动特性
// alpha: 攻角 [rad]
// Ma:    马赫数
// q:     动压 [Pa], q = 0.5 * rho * V^2
AeroResult calcAerodynamics(double alpha, double Ma, double q);

// 计算升力线斜率 CL_alpha [1/rad]
double calcCLalpha(double Ma);

// 计算零升阻力系数 CD0
double calcCD0(double Ma);

// 计算诱导阻力因子 K
double calcInducedDragFactor(double Ma);

// 静稳定性分析: 计算压力中心系数 x_cp/L (无量纲)
double calcPressureCenterCoeff(double Ma);

// 计算俯仰力矩系数 Cm (相对于质心)
double calcPitchMomentCoeff(double alpha, double Ma);

// 给定马赫数和动压，反算配平攻角 (L = W 时的攻角)
// weight: 导弹重力 [N]
double calcTrimAlpha(double Ma, double q, double weight);
