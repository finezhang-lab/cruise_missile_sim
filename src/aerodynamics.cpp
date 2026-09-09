// ============================================================================
// aerodynamics.cpp -- 气动力计算模块实现
// 基于 BGM-109 "战斧"巡航导弹外形进行估算
// 气动参数采用工程估算方法，参考有翼导弹飞行动力学
// ============================================================================

#include "aerodynamics.h"
#include "config.h"
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// 升力线斜率 CL_alpha [1/rad]
// 采用 Helmbold-Diederich 公式估算有限翼展升力线斜率:
//   CL_alpha = (2*pi*AR) / (2 + sqrt(4 + AR^2 * beta^2))
// 其中 beta = sqrt(1 - Ma^2)（亚声速）
// ---------------------------------------------------------------------------
double calcCLalpha(double Ma) {
    double beta;
    if (Ma < 1.0) {
        beta = std::sqrt(1.0 - Ma * Ma);
    } else {
        beta = std::sqrt(Ma * Ma - 1.0); // 超声速近似
    }

    double AR = config::AR;
    // Helmbold-Diederich 公式
    double cl_alpha = (2.0 * M_PI * AR) / (2.0 + std::sqrt(4.0 + AR * AR * beta * beta));

    // 弹体贡献修正（弹体升力约占全弹升力 10%~15%）
    cl_alpha *= 1.12;

    return cl_alpha; // [1/rad]
}

// ---------------------------------------------------------------------------
// 零升阻力系数 CD0
// 基于等效平板摩擦阻力 + 压差阻力 + 干扰阻力估算
// 亚声速: CD0 约 0.015~0.025 (典型巡航导弹量级)
// 跨声速区有阻力发散效应
// ---------------------------------------------------------------------------
double calcCD0(double Ma) {
    // 基准零升阻力系数 (亚声速巡航状态)
    double cd0_base = 0.018;

    if (Ma < 0.8) {
        // 亚声速段: 缓慢变化
        return cd0_base * (1.0 + 0.1 * Ma * Ma);
    } else if (Ma < 1.0) {
        // 跨声速段: 阻力发散
        double ratio = (Ma - 0.8) / 0.2; // 0 ~ 1
        return cd0_base * (1.0 + 0.1 * Ma * Ma + 0.8 * ratio * ratio);
    } else {
        // 超声速段
        return cd0_base * (1.0 + 0.1 + 0.8) / Ma; // 随 Ma 下降
    }
}

// ---------------------------------------------------------------------------
// 诱导阻力因子 K
// CD_induced = K * CL^2
// K = 1 / (pi * AR * e), e 为 Oswald 效率因子 (典型值 0.7~0.85)
// ---------------------------------------------------------------------------
double calcInducedDragFactor(double Ma) {
    double e = 0.78; // Oswald 效率因子（巡航导弹典型值）
    double AR = config::AR;

    double K = 1.0 / (M_PI * AR * e);

    // 跨声速修正: K 在跨声速区略有增大
    if (Ma > 0.8 && Ma < 1.0) {
        K *= 1.0 + 0.3 * (Ma - 0.8) / 0.2;
    }

    return K;
}

// ---------------------------------------------------------------------------
// 主函数: 计算气动力
// ---------------------------------------------------------------------------
AeroResult calcAerodynamics(double alpha, double Ma, double q) {
    AeroResult res;

    // 升力系数: CL = CL_alpha * alpha (线性段)
    // 失速模型: alpha > alpha_stall 时 CL 下降
    double cl_alpha = calcCLalpha(Ma);
    constexpr double alpha_stall = 14.0 * M_PI / 180.0; // 失速攻角 [rad]

    if (std::abs(alpha) <= alpha_stall) {
        res.CL = cl_alpha * alpha;
    } else {
        // 失速后近似: CL 下降
        double sign = (alpha >= 0) ? 1.0 : -1.0;
        double alpha_excess = std::abs(alpha) - alpha_stall;
        double CL_max = cl_alpha * alpha_stall;
        res.CL = sign * CL_max * std::cos(alpha_excess);
    }

    // 阻力系数: CD = CD0 + K * CL^2
    double cd0 = calcCD0(Ma);
    double K = calcInducedDragFactor(Ma);
    res.CD = cd0 + K * res.CL * res.CL;

    // 升力和阻力 (参考面积: 机翼面积)
    double S = config::S_WING;
    res.L = res.CL * q * S;
    res.D = res.CD * q * S;

    // 升阻比
    if (res.D > 1e-10) {
        res.LD = res.L / res.D;
    } else {
        res.LD = 0.0;
    }

    return res;
}

// ---------------------------------------------------------------------------
// 压力中心系数 x_cp / L (无量纲, 从弹头算起)
// 亚声速: 约 45%~55% 弹长
// 跨声速后移: 约 50%~60%
// ---------------------------------------------------------------------------
double calcPressureCenterCoeff(double Ma) {
    if (Ma < 0.8) {
        return 0.48; // 亚声速
    } else if (Ma < 1.0) {
        double ratio = (Ma - 0.8) / 0.2;
        return 0.48 + 0.10 * ratio; // 跨声速后移
    } else {
        return 0.58; // 超声速
    }
}

// ---------------------------------------------------------------------------
// 俯仰力矩系数 Cm (相对于质心)
// Cm = CL * (x_cg/L - x_cp/L)
// 质心位置假设在 45% 弹长处
// Cm < 0 表示静稳定
// ---------------------------------------------------------------------------
double calcPitchMomentCoeff(double alpha, double Ma) {
    double cl_alpha = calcCLalpha(Ma);
    double CL = cl_alpha * alpha;

    constexpr double x_cg_ratio = 0.45; // 质心位置 (从弹头)
    double x_cp_ratio = calcPressureCenterCoeff(Ma);

    // Cm = CL * (x_cg - x_cp) / L
    // 当 x_cp > x_cg (压力中心在质心之后) 时 Cm < 0 => 静稳定
    double Cm = CL * (x_cg_ratio - x_cp_ratio);
    return Cm;
}

// ---------------------------------------------------------------------------
// 配平攻角计算: 令 L = W, 即 CL * q * S = weight
// CL = weight / (q * S)
// alpha_trim = CL / CL_alpha
// ---------------------------------------------------------------------------
double calcTrimAlpha(double Ma, double q, double weight) {
    if (q < 1e-10) return 0.0;

    double S = config::S_WING;
    double CL_required = weight / (q * S);
    double cl_alpha = calcCLalpha(Ma);

    if (cl_alpha < 1e-10) return 0.0;

    double alpha_trim = CL_required / cl_alpha;
    return alpha_trim; // [rad]
}
