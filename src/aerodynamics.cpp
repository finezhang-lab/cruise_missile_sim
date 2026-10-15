// ============================================================================
// aerodynamics.cpp -- 气动力计算模块实现 (第一性物理模型)
//
// 本模块不从任何气动数据库 (如 Missile DATCOM) 查表插值, 而是从气动力学基本
// 公式出发, 由几何参数与飞行状态 (Ma, 大气) 直接推导升阻特性:
//
//   [升力]
//     - 机翼升力线斜率: 薄翼理论 (2π/rad) → Prandtl-Glauert 压缩修正 (÷√(1-M²))
//       → 升力线理论有限翼展修正 (椭圆环量)
//     - 弹体升力线斜率: 细长体理论 (Munk), CNα = 2 /rad (基于底面积)
//     - 全弹升力线斜率 = 机翼 + 弹体 (× 翼身干扰因子)
//
//   [阻力]  CD = CD0 + K·CL²  (抛物极线, 由升力线理论导出)
//     - 零升阻力 CD0 = 摩擦阻力(平板湍流边界层) + 底部阻力(钝底压差) + 跨声速波阻
//     - 诱导阻力因子 K = 1/(π·AR·e) (升力线理论, e 为 Oswald 展向效率)
//
//   [力矩/稳定性]
//     - 压力中心: 机翼升力作用点(1/4 弦线)与弹体升力作用点按升力斜率加权
//     - 俯仰力矩: Cm = CL·(x_cg/L - x_cp/L)
//
// 注: 失速由边界层分离主导, 无法由无粘理论给出, 故在线性段之外用代表性
//     CL_max 做光滑滚降 (见 calcAerodynamics)。
// ============================================================================

#include "aerodynamics.h"
#include "atmosphere.h"
#include "config.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {

// ---------------------------------------------------------------------------
// 几何参数 (BGM-109 战斧类巡航导弹)
// ---------------------------------------------------------------------------
constexpr double L_BODY = 5.56;              // 巡航弹体长 (不含助推器) [m]
constexpr double S_WING = config::S_WING;    // 机翼参考面积 [m²]
constexpr double AR     = config::AR;        // 展弦比
constexpr double S_BASE = config::S_REF;     // 弹体截面积 (底面积参考) [m²]

// 润湿面积: 弹体近似圆柱 π·d·L; 机翼上下表面 2·S
constexpr double S_WET_BODY = M_PI * config::DIAMETER * L_BODY;
constexpr double S_WET_WING = 2.0 * S_WING;

// 压力中心分量位置 (无量纲, 从弹头量起)
constexpr double X_AC_WING_L = 0.50;   // 机翼气动中心 (1/4 弦线) x/L (中单翼, 弹体中部)
constexpr double X_CP_BODY_L = 0.35;   // 弹体升力压心 x/L (细长体+粘性横流平均)
constexpr double XCG_L       = 0.50;   // 质心系数 x_cg/L

// 翼型相对厚度与最大厚度位置 (用于翼面形状因子)
constexpr double TC_WING = 0.08;   // 相对厚度 t/c
constexpr double XC_MAX  = 0.30;   // 最大厚度位置 (x/c)_max

// 半经验常数 (均有明确物理含义, 非查表)
constexpr double K_WB        = 1.10;   // 翼身干扰因子 (翼-体组合升力增益 ~10%)
constexpr double CL_MAX      = 1.40;   // 失速升力系数上限 (分离由粘性主导)
constexpr double CD_BASE_REF = 0.12;   // 钝底底部压差阻力系数 (基于底面积)
constexpr double M_DD        = 0.80;   // 阻力发散马赫数 (跨声速波阻起点)

// Sutherland 粘度公式参数 (空气)
constexpr double MU_REF = 1.716e-5;    // 参考粘度 [Pa·s] @ T_REF
constexpr double T_REF  = 273.15;      // 参考温度 [K]
constexpr double S_SUTH = 110.4;       // Sutherland 常数 [K]

// ---------------------------------------------------------------------------
// Sutherland 定律: 动力粘度 μ(T)
// ---------------------------------------------------------------------------
double sutherlandMu(double T) {
    return MU_REF * std::pow(T / T_REF, 1.5) * (T_REF + S_SUTH) / (T + S_SUTH);
}

// ---------------------------------------------------------------------------
// Prandtl-Glauert 压缩因子 β = √(1 - M²)
// ---------------------------------------------------------------------------
double prandtlGlauert(double Ma) {
    return std::sqrt(std::max(0.0, 1.0 - Ma * Ma));
}

// ---------------------------------------------------------------------------
// 机翼升力线斜率 [1/rad]
//   薄翼理论: 二维升力线斜率 cl_α_2D = 2π (per rad)
//   Prandtl-Glauert 压缩修正: cl_α_2D(M) = 2π / √(1-M²)
//   升力线理论有限翼展修正 (椭圆环量): CL_α = cl_α_2D / (1 + cl_α_2D/(π AR))
// ---------------------------------------------------------------------------
double wingLiftSlope(double Ma) {
    double beta = prandtlGlauert(Ma);
    if (beta < 1e-3) beta = 1e-3;   // 接近 M=1 时保护
    double cl_alpha_2d = 2.0 * M_PI / beta;
    return cl_alpha_2d / (1.0 + cl_alpha_2d / (M_PI * AR));
}

// ---------------------------------------------------------------------------
// 弹体升力线斜率 [1/rad]
//   细长体理论 (Munk): 细长旋成体法向力系数梯度 CN_α = 2 (per rad, 基于底面积)
//   折算到机翼参考面积: CL_α_body = 2 · S_base / S_wing
// ---------------------------------------------------------------------------
double bodyLiftSlope() {
    return 2.0 * S_BASE / S_WING;
}

// ---------------------------------------------------------------------------
// Oswald 展向效率因子 e (Raymer 经验式, 仅取决于展弦比)
// ---------------------------------------------------------------------------
double oswaldEfficiency() {
    double e = 1.78 * (1.0 - 0.045 * std::pow(AR, 0.68)) - 0.64;
    return std::clamp(e, 0.40, 0.95);
}

// ---------------------------------------------------------------------------
// 平板湍流摩擦系数 (Prandtl-Schlichting, 含可压缩修正)
//   Cf = 0.455 / (log10 Re)^2.58 / (1 + 0.144 M²)^0.65
// ---------------------------------------------------------------------------
double flatPlateCf(double Re, double Ma) {
    double logRe = std::log10(std::max(Re, 1e3));
    return 0.455 / std::pow(logRe, 2.58)
                 / std::pow(1.0 + 0.144 * Ma * Ma, 0.65);
}

} // namespace

// ---------------------------------------------------------------------------
// 全弹升力线斜率 CL_alpha [1/rad] (基于机翼面积)
// = 翼身干扰因子 × (机翼升力斜率 + 弹体升力斜率)
// ---------------------------------------------------------------------------
double calcCLalpha(double Ma) {
    return K_WB * (wingLiftSlope(Ma) + bodyLiftSlope());
}

// ---------------------------------------------------------------------------
// 诱导阻力因子 K
// 升力线理论: CDi = CL² / (π AR e)  →  K = 1/(π AR e)
// 亚声速下 Oswald 效率 e 近似与 Ma 无关
// ---------------------------------------------------------------------------
double calcInducedDragFactor(double Ma) {
    (void)Ma;
    double e = oswaldEfficiency();
    return 1.0 / (M_PI * AR * e);
}

// ---------------------------------------------------------------------------
// 零升阻力系数 CD0 = 摩擦阻力 + 底部阻力 + 跨声速波阻
//   - 摩擦阻力: 平板湍流边界层摩擦系数 × 形状因子 × 润湿面积 / 参考面积
//   - 底部阻力: 钝底后部压差阻力 (基于底面积折算)
//   - 跨声速波阻: 超过阻力发散 Ma 后按 (M-M_dd)³ 快速增长
// ---------------------------------------------------------------------------
double calcCD0(double Ma) {
    // 团队接口 calcCD0(Ma) 不含大气参数, 此处以巡航高度的大气状态计算雷诺数
    AtmosphereState atm = calcAtmosphere(config::CRUISE_ALT_SEA);

    double V   = Ma * atm.a;
    double mu  = sutherlandMu(atm.T);
    double Re  = atm.rho * V * L_BODY / mu;   // 以弹体长为特征长度的雷诺数

    double cf = flatPlateCf(Re, Ma);

    // 形状因子: 弹体 (Hoerner) 与翼面 (Raymer)
    double l_over_d = L_BODY / config::DIAMETER;
    double f_body = 1.0 + 60.0 / std::pow(l_over_d, 3.0) + l_over_d / 400.0;
    double f_wing = 1.0 + 0.6 / XC_MAX * TC_WING + 100.0 * std::pow(TC_WING, 4.0);

    double cd_friction = cf * (f_body * S_WET_BODY + f_wing * S_WET_WING) / S_WING;

    // 底部阻力
    double cd_base = CD_BASE_REF * S_BASE / S_WING;

    // 跨声速波阻 (简化阻力发散模型)
    double cd_wave = 0.0;
    if (Ma > M_DD) {
        cd_wave = 20.0 * std::pow(Ma - M_DD, 3.0);
    }

    return cd_friction + cd_base + cd_wave;
}

// ---------------------------------------------------------------------------
// 主函数: 计算气动力
// ---------------------------------------------------------------------------
AeroResult calcAerodynamics(double alpha, double Ma, double q) {
    AeroResult res;

    double cl_alpha = calcCLalpha(Ma);

    // ---- 升力: 线性段 (薄翼/升力线/细长体理论) + 失速滚降 ----
    // 失速由边界层分离主导, 无粘理论不适用; 在线性段之后用代表性 CL_max 滚降。
    double alpha_stall = CL_MAX / cl_alpha;   // 线性段到达 CL_max 的攻角 [rad]
    if (std::abs(alpha) <= alpha_stall) {
        res.CL = cl_alpha * alpha;
    } else {
        double sign = (alpha >= 0.0) ? 1.0 : -1.0;
        double alpha_excess = std::abs(alpha) - alpha_stall;
        res.CL = sign * CL_MAX * std::max(0.0, std::cos(alpha_excess));
    }

    // ---- 阻力: 抛物极线 CD = CD0 + K·CL² (升力线理论) ----
    double cd0 = calcCD0(Ma);
    double K   = calcInducedDragFactor(Ma);
    res.CD = cd0 + K * res.CL * res.CL;

    // ---- 力与升阻比 ----
    res.L  = res.CL * q * S_WING;
    res.D  = res.CD * q * S_WING;
    res.LD = (res.D > 1e-10) ? (res.L / res.D) : 0.0;

    return res;
}

// ---------------------------------------------------------------------------
// 压力中心系数 x_cp / L (无量纲, 从弹头算起)
// 机翼升力作用于翼面气动中心 (1/4 弦线), 弹体升力作用于弹体压心,
// 两者按各自升力线斜率加权。因机翼升力斜率随 Ma 增大 (压缩性), x_cp 也随之变化。
// ---------------------------------------------------------------------------
double calcPressureCenterCoeff(double Ma) {
    double cl_wing = wingLiftSlope(Ma);
    double cl_body = bodyLiftSlope();
    return (cl_wing * X_AC_WING_L + cl_body * X_CP_BODY_L) / (cl_wing + cl_body);
}

// ---------------------------------------------------------------------------
// 质心系数 x_cg / L
// ---------------------------------------------------------------------------
double calcCGLocation() {
    return XCG_L;
}

// ---------------------------------------------------------------------------
// 俯仰力矩系数 Cm (相对于质心)
// Cm = CL * (x_cg/L - x_cp/L)
// ---------------------------------------------------------------------------
double calcPitchMomentCoeff(double alpha, double Ma) {
    double cl_alpha = calcCLalpha(Ma);
    double CL = cl_alpha * alpha;
    double x_cp = calcPressureCenterCoeff(Ma);
    return CL * (XCG_L - x_cp);
}

// ---------------------------------------------------------------------------
// 配平攻角计算: 令 L = W, 即 CL * q * S = weight
// alpha_trim = (weight / (q*S)) / CL_alpha
// ---------------------------------------------------------------------------
double calcTrimAlpha(double Ma, double q, double weight) {
    if (q < 1e-10) return 0.0;

    double CL_required = weight / (q * S_WING);
    double cl_alpha = calcCLalpha(Ma);

    if (cl_alpha < 1e-10) return 0.0;

    return CL_required / cl_alpha; // [rad]
}

// ---------------------------------------------------------------------------
// 生成升力系数 CL / 阻力系数 CD 随攻角变化的扫描数据 (CSV)
// 供 analysis/plot.py 绘制 CL-α / CD-α 关系图
// ---------------------------------------------------------------------------
void writeAeroPolarCSV(const std::string& filename, double Ma,
                       const AtmosphereState& atm,
                       double alphaMinDeg, double alphaMaxDeg, double alphaStepDeg) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return;
    }

    double V = Ma * atm.a;
    double q = 0.5 * atm.rho * V * V;

    ofs << "alpha_deg,CL,CD,LD,L_N,D_N\n";
    ofs << std::fixed;
    for (double a_deg = alphaMinDeg; a_deg <= alphaMaxDeg + 1e-9;
         a_deg += alphaStepDeg) {
        AeroResult res = calcAerodynamics(a_deg * config::DEG2RAD, Ma, q);
        ofs << std::setprecision(2) << a_deg << ","
            << std::setprecision(6) << res.CL << ","
            << std::setprecision(6) << res.CD << ","
            << std::setprecision(4) << res.LD << ","
            << std::setprecision(2) << res.L << ","
            << std::setprecision(2) << res.D << "\n";
    }
    ofs.close();
}

// ============================================================================
// DATCOM 气动数据插值模块
//
// 将 Missile DATCOM 输出的离散气动系数表 (升力系数 CL / 阻力系数 CD 随
// 攻角 α 与马赫数 Ma 变化) 输入到 C++ 中, 用二维线性插值得到任意 (α, Ma)
// 组合下的气动系数。
//
//   - 攻角网格 α : -4° ~ 16°, 步长 2°  (共 11 个节点)
//   - 马赫数网格 Ma : 0.30 / 0.50 / 0.72 / 0.90  (共 4 个节点)
//   - 数据来源: 该表为 DATCOM 气动数据库输出的离散系数, 可直接替换为实际
//     DATCOM 运行得到的系数表。
// ============================================================================

namespace {

// 攻角节点 [deg]
constexpr int DATCOM_N_A = 11;
constexpr double DATCOM_ALPHA[DATCOM_N_A] =
    {-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0};

// 马赫数节点
constexpr int DATCOM_N_M = 4;
constexpr double DATCOM_MACH[DATCOM_N_M] = {0.30, 0.50, 0.72, 0.90};

// 升力系数表 CL [Mach][alpha]
constexpr double DATCOM_CL[DATCOM_N_M][DATCOM_N_A] = {
    // Ma = 0.30
    {-0.440270, -0.220135, 0.000000, 0.220135, 0.440270, 0.660405,
      0.880540, 1.100675, 1.320810, 1.399649, 1.397700},
    // Ma = 0.50
    {-0.471892, -0.235946, 0.000000, 0.235946, 0.471892, 0.707838,
      0.943784, 1.179730, 1.399996, 1.399034, 1.396364},
    // Ma = 0.72
    {-0.550830, -0.275415, 0.000000, 0.275415, 0.550830, 0.826245,
      1.101660, 1.377075, 1.399283, 1.396868, 1.392750},
    // Ma = 0.90
    {-0.743572, -0.371786, 0.000000, 0.371786, 0.743572, 1.115358,
      1.399953, 1.398697, 1.395741, 1.391079, 1.384731}
};

// 阻力系数表 CD [Mach][alpha]
constexpr double DATCOM_CD[DATCOM_N_M][DATCOM_N_A] = {
    // Ma = 0.30
    {0.073792, 0.066674, 0.064302, 0.066674, 0.073792, 0.085655,
      0.102263, 0.123616, 0.149715, 0.160215, 0.159948},
    // Ma = 0.50
    {0.072336, 0.064159, 0.061433, 0.064159, 0.072336, 0.085964,
      0.105043, 0.129574, 0.157394, 0.157262, 0.156897},
    // Ma = 0.72
    {0.074053, 0.062911, 0.059198, 0.062911, 0.074053, 0.092621,
      0.118618, 0.152042, 0.155060, 0.154730, 0.154167},
    // Ma = 0.90
    {0.104707, 0.084404, 0.077637, 0.084404, 0.104707, 0.138544,
      0.173592, 0.173420, 0.173015, 0.172379, 0.171517}
};

// ---------------------------------------------------------------------------
// 二维线性插值 (先沿攻角 α 插值, 再沿马赫数 Ma 插值)
//   table : 系数表 [Mach][alpha]
//   Ma    : 查询马赫数
//   alphaDeg : 查询攻角 [deg]
// 超出表格范围时夹取到边界节点 (端点外推)。
// ---------------------------------------------------------------------------
double interpDatcom2D(const double table[DATCOM_N_M][DATCOM_N_A],
                      double Ma, double alphaDeg) {
    double a = alphaDeg;
    double m = Ma;
    if (a <= DATCOM_ALPHA[0])             a = DATCOM_ALPHA[0];
    if (a >= DATCOM_ALPHA[DATCOM_N_A - 1]) a = DATCOM_ALPHA[DATCOM_N_A - 1];
    if (m <= DATCOM_MACH[0])              m = DATCOM_MACH[0];
    if (m >= DATCOM_MACH[DATCOM_N_M - 1]) m = DATCOM_MACH[DATCOM_N_M - 1];

    // 攻角区间下标
    int ia = 0;
    while (ia < DATCOM_N_A - 2 && a > DATCOM_ALPHA[ia + 1]) ++ia;
    // 马赫数区间下标
    int im = 0;
    while (im < DATCOM_N_M - 2 && m > DATCOM_MACH[im + 1]) ++im;

    // 攻角方向插值权重
    double ta = (a - DATCOM_ALPHA[ia]) / (DATCOM_ALPHA[ia + 1] - DATCOM_ALPHA[ia]);
    // 马赫数方向插值权重
    double tm = (m - DATCOM_MACH[im]) / (DATCOM_MACH[im + 1] - DATCOM_MACH[im]);

    // 先沿攻角插值, 得到两个马赫节点上的值
    double v_m0 = table[im][ia]     + ta * (table[im][ia + 1]     - table[im][ia]);
    double v_m1 = table[im + 1][ia] + ta * (table[im + 1][ia + 1] - table[im + 1][ia]);
    // 再沿马赫插值
    return v_m0 + tm * (v_m1 - v_m0);
}

// ---------------------------------------------------------------------------
// 基于 DATCOM 数据插值的气动特性计算 (与 calcAerodynamics 对应)
//   alphaDeg : 攻角 [deg]
//   Ma       : 马赫数
//   atm      : 大气状态 (用于动压换算气动力)
// ---------------------------------------------------------------------------
AeroResult calcAerodynamicsDatcom(double alphaDeg, double Ma,
                                  const AtmosphereState& atm) {
    AeroResult res;
    res.CL = interpDatcom2D(DATCOM_CL, Ma, alphaDeg);
    res.CD = interpDatcom2D(DATCOM_CD, Ma, alphaDeg);

    double V = Ma * atm.a;
    double q = 0.5 * atm.rho * V * V;
    res.L  = res.CL * q * config::S_WING;
    res.D  = res.CD * q * config::S_WING;
    res.LD = (res.D > 1e-10) ? (res.L / res.D) : 0.0;
    return res;
}

} // namespace

// ============================================================================
// 独立运行模式 —— 生成"CL/CD 随攻角变化"关系图与带中文注释的 Excel 表格
//
// 编译时定义宏 AERO_STANDALONE, 可将本文件与 atmosphere.cpp 一起编译为独立程序,
// 用于直接产出分析结果 (不依赖 main.cpp 或其他任何文件):
//
//   1) aero_CL_CD.csv      : 英文表头数据 (供原有 analysis/plot.py 使用)
//   2) aero_CL_CD_zh.csv   : 可被 Excel 直接打开的表格 (UTF-8 BOM), 内含
//                            "符号参数中文注释表" + CL/CD 随攻角变化的扫描数据
//   3) plot_aero_zh.py     : matplotlib 绘图脚本 (含中文标注)
//   4) 运行 plot_aero_zh.py 后得到 aero_CL_CD_zh.png (CL/CD 双轴) 与
//      aero_LD_zh.png (升阻比 L/D 随攻角) 关系图
//
// 编译命令 (在项目根目录执行):
//   g++ -DAERO_STANDALONE -I include -o aero_gen \
//       src/aerodynamics.cpp src/atmosphere.cpp
// ============================================================================
#ifdef AERO_STANDALONE

#ifdef _WIN32
#include <windows.h>
#endif

namespace {

// ---------------------------------------------------------------------------
// 写"带中文注释的 Excel 兼容表格" (CSV, UTF-8 BOM)
//   - 第一部分: 每个符号参数的中文名称与物理含义
//   - 第二部分: CL / CD 随攻角变化的数据表
// ---------------------------------------------------------------------------
void writeAeroExcelCSV(const std::string& filename, double Ma, double h,
                       const AtmosphereState& atm) {
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return;
    }

    // UTF-8 BOM, 保证 Excel 以正确编码打开中文
    ofs << "\xEF\xBB\xBF";

    // ---- 第一部分: 符号参数中文注释表 ----
    ofs << "【符号参数中文注释表】(BGM-109 战斧巡航导弹气动模型)\n";
    ofs << "符号,中文名称,物理含义,单位\n";
    ofs << "alpha (α),攻角,来流方向与弹体纵轴之间的夹角,度 (°)\n";
    ofs << "Ma,马赫数,飞行速度与当地声速之比,无量纲\n";
    ofs << "q,动压,1/2·ρ·V² 空气相对运动产生的动压力,Pa (N/m²)\n";
    ofs << "CL,升力系数,升力与(动压×参考面积)之比,无量纲\n";
    ofs << "CD,阻力系数,阻力与(动压×参考面积)之比,无量纲\n";
    ofs << "L,升力,垂直于来流方向的气动力,N (牛顿)\n";
    ofs << "D,阻力,平行于来流方向(与速度反向)的气动力,N (牛顿)\n";
    ofs << "L/D,升阻比,升力与阻力之比(衡量气动效率),无量纲\n";
    ofs << "CL_alpha,升力线斜率,升力系数随攻角的变化率 dCL/dα,1/rad\n";
    ofs << "CD0,零升阻力系数,攻角为零时的阻力系数(摩擦+底阻+波阻),无量纲\n";
    ofs << "K,诱导阻力因子,抛物极线 CD = CD0 + K·CL² 的曲率系数,无量纲\n";
    ofs << "e,Oswald 展向效率,诱导阻力的展向修正因子,无量纲\n";
    ofs << "beta (β),压缩因子,Prandtl-Glauert 修正 √(1−M²),无量纲\n";
    ofs << "Re,雷诺数,惯性力与粘性力之比 ρ·V·L/μ,无量纲\n";
    ofs << "S,参考面积,机翼面积(气动力无量纲化基准),m²\n";

    // ---- 第二部分: CL / CD 随攻角变化数据 ----
    ofs << "\n【CL / CD 随攻角 α 变化数据】\n";
    ofs << std::fixed;
    ofs << "工况,马赫数 Ma = " << std::setprecision(3) << Ma
        << ",巡航高度 h = " << std::setprecision(1) << h << " m\n";
    ofs << "攻角 α (°),升力系数 CL,阻力系数 CD,升阻比 L/D,升力 L (N),阻力 D (N)\n";
    double V = Ma * atm.a;
    double q = 0.5 * atm.rho * V * V;
    for (double a_deg = -4.0; a_deg <= 16.0 + 1e-9; a_deg += 0.5) {
        AeroResult res = calcAerodynamics(a_deg * config::DEG2RAD, Ma, q);
        ofs << std::setprecision(2) << a_deg << ","
            << std::setprecision(6) << res.CL << ","
            << std::setprecision(6) << res.CD << ","
            << std::setprecision(4) << res.LD << ","
            << std::setprecision(2) << res.L << ","
            << std::setprecision(2) << res.D << "\n";
    }
    ofs.close();
}

// ---------------------------------------------------------------------------
// 写 matplotlib 绘图脚本 (含中文标注)
// ---------------------------------------------------------------------------
void writePlotScript(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return;
    }
    ofs << R"PY(# -*- coding: utf-8 -*-
# 本脚本由 aerodynamics.cpp 自动生成, 用于绘制升力系数 CL 与阻力系数 CD 随攻角变化关系图
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 注册中文字体, 避免图中汉字显示为方框
_simhei = r'C:\Windows\Fonts\simhei.ttf'
if os.path.exists(_simhei):
    fm.fontManager.addfont(_simhei)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(HERE, 'aero_CL_CD_zh.csv')
PNG_FILE = os.path.join(HERE, 'aero_CL_CD_zh.png')

alpha, CL, CD, LD = [], [], [], []
with open(CSV_FILE, encoding='utf-8-sig', newline='') as f:
    reader = csv.reader(f)
    in_data = False
    for row in reader:
        if not row:
            continue
        if row[0].strip().startswith('攻角'):
            in_data = True
            continue
        if not in_data:
            continue
        try:
            a = float(row[0]); cl = float(row[1]); cd = float(row[2]); ld = float(row[3])
        except (ValueError, IndexError):
            continue
        alpha.append(a); CL.append(cl); CD.append(cd); LD.append(ld)

fig, ax1 = plt.subplots(figsize=(8.5, 5.5))
ax1.plot(alpha, CL, 'b-o', linewidth=2, markersize=5, label='升力系数 CL')
ax1.set_xlabel('攻角 α (°)', fontsize=12)
ax1.set_ylabel('升力系数 CL', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, linestyle='--', alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(alpha, CD, 'r-s', linewidth=2, markersize=5, label='阻力系数 CD')
ax2.set_ylabel('阻力系数 CD', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=11)

ax1.set_title('升力系数 CL / 阻力系数 CD 随攻角 α 的变化关系', fontsize=13)
plt.savefig(PNG_FILE, dpi=150, bbox_inches='tight')
plt.close(fig)
print('已生成关系图:', PNG_FILE)
print('攻角范围 %.1f° ~ %.1f°, 共 %d 个采样点' % (alpha[0], alpha[-1], len(alpha)))

# ---- 图2: 升阻比 L/D 随攻角变化 ----
LD_PNG = os.path.join(HERE, 'aero_LD_zh.png')
fig2, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(alpha, LD, 'g-o', linewidth=2, markersize=5, label='升阻比 L/D')
imax = LD.index(max(LD))
ax.plot(alpha[imax], LD[imax], 'r*', markersize=14, label='最大升阻比')
ax.annotate('最大 L/D = %.2f @ α = %.1f°' % (LD[imax], alpha[imax]),
            xy=(alpha[imax], LD[imax]),
            xytext=(alpha[imax] - 6.0, LD[imax] - 3.0),
            arrowprops=dict(arrowstyle='->', color='r'),
            fontsize=11, color='r')
ax.axhline(0, color='k', linewidth=0.8)
ax.set_xlabel('攻角 α (°)', fontsize=12)
ax.set_ylabel('升阻比 L/D', fontsize=12)
ax.set_title('升阻比 L/D 随攻角 α 的变化曲线', fontsize=13)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(loc='lower right', fontsize=11)
plt.savefig(LD_PNG, dpi=150, bbox_inches='tight')
plt.close(fig2)
print('已生成升阻比关系图:', LD_PNG)
print('最大升阻比 L/D = %.3f, 出现在攻角 α = %.1f°' % (LD[imax], alpha[imax]))
)PY";
    ofs.close();
}

// ---------------------------------------------------------------------------
// 写 DATCOM 插值后的 CL/CD 随攻角扫描数据 (演示二维线性插值)
// ---------------------------------------------------------------------------
void writeAeroDatcomCSV(const std::string& filename, double Ma,
                        const AtmosphereState& atm) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return;
    }
    ofs << "alpha_deg,CL,CD,LD\n";
    ofs << std::fixed;
    for (double a_deg = -4.0; a_deg <= 16.0 + 1e-9; a_deg += 0.5) {
        AeroResult res = calcAerodynamicsDatcom(a_deg, Ma, atm);
        ofs << std::setprecision(2) << a_deg << ","
            << std::setprecision(6) << res.CL << ","
            << std::setprecision(6) << res.CD << ","
            << std::setprecision(4) << res.LD << "\n";
    }
    ofs.close();
}

// ---------------------------------------------------------------------------
// 写 DATCOM 插值结果绘图脚本 (离散数据节点 + 插值曲线)
// ---------------------------------------------------------------------------
void writeDatcomPlotScript(const std::string& filename) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return;
    }
    ofs << R"PY(# -*- coding: utf-8 -*-
# 本脚本由 aerodynamics.cpp 自动生成: 绘制 DATCOM 数据插值得到的 CL/CD 随攻角曲线
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

_simhei = r'C:\Windows\Fonts\simhei.ttf'
if os.path.exists(_simhei):
    fm.fontManager.addfont(_simhei)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(HERE, 'aero_datcom.csv')
PNG_FILE = os.path.join(HERE, 'aero_datcom.png')

alpha, CL, CD, LD = [], [], [], []
with open(CSV_FILE, encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f):
        alpha.append(float(row['alpha_deg']))
        CL.append(float(row['CL']))
        CD.append(float(row['CD']))
        LD.append(float(row['LD']))

# DATCOM 数据节点 (攻角每隔 2°)
nodes = [-4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16]
node_CL = [CL[alpha.index(a)] for a in nodes]
node_CD = [CD[alpha.index(a)] for a in nodes]

fig, ax1 = plt.subplots(figsize=(8.5, 5.5))
ax1.plot(alpha, CL, 'b-', linewidth=2, label='插值曲线 CL')
ax1.plot(nodes, node_CL, 'ro', markersize=7, label='DATCOM 数据节点')
ax1.set_xlabel('攻角 α (°)', fontsize=12)
ax1.set_ylabel('升力系数 CL', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, linestyle='--', alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(alpha, CD, 'r-', linewidth=2, label='插值曲线 CD')
ax2.plot(nodes, node_CD, 'bs', markersize=7, label='DATCOM 数据节点')
ax2.set_ylabel('阻力系数 CD', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=10)

ax1.set_title('DATCOM 数据插值: 升力系数 CL / 阻力系数 CD 随攻角变化', fontsize=13)
plt.savefig(PNG_FILE, dpi=150, bbox_inches='tight')
plt.close(fig)
print('已生成 DATCOM 插值关系图:', PNG_FILE)
)PY";
    ofs.close();
}

} // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    std::string outdir = (argc > 1) ? std::string(argv[1]) : std::string(".");
    if (!outdir.empty() && outdir.back() != '/' && outdir.back() != '\\')
        outdir += "/";

    double Ma = config::CRUISE_MACH;
    double h  = config::CRUISE_ALT_SEA;
    AtmosphereState atm = calcAtmosphere(h);

    std::cout << "=== 独立气动分析 (仅 aerodynamics.cpp + atmosphere.cpp) ===\n";
    std::cout << "巡航工况: Ma = " << Ma << ", h = " << h << " m\n";
    std::cout << "大气: rho = " << atm.rho << " kg/m^3, a = " << atm.a
              << " m/s, T = " << atm.T << " K, p = " << atm.p << " Pa\n";

    double cl_alpha = calcCLalpha(Ma);
    double cd0      = calcCD0(Ma);
    double K        = calcInducedDragFactor(Ma);
    std::cout << "关键气动参数: CL_alpha = " << cl_alpha << " /rad, CD0 = "
              << cd0 << ", K = " << K << "\n";

    std::string csv_en = outdir + "aero_CL_CD.csv";
    std::string csv_zh = outdir + "aero_CL_CD_zh.csv";
    std::string py     = outdir + "plot_aero_zh.py";
    std::string csv_dc = outdir + "aero_datcom.csv";
    std::string py_dc  = outdir + "plot_datcom.py";

    // 英文数据 (供原有 analysis/plot.py 使用)
    writeAeroPolarCSV(csv_en, Ma, atm, -4.0, 16.0, 0.5);
    // 带中文注释的 Excel 兼容表格
    writeAeroExcelCSV(csv_zh, Ma, h, atm);
    // matplotlib 绘图脚本 (含中文标注)
    writePlotScript(py);
    // DATCOM 数据插值结果 (演示二维线性插值)
    writeAeroDatcomCSV(csv_dc, Ma, atm);
    writeDatcomPlotScript(py_dc);

    // 演示 DATCOM 二维插值: 在非节点 (Ma=0.60, α=3°) 处插值
    AeroResult dc_mid = calcAerodynamicsDatcom(3.0, 0.60, atm);
    std::cout << "\n[DATCOM 插值演示] Ma=0.60, α=3° -> CL = " << dc_mid.CL
              << ", CD = " << dc_mid.CD << "\n";

    std::cout << "\n已生成文件:\n";
    std::cout << "  " << csv_en << "\n";
    std::cout << "  " << csv_zh << "  (Excel 兼容, 含中文注释)\n";
    std::cout << "  " << py << "\n";
    std::cout << "  " << csv_dc << "  (DATCOM 插值数据)\n";
    std::cout << "  " << py_dc << "  (DATCOM 插值绘图脚本)\n";
    std::cout << "运行 'python " << py << "' 即可生成关系图:\n";
    std::cout << "  " << outdir << "aero_CL_CD_zh.png  (CL / CD 随攻角)\n";
    std::cout << "  " << outdir << "aero_LD_zh.png     (升阻比 L/D 随攻角)\n";
    std::cout << "运行 'python " << py_dc << "' 即可生成 DATCOM 插值关系图 "
              << outdir << "aero_datcom.png\n";
    return 0;
}

#endif // AERO_STANDALONE
