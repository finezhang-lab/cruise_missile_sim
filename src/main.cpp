// ============================================================================
// main.cpp -- BGM-109 "战斧"巡航导弹质点弹道仿真主程序
//
// 功能:
//   1. 气动特性分析 (阶段一)
//   2. 动力特性分析 (阶段二)
//   3. 基准弹道计算 (阶段三/四)
//   4. 弹道特性参数影响分析 (变参数扫描)
//
// 使用方法:
//   ./cruise_missile_sim              # 运行全部分析
//   ./cruise_missile_sim aero         # 仅运行气动分析
//   ./cruise_missile_sim propulsion   # 仅运行动力分析
//   ./cruise_missile_sim trajectory   # 仅运行弹道计算
//   ./cruise_missile_sim sweep        # 仅运行弹道特性分析 (§3-三)
//   ./cruise_missile_sim baseline      # 仅运行基准工况分析
// ============================================================================

#include "config.h"
#include "atmosphere.h"
#include "aerodynamics.h"
#include "propulsion.h"
#include "trajectory.h"
#include "integrator.h"
#include "data_output.h"
#include "parametric_study.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>

// 跨平台目录创建
#ifdef _WIN32
#include <direct.h>
static void ensureDirectory(const std::string& path) { _mkdir(path.c_str()); }
#else
#include <sys/stat.h>
static void ensureDirectory(const std::string& path) { mkdir(path.c_str(), 0755); }
#endif

// ---------------------------------------------------------------------------
// 阶段一: 气动特性分析
// ---------------------------------------------------------------------------
static void runAeroAnalysis() {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "  阶段一: 升阻特性分析\n";
    std::cout << "============================================================\n";

    std::string outdir = "output";
    ensureDirectory(outdir);

    // --- 1.1 升力系数和阻力系数随攻角变化 ---
    std::cout << "\n--- 1.1 CL, CD 随攻角变化 (Ma=0.72, H=15m) ---\n";
    std::cout << std::setw(10) << "alpha_deg"
              << std::setw(10) << "CL"
              << std::setw(10) << "CD"
              << std::setw(10) << "L/D" << "\n";

    AtmosphereState atm_cruise = calcAtmosphere(config::CRUISE_ALT_SEA);
    double V_cruise = config::CRUISE_MACH * atm_cruise.a;
    double q_cruise = 0.5 * atm_cruise.rho * V_cruise * V_cruise;

    std::ofstream ofs_clcd(outdir + "/aero_CL_CD.csv");
    ofs_clcd << "alpha_deg,CL,CD,LD,L_N,D_N\n";

    for (double alpha_deg = -4.0; alpha_deg <= 16.0; alpha_deg += 0.5) {
        double alpha = alpha_deg * config::DEG2RAD;
        AeroResult res = calcAerodynamics(alpha, config::CRUISE_MACH, q_cruise);

        if (std::abs(alpha_deg - std::round(alpha_deg * 2) / 2.0) < 0.01
            && std::fmod(alpha_deg, 2.0) < 0.01) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(1) << alpha_deg
                      << std::setw(10) << std::setprecision(4) << res.CL
                      << std::setw(10) << std::setprecision(4) << res.CD
                      << std::setw(10) << std::setprecision(2) << res.LD << "\n";
        }
        ofs_clcd << std::fixed << std::setprecision(2) << alpha_deg << ","
                 << std::setprecision(6) << res.CL << ","
                 << std::setprecision(6) << res.CD << ","
                 << std::setprecision(4) << res.LD << ","
                 << std::setprecision(2) << res.L << ","
                 << std::setprecision(2) << res.D << "\n";
    }
    ofs_clcd.close();
    std::cout << "  数据已保存至: " << outdir << "/aero_CL_CD.csv\n";

    // --- 1.2 不同马赫数下的升阻特性 ---
    std::cout << "\n--- 1.2 不同马赫数下的 CL_alpha 和 CD0 ---\n";
    std::cout << std::setw(8) << "Ma"
              << std::setw(12) << "CL_alpha"
              << std::setw(10) << "CD0"
              << std::setw(10) << "K" << "\n";

    std::ofstream ofs_mach(outdir + "/aero_mach_sweep.csv");
    ofs_mach << "Ma,CL_alpha,CD0,K\n";
    double mach_values[] = {0.3, 0.5, 0.6, 0.65, 0.72, 0.75, 0.8, 0.85, 0.9, 0.95};
    for (double Ma : mach_values) {
        double cla = calcCLalpha(Ma);
        double cd0 = calcCD0(Ma);
        double K = calcInducedDragFactor(Ma);
        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << Ma
                  << std::setw(12) << std::setprecision(4) << cla
                  << std::setw(10) << std::setprecision(5) << cd0
                  << std::setw(10) << std::setprecision(4) << K << "\n";
        ofs_mach << Ma << "," << cla << "," << cd0 << "," << K << "\n";
    }
    ofs_mach.close();

    // --- 1.3 静稳定性分析 ---
    std::cout << "\n--- 1.3 静稳定性分析 ---\n";
    std::cout << std::setw(8) << "Ma"
              << std::setw(12) << "x_cp/L"
              << std::setw(12) << "x_cg/L"
              << std::setw(10) << "Cm(3deg)" << "\n";

    for (double Ma : mach_values) {
        double x_cp = calcPressureCenterCoeff(Ma);
        double Cm = calcPitchMomentCoeff(3.0 * config::DEG2RAD, Ma);
        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << Ma
                  << std::setw(12) << std::setprecision(4) << x_cp
                  << std::setw(12) << std::setprecision(4) << 0.45
                  << std::setw(10) << std::setprecision(5) << Cm;
        if (Cm < 0) std::cout << " (稳定)";
        else std::cout << " (不稳定)";
        std::cout << "\n";
    }

    // --- 1.4 配平攻角 ---
    std::cout << "\n--- 1.4 巡航配平攻角估算 ---\n";
    double mass_cruise = config::MASS_LAUNCH - config::MASS_BOOSTER;
    double W_cruise = mass_cruise * config::g0;
    double alpha_trim = calcTrimAlpha(config::CRUISE_MACH, q_cruise, W_cruise);
    std::cout << "  巡航质量: " << mass_cruise << " kg\n";
    std::cout << "  巡航重力: " << W_cruise << " N\n";
    std::cout << "  动压 (Ma=0.72, H=15m): " << q_cruise << " Pa\n";
    std::cout << "  配平攻角: " << alpha_trim * config::RAD2DEG << " deg\n";
}

// ---------------------------------------------------------------------------
// 阶段二: 动力特性分析
// ---------------------------------------------------------------------------
static void runPropulsionAnalysis() {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "  阶段二: 动力特性分析\n";
    std::cout << "============================================================\n";

    std::string outdir = "output";
    ensureDirectory(outdir);

    // --- 2.1 固体助推器推力特性 ---
    std::cout << "\n--- 2.1 固体助推器推力-时间曲线 ---\n";
    std::cout << std::setw(8) << "t_s"
              << std::setw(12) << "T_booster_N" << "\n";

    std::ofstream ofs_boost(outdir + "/propulsion_booster.csv");
    ofs_boost << "time_s,thrust_N\n";

    for (double t = 0.0; t <= config::BOOSTER_BURN_TIME + 0.5; t += 0.1) {
        double T = calcBoosterThrust(t);
        ofs_boost << std::fixed << std::setprecision(2) << t << ","
                  << std::setprecision(1) << T << "\n";
        if (std::fmod(t * 10, 5.0) < 0.01) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(1) << t
                      << std::setw(12) << std::setprecision(0) << T << "\n";
        }
    }
    ofs_boost.close();

    // 助推器总冲计算
    double total_impulse = 0.0;
    double dt = 0.001;
    for (double t = 0.0; t <= config::BOOSTER_BURN_TIME; t += dt) {
        total_impulse += calcBoosterThrust(t) * dt;
    }
    double Isp_booster = total_impulse / (config::MASS_BOOSTER_PROPELLANT * config::g0);
    std::cout << "\n  助推器总冲: " << total_impulse << " N*s\n";
    std::cout << "  助推器比冲: " << Isp_booster << " s\n";

    // --- 2.2 涡扇发动机推力特性 ---
    std::cout << "\n--- 2.2 涡扇发动机推力随高度变化 (Ma=0.72) ---\n";
    std::cout << std::setw(10) << "H_m"
              << std::setw(12) << "T_engine_N"
              << std::setw(14) << "SFC_kg_N_h" << "\n";

    std::ofstream ofs_engine(outdir + "/propulsion_engine.csv");
    ofs_engine << "altitude_m,thrust_N,sfc_kg_N_h\n";

    double alt_values[] = {0, 15, 50, 60, 100, 200, 500, 1000, 3000, 5000, 10000};
    for (double h : alt_values) {
        AtmosphereState atm = calcAtmosphere(h);
        double V = config::CRUISE_MACH * atm.a;
        double T = calcEngineThrust(h, V, config::CRUISE_MACH);
        double sfc = calcEngineSFC(h, V, config::CRUISE_MACH);
        double sfc_hr = sfc * 3600.0; // 转换为 kg/(N*h)

        ofs_engine << h << "," << T << "," << sfc_hr << "\n";
        std::cout << std::setw(10) << std::fixed << std::setprecision(0) << h
                  << std::setw(12) << std::setprecision(1) << T
                  << std::setw(14) << std::setprecision(4) << sfc_hr << "\n";
    }
    ofs_engine.close();

    // --- 2.3 不同速度下的推力 ---
    std::cout << "\n--- 2.3 涡扇发动机推力随速度变化 (H=15m) ---\n";
    double h_ref = config::CRUISE_ALT_SEA;
    std::cout << std::setw(8) << "Ma"
              << std::setw(12) << "V_m_s"
              << std::setw(12) << "T_engine_N" << "\n";

    double machs[] = {0.3, 0.5, 0.6, 0.65, 0.72, 0.8, 0.85, 0.9};
    for (double Ma : machs) {
        AtmosphereState atm = calcAtmosphere(h_ref);
        double V = Ma * atm.a;
        double T = calcEngineThrust(h_ref, V, Ma);
        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << Ma
                  << std::setw(12) << std::setprecision(1) << V
                  << std::setw(12) << std::setprecision(1) << T << "\n";
    }
}

// ---------------------------------------------------------------------------
// 辅助: 运行工况 → 写 CSV → 打印摘要
// ---------------------------------------------------------------------------
static void runAndSaveCase(const std::string& caseName,
                           const TrajectoryParams& params)
{
    auto records = runCase(params);
    std::string outdir = "output";
    ensureDirectory(outdir);
    writeTrajectoryCSV(outdir + "/trajectory_" + caseName + ".csv", records);
    printTrajectorySummary(caseName, records);
}

// ---------------------------------------------------------------------------
// 阶段三/四: 弹道计算 (C 系列工况)
// ---------------------------------------------------------------------------
static void runTrajectoryCases() {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "  阶段三/四: 弹道计算 (C 系列工况)\n";
    std::cout << "============================================================\n";

    // C-1: 基准工况 (海上巡航)
    runAndSaveCase("C1_baseline_sea", defaultTrajectoryParams());

    // C-2: 陆上巡航
    { TrajectoryParams p = defaultTrajectoryParams(); p.cruise_alt = 60.0;
      runAndSaveCase("C2_land_60m", p); }

    // C-3: 低速巡航
    { TrajectoryParams p = defaultTrajectoryParams(); p.cruise_mach = 0.65;
      runAndSaveCase("C3_Ma065", p); }

    // C-4: 高速巡航
    { TrajectoryParams p = defaultTrajectoryParams(); p.cruise_mach = 0.80;
      runAndSaveCase("C4_Ma080", p); }

    // C-5: 中高度巡航
    { TrajectoryParams p = defaultTrajectoryParams(); p.cruise_alt = 30.0;
      runAndSaveCase("C5_alt30m", p); }

    // C-6: 减重方案
    { TrajectoryParams p = defaultTrajectoryParams(); p.initial_mass = 1350.0;
      p.fuel_mass = 350.0; runAndSaveCase("C6_mass1350", p); }

    // C-7: 高空巡航
    { TrajectoryParams p = defaultTrajectoryParams(); p.cruise_alt = 100.0;
      runAndSaveCase("C7_alt100m", p); }
}

// ---------------------------------------------------------------------------
// §3-三: 弹道特性分析 (D 系列 + 基准工况表现分析)
// ---------------------------------------------------------------------------
static void runParametricStudy() {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "  弹道特性分析 (§3-三)\n";
    std::cout << "============================================================\n";

    std::string outdir = "output";
    ensureDirectory(outdir);

    // ---- 全弹道气动和动力表现分析 (基准工况 C-1) ----
    std::cout << "\n--- 全弹道气动和动力表现分析 ---\n";
    auto c1_records = runCase(defaultTrajectoryParams());
    PhaseAnalysis pa = analyzeBaseline(c1_records);
    printPhaseAnalysis(pa);
    writeTrajectoryCSV(outdir + "/trajectory_C1_baseline_sea.csv", c1_records);
    printTrajectorySummary("C1_baseline_sea", c1_records);

    // ---- D-1: 巡航速度扫描 ----
    std::cout << "\n--- D-1: 巡航速度对射程的影响 (H=15m, m0=1450kg) ---\n";
    auto vel_results = sweepVelocity(0.55, 0.85, 0.05, 15.0, config::MASS_LAUNCH);
    writeSweepCSV(outdir + "/sweep_velocity.csv", vel_results, "Ma,V_m_s");

    // D-1 结果表
    std::cout << "\n  " << std::string(60, '-') << "\n";
    std::cout << "  " << std::left << std::setw(6) << "Ma"
              << std::setw(8) << "V(m/s)"
              << std::setw(10) << "射程(km)"
              << std::setw(10) << "时间(min)"
              << std::setw(8) << "L/D" << "\n";
    for (const auto& sr : vel_results) {
        std::cout << "  " << std::fixed << std::setprecision(2) << sr.param_value
                  << "  " << std::setprecision(1) << sr.V_cruise_ms
                  << "  " << std::setprecision(1) << sr.range_km
                  << "    " << std::setprecision(1) << sr.flight_time_min
                  << "    " << std::setprecision(2) << sr.avg_LD << "\n";
    }

    // ---- D-2: 巡航高度扫描 ----
    std::cout << "\n--- D-2: 巡航高度对射程的影响 (Ma=0.72, m0=1450kg) ---\n";
    std::vector<double> altitudes = {7, 15, 30, 60, 100, 200};
    auto alt_results = sweepAltitude(altitudes, config::CRUISE_MACH, config::MASS_LAUNCH);
    writeSweepCSV(outdir + "/sweep_altitude.csv", alt_results, "altitude_m");

    std::cout << "\n  " << std::string(55, '-') << "\n";
    std::cout << "  " << std::setw(8) << "H(m)"
              << std::setw(10) << "射程(km)"
              << std::setw(10) << "时间(min)"
              << std::setw(8) << "L/D" << "\n";
    for (const auto& sr : alt_results) {
        std::cout << "  " << std::fixed << std::setprecision(0) << sr.param_value
                  << "    " << std::setprecision(1) << sr.range_km
                  << "    " << std::setprecision(1) << sr.flight_time_min
                  << "    " << std::setprecision(2) << sr.avg_LD << "\n";
    }

    // ---- D-3: 初始质量扫描 ----
    std::cout << "\n--- D-3: 初始质量对射程的影响 (Ma=0.72, H=15m) ---\n";
    auto mass_results = sweepMass(1200.0, 1600.0, 50.0,
                                   config::CRUISE_MACH, 15.0);
    writeSweepCSV(outdir + "/sweep_mass.csv", mass_results, "mass_kg");

    std::cout << "\n  " << std::string(55, '-') << "\n";
    std::cout << "  " << std::setw(10) << "m0(kg)"
              << std::setw(10) << "射程(km)"
              << std::setw(10) << "时间(min)"
              << std::setw(8) << "L/D" << "\n";
    for (const auto& sr : mass_results) {
        std::cout << "  " << std::fixed << std::setprecision(0) << sr.param_value
                  << "      " << std::setprecision(1) << sr.range_km
                  << "    " << std::setprecision(1) << sr.flight_time_min
                  << "    " << std::setprecision(2) << sr.avg_LD << "\n";
    }

    std::cout << "\n  扫描结果已保存至 " << outdir << "/sweep_*.csv\n";
}

// ---------------------------------------------------------------------------
// 主函数
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    std::cout << "============================================================\n";
    std::cout << "  BGM-109 \"战斧\" 巡航导弹质点弹道仿真程序\n";
    std::cout << "  四阶龙格-库塔法 | 二维纵向质点弹道模型\n";
    std::cout << "============================================================\n";

    // 打印基本参数
    std::cout << "\n--- 导弹基本参数 ---\n";
    std::cout << "  全长:       " << config::LENGTH << " m\n";
    std::cout << "  弹体直径:   " << config::DIAMETER << " m\n";
    std::cout << "  翼展:       " << config::WINGSPAN << " m\n";
    std::cout << "  机翼面积:   " << config::S_WING << " m^2\n";
    std::cout << "  展弦比:     " << config::AR << "\n";
    std::cout << "  发射质量:   " << config::MASS_LAUNCH << " kg\n";
    std::cout << "  巡航Ma:     " << config::CRUISE_MACH << "\n";
    std::cout << "  涡扇推力:   " << config::ENGINE_THRUST << " N\n";
    std::cout << "  助推器推力: " << config::BOOSTER_THRUST << " N\n";

    ensureDirectory("output");

    // 根据命令行参数选择运行模式
    std::string mode = (argc > 1) ? argv[1] : "all";

    if (mode == "all" || mode == "aero") {
        runAeroAnalysis();
    }
    if (mode == "all" || mode == "propulsion") {
        runPropulsionAnalysis();
    }
    if (mode == "all" || mode == "trajectory") {
        runTrajectoryCases();
    }
    if (mode == "all" || mode == "sweep" || mode == "baseline") {
        runParametricStudy();
    }

    std::cout << "\n\n============================================================\n";
    std::cout << "  全部计算完成! 结果保存在 output/ 目录下\n";
    std::cout << "============================================================\n";

    return 0;
}
