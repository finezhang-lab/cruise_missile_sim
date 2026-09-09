#pragma once

// ============================================================================
// parametric_study.h -- 弹道特性分析模块 (§3-三)
//
// 职责:
//   1. 控制参数分析: 巡航速度/高度/质量对射程的影响
//   2. 全弹道气动和动力表现分析: 基准工况下各物理量统计
//   3. 结合升阻特性和动力特性解释参数影响规律
//
// 依赖: integrator, atmosphere, aerodynamics (只读)
// ============================================================================

#include "trajectory.h"
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// 参数扫描单点结果
// ---------------------------------------------------------------------------
struct SweepResult {
    double param_value;    // 扫描参数值 (Ma / H / m0)
    double V_cruise_ms;   // 巡航速度 [m/s]  (仅 D-1 有意义)
    double range_km;       // 总射程 [km]
    double flight_time_min;// 飞行时间 [min]
    double fuel_consumed_kg;// 燃油消耗 [kg]
    double avg_LD;         // 巡航段平均升阻比
    double avg_sfc_hr;     // 巡航段平均耗油率 [kg/(N·h)]
};

// ---------------------------------------------------------------------------
// 基准工况全弹道表现分析
// ---------------------------------------------------------------------------
struct PhaseAnalysis {
    // 助推段
    double boost_time_s;       // 助推时间 [s]
    double boost_end_V_ms;     // 助推段末速度 [m/s]
    double boost_max_accel_g;  // 最大加速度 [g]

    // 巡航段
    double cruise_avg_alpha_deg;// 平均配平攻角 [deg]
    double cruise_avg_LD;      // 平均升阻比
    double cruise_avg_T_N;     // 平均推力 [N]
    double cruise_avg_D_N;     // 平均阻力 [N]
    double cruise_avg_V_ms;    // 平均速度 [m/s]
    double cruise_avg_h_m;     // 平均高度 [m]

    // 全局
    double total_range_km;     // 总射程 [km]
    double total_flight_min;   // 总飞行时间 [min]
    double total_fuel_kg;      // 总燃油消耗 [kg]
};

// ---------------------------------------------------------------------------
// 运行单个弹道工况 (不写文件，返回记录)
// ---------------------------------------------------------------------------
std::vector<TrajectoryRecord> runCase(const TrajectoryParams& params);

// ---------------------------------------------------------------------------
// 基准工况全弹道表现分析 (对 C-1 结果做统计)
// ---------------------------------------------------------------------------
PhaseAnalysis analyzeBaseline(const std::vector<TrajectoryRecord>& records);

// ---------------------------------------------------------------------------
// 参数扫描 (D 系列)
// ---------------------------------------------------------------------------

// D-1: 巡航速度扫描
std::vector<SweepResult> sweepVelocity(double Ma_min, double Ma_max, double step,
                                       double H, double m0);

// D-2: 巡航高度扫描
std::vector<SweepResult> sweepAltitude(const std::vector<double>& altitudes,
                                       double Ma, double m0);

// D-3: 初始质量扫描
std::vector<SweepResult> sweepMass(double m0_min, double m0_max, double step,
                                    double Ma, double H);

// ---------------------------------------------------------------------------
// 输出函数
// ---------------------------------------------------------------------------
void printPhaseAnalysis(const PhaseAnalysis& pa);
void writeSweepCSV(const std::string& filename, const std::vector<SweepResult>& results,
                   const std::string& param_header);
