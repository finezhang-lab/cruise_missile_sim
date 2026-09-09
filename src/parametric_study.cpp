// ============================================================================
// parametric_study.cpp -- 弹道特性分析模块实现
//
// 对应课题要求 §3-三:
//   - 控制参数分析 (D-1/D-2/D-3 扫描)
//   - 全弹道气动和动力表现分析 (基准工况统计)
//   - 不同巡航速度与高度对射程的影响
// ============================================================================

#include "parametric_study.h"
#include "integrator.h"
#include "atmosphere.h"
#include "propulsion.h"
#include "config.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// 运行单个弹道工况 (仅积分，不写文件)
// ---------------------------------------------------------------------------
std::vector<TrajectoryRecord> runCase(const TrajectoryParams& params) {
    StateVec x0{};
    // 弹道计算从助推器分离后的巡航状态开始，避免 V=0 奇点
    x0[IDX_V]     = 120.0;                              // 助推段末速度 [m/s]
    x0[IDX_THETA] = 0.0;                                // 水平飞行
    x0[IDX_X]     = 0.0;
    x0[IDX_Y]     = params.cruise_alt;                  // 直接到达巡航高度
    x0[IDX_M]     = config::MASS_LAUNCH - config::MASS_BOOSTER;  // 1150 kg

    return integrateTrajectory(x0, 0.0, params, trajectoryDerivatives);
}

// ---------------------------------------------------------------------------
// 基准工况全弹道表现分析
// ---------------------------------------------------------------------------
PhaseAnalysis analyzeBaseline(const std::vector<TrajectoryRecord>& records) {
    PhaseAnalysis pa{};

    if (records.empty()) return pa;

    // ---- 分段：助推段 (phase==0) vs 巡航段 (phase==1) ----
    int n_boost = 0, n_cruise = 0;
    double sum_alpha = 0, sum_LD = 0, sum_T = 0, sum_D = 0;
    double sum_V = 0, sum_h = 0;
    double boost_max_dVdt = 0.0;
    double boost_end_V = 0.0;

    for (size_t i = 0; i < records.size(); ++i) {
        const auto& r = records[i];

        if (r.phase == 0) {
            // 助推段
            ++n_boost;
            boost_end_V = r.V;

            // 用相邻记录估算瞬时加速度
            if (i > 0) {
                double dt = r.t - records[i - 1].t;
                if (dt > 0 && records[i - 1].phase == 0) {
                    double accel = (r.V - records[i - 1].V) / dt;
                    if (accel > boost_max_dVdt) {
                        boost_max_dVdt = accel;
                    }
                }
            }
        } else {
            // 巡航段
            ++n_cruise;
            sum_alpha += r.alpha;

            double ld = (r.D > 0.1) ? (r.L / r.D) : 0.0;
            sum_LD += ld;
            sum_T  += r.T;
            sum_D  += r.D;
            sum_V  += r.V;
            sum_h  += r.y;
        }
    }

    // ---- 助推段统计 ----
    pa.boost_end_V_ms    = boost_end_V;
    pa.boost_max_accel_g = boost_max_dVdt / config::g0;
    if (n_boost > 0) {
        // 找助推段最后一条记录的时间
        for (auto rit = records.rbegin(); rit != records.rend(); ++rit) {
            if (rit->phase == 0) {
                pa.boost_time_s = rit->t;
                break;
            }
        }
    }

    // ---- 巡航段统计 ----
    if (n_cruise > 0) {
        pa.cruise_avg_alpha_deg = (sum_alpha / n_cruise);
        pa.cruise_avg_LD       = sum_LD  / n_cruise;
        pa.cruise_avg_T_N     = sum_T    / n_cruise;
        pa.cruise_avg_D_N     = sum_D    / n_cruise;
        pa.cruise_avg_V_ms    = sum_V    / n_cruise;
        pa.cruise_avg_h_m     = sum_h    / n_cruise;
    }

    // ---- 全局 ----
    const auto& last = records.back();
    pa.total_range_km    = last.x / 1000.0;
    pa.total_flight_min  = last.t / 60.0;
    pa.total_fuel_kg     = (config::MASS_LAUNCH - config::MASS_BOOSTER) - last.m;

    return pa;
}

// ---------------------------------------------------------------------------
// 从弹道记录中提取扫描结果 (公用逻辑)
// ---------------------------------------------------------------------------
static SweepResult extractSweepResult(double param_value,
                                      const std::vector<TrajectoryRecord>& records,
                                      double V_cruise_ms)
{
    SweepResult sr{};
    sr.param_value     = param_value;
    sr.V_cruise_ms     = V_cruise_ms;

    if (records.empty()) {
        sr.range_km = sr.flight_time_min = sr.fuel_consumed_kg = 0;
        sr.avg_LD = sr.avg_sfc_hr = 0;
        return sr;
    }

    const auto& last = records.back();
    sr.range_km        = last.x / 1000.0;
    sr.flight_time_min = last.t / 60.0;
    sr.fuel_consumed_kg = (config::MASS_LAUNCH - config::MASS_BOOSTER) - last.m;

    // 巡航段平均 L/D 和 sfc
    int n_cruise = 0;
    double sum_LD = 0, sum_sfc_hr = 0;
    for (const auto& r : records) {
        if (r.phase == 1) {
            ++n_cruise;
            double ld = (r.D > 0.1) ? (r.L / r.D) : 0.0;
            sum_LD += ld;

            // 从记录中还原 sfc
            double sfc = calcEngineSFC(r.y, r.V, r.Ma);
            sum_sfc_hr += sfc * 3600.0; // kg/(N·s) → kg/(N·h)
        }
    }
    sr.avg_LD     = (n_cruise > 0) ? sum_LD / n_cruise : 0.0;
    sr.avg_sfc_hr = (n_cruise > 0) ? sum_sfc_hr / n_cruise : 0.0;

    return sr;
}

// ---------------------------------------------------------------------------
// D-1: 巡航速度扫描
// ---------------------------------------------------------------------------
std::vector<SweepResult> sweepVelocity(double Ma_min, double Ma_max, double step,
                                       double H, double m0)
{
    std::vector<SweepResult> results;

    for (double Ma = Ma_min; Ma <= Ma_max + 1e-9; Ma += step) {
        AtmosphereState atm = calcAtmosphere(H);
        double V = Ma * atm.a;

        TrajectoryParams p = defaultTrajectoryParams();
        p.cruise_mach = Ma;
        p.initial_mass = m0;

        std::cout << "  [D-1] Ma=" << std::fixed << std::setprecision(2) << Ma << " ... " << std::flush;
        auto records = runCase(p);
        auto sr = extractSweepResult(Ma, records, V);
        results.push_back(sr);

        std::cout << "射程=" << std::setprecision(1) << sr.range_km << " km"
                  << ", 时间=" << std::setprecision(1) << sr.flight_time_min << " min\n";
    }
    return results;
}

// ---------------------------------------------------------------------------
// D-2: 巡航高度扫描
// ---------------------------------------------------------------------------
std::vector<SweepResult> sweepAltitude(const std::vector<double>& altitudes,
                                       double Ma, double m0)
{
    std::vector<SweepResult> results;

    for (double h : altitudes) {
        AtmosphereState atm = calcAtmosphere(h);
        double V = Ma * atm.a;

        TrajectoryParams p = defaultTrajectoryParams();
        p.cruise_alt = h;
        p.initial_mass = m0;

        std::cout << "  [D-2] H=" << std::fixed << std::setprecision(0) << h << " m ... " << std::flush;
        auto records = runCase(p);
        auto sr = extractSweepResult(h, records, V);
        results.push_back(sr);

        std::cout << "射程=" << std::setprecision(1) << sr.range_km << " km"
                  << ", 时间=" << std::setprecision(1) << sr.flight_time_min << " min\n";
    }
    return results;
}

// ---------------------------------------------------------------------------
// D-3: 初始质量扫描
// ---------------------------------------------------------------------------
std::vector<SweepResult> sweepMass(double m0_min, double m0_max, double step,
                                    double Ma, double H)
{
    std::vector<SweepResult> results;

    for (double m0 = m0_min; m0 <= m0_max + 1e-9; m0 += step) {
        TrajectoryParams p = defaultTrajectoryParams();
        p.initial_mass = m0;
        // 燃油质量联动调整: 增减的质量全部算在燃油上
        p.fuel_mass = config::MASS_FUEL - (config::MASS_LAUNCH - m0);
        if (p.fuel_mass < 50.0) p.fuel_mass = 50.0;

        AtmosphereState atm = calcAtmosphere(H);
        double V = Ma * atm.a;

        std::cout << "  [D-3] m0=" << std::fixed << std::setprecision(0) << m0 << " kg ... " << std::flush;
        auto records = runCase(p);
        auto sr = extractSweepResult(m0, records, V);
        results.push_back(sr);

        std::cout << "射程=" << std::setprecision(1) << sr.range_km << " km"
                  << ", 时间=" << std::setprecision(1) << sr.flight_time_min << " min\n";
    }
    return results;
}

// ---------------------------------------------------------------------------
// 输出: 基准工况分析打印
// ---------------------------------------------------------------------------
void printPhaseAnalysis(const PhaseAnalysis& pa) {
    std::cout << "\n";
    std::cout << "  === 全弹道气动和动力表现分析 (基准工况 C-1) ===\n\n";

    std::cout << "  --- 助推段 ---\n";
    std::cout << "    助推时间:       " << std::fixed << std::setprecision(1)
              << pa.boost_time_s << " s\n";
    std::cout << "    助推段末速度:   " << std::setprecision(1)
              << pa.boost_end_V_ms << " m/s  ("
              << std::setprecision(2) << pa.boost_end_V_ms / 343.0 << " Ma)\n";
    std::cout << "    最大加速度:     " << std::setprecision(2)
              << pa.boost_max_accel_g << " g\n";

    std::cout << "\n  --- 巡航段 ---\n";
    std::cout << "    平均配平攻角:   " << std::setprecision(2)
              << pa.cruise_avg_alpha_deg << " deg\n";
    std::cout << "    平均升阻比:     " << std::setprecision(2)
              << pa.cruise_avg_LD << "\n";
    std::cout << "    平均推力:       " << std::setprecision(0)
              << pa.cruise_avg_T_N << " N\n";
    std::cout << "    平均阻力:       " << std::setprecision(0)
              << pa.cruise_avg_D_N << " N\n";
    std::cout << "    平均巡航速度:   " << std::setprecision(1)
              << pa.cruise_avg_V_ms << " m/s\n";
    std::cout << "    平均巡航高度:   " << std::setprecision(0)
              << pa.cruise_avg_h_m << " m\n";

    std::cout << "\n  --- 全局 ---\n";
    std::cout << "    总射程:         " << std::setprecision(1)
              << pa.total_range_km << " km\n";
    std::cout << "    总飞行时间:     " << std::setprecision(1)
              << pa.total_flight_min << " min  ("
              << std::setprecision(2) << pa.total_flight_min / 60.0 << " h)\n";
    std::cout << "    总燃油消耗:     " << std::setprecision(1)
              << pa.total_fuel_kg << " kg\n";

    // 推阻平衡分析
    double TD_ratio = (pa.cruise_avg_D_N > 0.1)
        ? pa.cruise_avg_T_N / pa.cruise_avg_D_N : 0.0;
    std::cout << "\n    推阻比 T/D:     " << std::setprecision(3) << TD_ratio;
    if (std::abs(TD_ratio - 1.0) < 0.05)
        std::cout << " (推力≈阻力，巡航平衡良好)\n";
    else if (TD_ratio > 1.05)
        std::cout << " (推力>阻力，有加速余量)\n";
    else
        std::cout << " (推力<阻力，巡航可能减速)\n";
}

// ---------------------------------------------------------------------------
// 输出: 扫描结果写入 CSV
// ---------------------------------------------------------------------------
void writeSweepCSV(const std::string& filename,
                   const std::vector<SweepResult>& results,
                   const std::string& param_header)
{
    std::ofstream ofs(filename);
    ofs << param_header << ",range_km,flight_time_min,fuel_consumed_kg,avg_LD,avg_sfc_hr\n";

    for (const auto& sr : results) {
        ofs << std::fixed << std::setprecision(3) << sr.param_value << ","
            << std::setprecision(2) << sr.range_km << ","
            << std::setprecision(2) << sr.flight_time_min << ","
            << std::setprecision(1) << sr.fuel_consumed_kg << ","
            << std::setprecision(3) << sr.avg_LD << ","
            << std::setprecision(4) << sr.avg_sfc_hr << "\n";
    }
    ofs.close();
}
