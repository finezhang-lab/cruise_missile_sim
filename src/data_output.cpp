// ============================================================================
// data_output.cpp -- 数据输出模块实现
// ============================================================================

#include "data_output.h"
#include "config.h"
#include <fstream>
#include <iomanip>
#include <iostream>

// ---------------------------------------------------------------------------
// 写入 CSV 文件
// ---------------------------------------------------------------------------
bool writeTrajectoryCSV(const std::string& filename,
                        const std::vector<TrajectoryRecord>& records) {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cerr << "错误: 无法打开文件 " << filename << " 进行写入\n";
        return false;
    }

    // 写入表头
    ofs << "time_s"
        << ",V_m_s"
        << ",theta_deg"
        << ",x_m"
        << ",y_m"
        << ",mass_kg"
        << ",L_N"
        << ",D_N"
        << ",T_N"
        << ",Ma"
        << ",alpha_deg"
        << ",phase"
        << "\n";

    // 写入数据
    ofs << std::fixed;
    for (const auto& r : records) {
        ofs << std::setprecision(3) << r.t << ","
            << std::setprecision(2) << r.V << ","
            << std::setprecision(4) << r.theta << ","
            << std::setprecision(2) << r.x << ","
            << std::setprecision(2) << r.y << ","
            << std::setprecision(2) << r.m << ","
            << std::setprecision(2) << r.L << ","
            << std::setprecision(2) << r.D << ","
            << std::setprecision(2) << r.T << ","
            << std::setprecision(4) << r.Ma << ","
            << std::setprecision(4) << r.alpha << ","
            << std::setprecision(0) << r.phase
            << "\n";
    }

    ofs.close();
    return true;
}

// ---------------------------------------------------------------------------
// 打印弹道摘要
// ---------------------------------------------------------------------------
void printTrajectorySummary(const std::string& caseName,
                            const std::vector<TrajectoryRecord>& records) {
    if (records.empty()) {
        std::cout << "  [" << caseName << "] 无数据\n";
        return;
    }

    const auto& first = records.front();
    const auto& last  = records.back();

    // 寻找最大高度和最大速度
    double max_alt = 0.0, max_V = 0.0;
    for (const auto& r : records) {
        if (r.y > max_alt) max_alt = r.y;
        if (r.V > max_V) max_V = r.V;
    }

    double range_km = last.x / 1000.0;
    double flight_time_min = last.t / 60.0;

    std::cout << "\n========== 工况: " << caseName << " ==========\n";
    std::cout << "  飞行时间:     " << last.t << " s (" << flight_time_min << " min)\n";
    std::cout << "  射程:         " << range_km << " km\n";
    std::cout << "  最大高度:     " << max_alt << " m\n";
    std::cout << "  最大速度:     " << max_V << " m/s (Ma "
              << max_V / 340.3 << ")\n";
    std::cout << "  终止速度:     " << last.V << " m/s\n";
    std::cout << "  终止质量:     " << last.m << " kg\n";
    std::cout << "  数据点数:     " << records.size() << "\n";
    std::cout << "=============================================\n";
}
