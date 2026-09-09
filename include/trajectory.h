#pragma once

// ============================================================================
// trajectory.h -- 导弹质点弹道运动方程 (二维纵向模型)
// 状态向量: x = [V, theta, x_pos, y_pos, m]
//   V     : 飞行速度 [m/s]
//   theta : 弹道倾角 [rad] (相对于水平面)
//   x_pos : 水平距离 [m]
//   y_pos : 高度 [m]
//   m     : 导弹质量 [kg]
// ============================================================================

#include <array>
#include <vector>

// 状态向量维度
constexpr int STATE_DIM = 5;
using StateVec = std::array<double, STATE_DIM>;

// 状态分量索引
enum StateIdx {
    IDX_V     = 0, // 速度
    IDX_THETA = 1, // 弹道倾角
    IDX_X     = 2, // 水平距离
    IDX_Y     = 3, // 高度
    IDX_M     = 4  // 质量
};

// 单步记录 (用于输出)
struct TrajectoryRecord {
    double t;      // 时间 [s]
    double V;      // 速度 [m/s]
    double theta;  // 弹道倾角 [deg]
    double x;      // 水平距离 [m]
    double y;      // 高度 [m]
    double m;      // 质量 [kg]
    double L;      // 升力 [N]
    double D;      // 阻力 [N]
    double T;      // 推力 [N]
    double Ma;     // 马赫数
    double alpha;  // 攻角 [deg]
    double phase;  // 飞行阶段 (0=助推, 1=巡航)
};

// 弹道计算参数 (工况配置)
struct TrajectoryParams {
    double cruise_alt;     // 巡航高度 [m]
    double cruise_mach;    // 巡航马赫数
    double initial_mass;   // 发射质量 [kg]
    double fuel_mass;      // 巡航段可用燃油 [kg]
    double cruise_alpha;   // 巡航攻角 [rad] (若为 0 则自动计算配平攻角)
};

// 计算状态向量的时间导数 (右端函数)
// state: 当前状态 [V, theta, x, y, m]
// t: 当前时间 [s]
// params: 工况参数
// 返回: dstate/dt
StateVec trajectoryDerivatives(const StateVec& state, double t, const TrajectoryParams& params);

// 检查终止条件
// 返回 true 表示仿真应终止 (燃油耗尽 / 触地 / 超时)
bool checkTermination(const StateVec& state, double t, const TrajectoryParams& params,
                      double initial_fuel_mass);

// 生成默认工况参数
TrajectoryParams defaultTrajectoryParams();
