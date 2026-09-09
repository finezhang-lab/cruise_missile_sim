#pragma once

// ============================================================================
// integrator.h -- 四阶龙格-库塔 (RK4) 数值积分器
// 通用实现，适用于任意 ODE 系统
// ============================================================================

#include "trajectory.h"
#include <vector>
#include <functional>

// 右端函数类型: f(state, t, params) -> dstate/dt
using DerivativeFunc = std::function<StateVec(const StateVec&, double, const TrajectoryParams&)>;

// 单步 RK4 积分
// state: 当前状态
// t: 当前时间
// dt: 时间步长
// f: 右端函数
// params: 弹道参数
// 返回: 下一步状态
StateVec rk4Step(const StateVec& state, double t, double dt,
                 const TrajectoryParams& params,
                 DerivativeFunc f);

// 全弹道积分
// 从 t_start 积分到终止条件满足
// initial_state: 初始状态
// t_start: 起始时间
// dt: 时间步长 (助推段/巡航段自动切换)
// params: 工况参数
// 返回: 完整弹道记录序列
std::vector<TrajectoryRecord> integrateTrajectory(
    const StateVec& initial_state,
    double t_start,
    const TrajectoryParams& params,
    DerivativeFunc f);
