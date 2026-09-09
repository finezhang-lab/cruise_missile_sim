#pragma once

// ============================================================================
// propulsion.h -- 动力系统模型
// 包含: 固体火箭助推器推力模型 + 涡扇发动机推力/耗油模型
// 飞行阶段: 助推段 -> 助推器分离 -> 巡航段
// ============================================================================

// 飞行阶段枚举
enum class FlightPhase {
    BOOST,       // 助推段（固体助推器工作）
    SEPARATION,  // 助推器分离瞬间
    CRUISE       // 巡航段（涡扇发动机工作）
};

struct PropulsionResult {
    double thrust;       // 总推力 [N]
    double massFlowRate; // 质量流率 [kg/s] (正值, dm/dt = -massFlowRate)
    FlightPhase phase;   // 当前飞行阶段
};

// 获取当前飞行阶段
FlightPhase getFlightPhase(double t, double mass);

// 计算助推器推力 (恒推力近似)
// t: 从发射起算时间 [s]
double calcBoosterThrust(double t);

// 计算涡扇发动机推力
// h: 飞行高度 [m]
// V: 飞行速度 [m/s]
// Ma: 马赫数
double calcEngineThrust(double h, double V, double Ma);

// 计算涡扇发动机耗油率 [kg/(N*s)]
// 可根据高度和速度做修正
double calcEngineSFC(double h, double V, double Ma);

// 综合推力计算: 根据飞行阶段自动选择动力源
PropulsionResult calcPropulsion(double t, double h, double V, double Ma, double mass);
