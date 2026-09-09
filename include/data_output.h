#pragma once

// ============================================================================
// data_output.h -- 数据输出模块 (CSV 格式)
// 将弹道计算结果输出为 CSV 文件，便于后续用 Python/MATLAB/Excel 分析
// ============================================================================

#include "trajectory.h"
#include <vector>
#include <string>

// 将弹道记录写入 CSV 文件
// filename: 输出文件路径
// records: 弹道记录序列
// 返回: 是否写入成功
bool writeTrajectoryCSV(const std::string& filename,
                        const std::vector<TrajectoryRecord>& records);

// 打印弹道摘要信息到控制台
void printTrajectorySummary(const std::string& caseName,
                            const std::vector<TrajectoryRecord>& records);
