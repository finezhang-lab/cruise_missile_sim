# BGM-109 "战斧" 巡航导弹质点弹道仿真

> 基于二维纵向质点运动方程，采用经典四阶龙格-库塔法（RK4）数值积分，模拟战斧巡航导弹从发射到燃油耗尽的完整飞行过程，分析气动特性、动力特性及弹道参数影响规律。

## 项目结构

```
cruise_missile_sim/
├── CMakeLists.txt                # CMake 构建配置
├── include/                      # 头文件（模块接口）
│   ├── config.h                  # 全局参数定义（BGM-109 导弹参数、物理常数）
│   ├── atmosphere.h              # ISA 1976 标准大气模型
│   ├── aerodynamics.h            # 气动力计算（CL/CD/LD/静稳/配平）
│   ├── propulsion.h              # 动力系统（助推器 + 涡扇发动机）
│   ├── trajectory.h              # 质点运动方程（5 状态量右端函数）
│   ├── integrator.h              # RK4 数值积分器
│   ├── parametric_study.h        # 弹道特性分析（参数扫描 + 基准分析）
│   └── data_output.h             # CSV 数据输出
├── src/                          # 实现文件
│   ├── atmosphere.cpp
│   ├── aerodynamics.cpp
│   ├── propulsion.cpp
│   ├── trajectory.cpp
│   ├── integrator.cpp
│   ├── parametric_study.cpp
│   ├── data_output.cpp
│   └── main.cpp                  # 主程序（仅调度，不含计算逻辑）
└── analysis/
    └── plot.py                   # Python 后处理绘图（matplotlib）
```

## 构建与运行

### 环境要求

- C++17 编译器（GCC ≥ 7 / MSVC ≥ 2019 / Clang ≥ 5）
- CMake ≥ 3.14
- Python 3（绘图可选，需 matplotlib）

### 编译

```bash
cd cruise_missile_sim
mkdir build && cd build
cmake ..
cmake --build .
```

### 运行

```bash
./bin/cruise_missile_sim              # 运行全部分析
./bin/cruise_missile_sim aero         # 仅阶段一：升阻特性分析
./bin/cruise_missile_sim propulsion   # 仅阶段二：动力特性分析
./bin/cruise_missile_sim trajectory   # 仅阶段三四：弹道计算（C 系列）
./bin/cruise_missile_sim sweep        # 仅弹道特性分析（D 系列 + 基准分析）
```

计算结果输出至 `output/` 目录（CSV 格式）。

### 绘图

```bash
python analysis/plot.py
```

## 模块说明与分工

### 模块依赖关系

```
config.h （全局参数）
    ↓
atmosphere ──→ propulsion ──→ trajectory ──→ integrator
    ↓              ↓               ↓              ↓
aerodynamics ──→ trajectory    parametric_study ←──┘
                                      ↓
                                data_output
                                      ↓
                                   main.cpp
```

**原则**：每个模块只通过 `.h` 头文件接口通信，禁止跨模块直接访问 `.cpp` 内部实现。

### 团队分工

| 成员 | 学院 / 专业 / 年级 | 负责模块 | 对应课题阶段 | 主要职责 |
|---|---|---|---|---|
| **张帆** | 航空航天 / 机械 / 博二 | `config.h` + 全部 `.h` 接口 | 全部 | 程序架构设计、代码合并与审查、报告修订 |
| **杜立成** | 航空航天 / 力学 / 博三 | `trajectory.h/cpp` + `integrator.h/cpp` | 阶段三、四 | 质点运动方程实现、RK4 积分器、飞行阶段切换 |
| **姚飞** | 航空航天 / 工程力学 / 大二 | `propulsion.h/cpp` | 阶段二 | 固体助推器推力模型、涡扇推力/SFC 修正模型 |
| **曾刘彰文** | 竺可桢学院 / 工程力学 / 大二 | `atmosphere.h/cpp` + `aerodynamics.h/cpp` | 阶段一 | ISA 标准大气、升力/阻力系数、静稳定性、配平攻角 |
| **楼宇凡** | 竺可桢学院 / 工程力学 / 大二 | `parametric_study.h/cpp` + `data_output.h/cpp` | §3-三 | 基准工况分析、D 系列参数扫描、CSV 输出 |

### 公共模块

| 文件 | 负责人 | 说明 |
|---|---|---|
| `config.h` | 张帆 | 导弹参数、物理常数，其他人只读 |
| `main.cpp` | 张帆 | 调度入口，不含计算逻辑 |
| `CMakeLists.txt` | 张帆 | 构建配置 |
| `analysis/plot.py` | 楼宇凡 | 后处理绘图 |

## 开发规范

### Git 工作流

```
日常开发 → dev 分支（自由推送）
发版审核 → Pull Request: dev → main（需审核合并）
```

```bash
# 首次克隆
git clone https://github.com/finezhang-lab/cruise_missile_sim.git
cd cruise_missile_sim
git checkout dev

# 日常开发
git pull origin dev          # 拉取最新代码
# ... 编辑代码 ...
git add .
git commit -m "描述: 简要说明本次改动"
git push origin dev          # 推送到 dev
```

**禁止直接推送到 `main` 分支**，所有合并必须通过 Pull Request。

### Commit 规范

采用以下前缀格式：

```
feat: 新增 xxx 功能
fix: 修复 xxx 问题
refactor: 重构 xxx 模块
docs: 更新文档
test: 添加测试
```

### 命名规范

| 类别 | 风格 | 示例 |
|---|---|---|
| 函数名 | camelCase | `calcAerodynamics()`, `calcBoosterThrust()` |
| 结构体名 | PascalCase | `AeroResult`, `TrajectoryParams`, `StateVec` |
| 常量 | UPPER_SNAKE_CASE | `config::MASS_LAUNCH`, `config::CRUISE_MACH` |
| 文件名 | snake_case | `parametric_study.cpp`, `data_output.h` |

### 单位制

全部采用 **SI 国际单位制**，内部计算一律使用弧度：

| 物理量 | 单位 |
|---|---|
| 长度 / 距离 / 高度 | m |
| 质量 | kg |
| 力（升力/阻力/推力/重力） | N |
| 速度 | m/s |
| 时间 | s |
| 角度（攻角/弹道倾角） | rad（内部），deg（输出/显示） |
| 耗油率 | kg/(N·s)（内部），kg/(N·h)（输出/显示） |

### 模块开发要求

1. **接口先行**：所有 `.h` 头文件由架构负责人统一定义，开发过程中不得修改函数签名
2. **单向依赖**：模块之间只通过头文件调用，禁止 `#include` 对方的 `.cpp` 文件
3. **独立可测**：每个模块应能单独编写单元测试验证基本正确性
4. **桩函数联调**：依赖模块未完成时，可用返回常数的桩函数（stub）进行本地测试，联调时替换为真实实现

### 调试检查

计算结果应通过以下合理性检验：

| 检查项 | 预期范围 |
|---|---|
| 巡航配平攻角 | 2° ~ 5° |
| 巡航升阻比 L/D | 10 ~ 25 |
| 助推段末速度 | 60 ~ 120 m/s |
| 总射程 | 500 ~ 2500 km |
| 总飞行时间 | 30 ~ 240 min |

若计算结果偏离上述范围，优先检查单位换算和物理模型实现。

## 许可证

本项目为浙江大学航空航天学院课程课题，仅供学术使用。
