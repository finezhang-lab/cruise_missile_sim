#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot.py -- BGM-109 巡航导弹弹道仿真后处理绘图脚本
读取 C++ 程序输出的 CSV 文件，绘制各类弹道曲线图

使用方法:
  python plot.py               # 绘制全部图表
  python plot.py trajectory    # 仅绘制弹道曲线
  python plot.py aero          # 仅绘制气动特性
  python plot.py sweep         # 仅绘制参数扫描结果

依赖: matplotlib, pandas, numpy
安装: pip install matplotlib pandas numpy
"""

import os
import sys
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# 中文字体设置 (Windows)
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

# 输出目录: 兼容从项目根目录运行和从 build/bin 运行
# 优先选择包含 CSV 数据文件的目录
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_SCRIPT_DIR)
OUTPUT_DIR_CANDIDATES = [
    os.path.join(_PROJECT_ROOT, 'output'),
    os.path.join(_PROJECT_ROOT, 'build', 'bin', 'output'),
]

def _pick_output_dir():
    candidates_with_data = [p for p in OUTPUT_DIR_CANDIDATES
                            if os.path.isdir(p) and glob.glob(os.path.join(p, '*.csv'))]
    if candidates_with_data:
        return candidates_with_data[0]
    existing = [p for p in OUTPUT_DIR_CANDIDATES if os.path.isdir(p)]
    return existing[0] if existing else OUTPUT_DIR_CANDIDATES[0]

OUTPUT_DIR = _pick_output_dir()
FIGURE_DIR = os.path.join(os.path.dirname(__file__), 'figures')

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

# ============================================================================
# 1. 气动特性绘图
# ============================================================================
def plot_aero():
    ensure_dir(FIGURE_DIR)
    print("绘制气动特性图...")

    # 1.1 CL, CD 随攻角变化
    fpath = os.path.join(OUTPUT_DIR, 'aero_CL_CD.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)

        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        # CL vs alpha
        axes[0].plot(df['alpha_deg'], df['CL'], 'b-', linewidth=2)
        axes[0].set_xlabel('攻角 α (°)', fontsize=12)
        axes[0].set_ylabel('升力系数 CL', fontsize=12)
        axes[0].set_title('升力系数随攻角变化', fontsize=13)
        axes[0].grid(True, alpha=0.3)

        # CD vs alpha
        axes[1].plot(df['alpha_deg'], df['CD'], 'r-', linewidth=2)
        axes[1].set_xlabel('攻角 α (°)', fontsize=12)
        axes[1].set_ylabel('阻力系数 CD', fontsize=12)
        axes[1].set_title('阻力系数随攻角变化', fontsize=13)
        axes[1].grid(True, alpha=0.3)

        # L/D vs alpha
        axes[2].plot(df['alpha_deg'], df['LD'], 'g-', linewidth=2)
        axes[2].set_xlabel('攻角 α (°)', fontsize=12)
        axes[2].set_ylabel('升阻比 L/D', fontsize=12)
        axes[2].set_title('升阻比随攻角变化', fontsize=13)
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'aero_CL_CD.png'), dpi=150)
        plt.close()
        print("  -> aero_CL_CD.png")

        # 阻力极曲线 (CD vs CL)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.plot(df['CL'], df['CD'], 'b-', linewidth=2)
        ax.set_xlabel('升力系数 CL', fontsize=12)
        ax.set_ylabel('阻力系数 CD', fontsize=12)
        ax.set_title('阻力极曲线 (Ma=0.72)', fontsize=13)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'aero_drag_polar.png'), dpi=150)
        plt.close()
        print("  -> aero_drag_polar.png")

    # 1.2 马赫数扫描
    fpath = os.path.join(OUTPUT_DIR, 'aero_mach_sweep.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))

        axes[0].plot(df['Ma'], df['CL_alpha'], 'b-o', linewidth=2, markersize=6)
        axes[0].set_xlabel('马赫数 Ma', fontsize=12)
        axes[0].set_ylabel('CLα (1/rad)', fontsize=12)
        axes[0].set_title('升力线斜率随Ma变化', fontsize=13)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(df['Ma'], df['CD0'], 'r-o', linewidth=2, markersize=6)
        axes[1].set_xlabel('马赫数 Ma', fontsize=12)
        axes[1].set_ylabel('CD0', fontsize=12)
        axes[1].set_title('零升阻力系数随Ma变化', fontsize=13)
        axes[1].grid(True, alpha=0.3)

        axes[2].plot(df['Ma'], df['K'], 'g-o', linewidth=2, markersize=6)
        axes[2].set_xlabel('马赫数 Ma', fontsize=12)
        axes[2].set_ylabel('K', fontsize=12)
        axes[2].set_title('诱导阻力因子随Ma变化', fontsize=13)
        axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'aero_mach_sweep.png'), dpi=150)
        plt.close()
        print("  -> aero_mach_sweep.png")
    else:
        print(f"  未找到 {fpath}")


# ============================================================================
# 2. 动力特性绘图
# ============================================================================
def plot_propulsion():
    ensure_dir(FIGURE_DIR)
    print("绘制动力特性图...")

    # 2.1 助推器推力曲线
    fpath = os.path.join(OUTPUT_DIR, 'propulsion_booster.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(df['time_s'], df['thrust_N'] / 1000, 'r-', linewidth=2)
        ax.set_xlabel('时间 (s)', fontsize=12)
        ax.set_ylabel('推力 (kN)', fontsize=12)
        ax.set_title('固体助推器推力-时间曲线', fontsize=13)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'propulsion_booster.png'), dpi=150)
        plt.close()
        print("  -> propulsion_booster.png")
    else:
        print(f"  未找到 {fpath}")

    # 2.2 涡扇推力随高度变化
    fpath = os.path.join(OUTPUT_DIR, 'propulsion_engine.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].plot(df['altitude_m'], df['thrust_N'], 'b-o', linewidth=2, markersize=6)
        axes[0].set_xlabel('高度 (m)', fontsize=12)
        axes[0].set_ylabel('推力 (N)', fontsize=12)
        axes[0].set_title('涡扇推力随高度变化 (Ma=0.72)', fontsize=13)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(df['altitude_m'], df['sfc_kg_N_h'], 'r-o', linewidth=2, markersize=6)
        axes[1].set_xlabel('高度 (m)', fontsize=12)
        axes[1].set_ylabel('耗油率 (kg/(N·h))', fontsize=12)
        axes[1].set_title('耗油率随高度变化', fontsize=13)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'propulsion_engine.png'), dpi=150)
        plt.close()
        print("  -> propulsion_engine.png")


# ============================================================================
# 3. 弹道曲线绘图
# ============================================================================
def plot_trajectory():
    ensure_dir(FIGURE_DIR)
    print("绘制弹道曲线...")

    traj_files = glob.glob(os.path.join(OUTPUT_DIR, 'trajectory_C*.csv'))
    if not traj_files:
        print("  未找到弹道数据文件")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    colors = ['b', 'r', 'g', 'm', 'c', 'orange', 'k']
    for idx, fpath in enumerate(sorted(traj_files)):
        df = pd.read_csv(fpath)
        case_name = os.path.basename(fpath).replace('trajectory_', '').replace('.csv', '')
        color = colors[idx % len(colors)]

        # 弹道曲线 x-y
        axes[0, 0].plot(df['x_m'] / 1000, df['y_m'], color=color, linewidth=1.5, label=case_name)

        # 速度-时间
        axes[0, 1].plot(df['time_s'] / 60, df['V_m_s'], color=color, linewidth=1.5, label=case_name)

        # 高度-时间
        axes[1, 0].plot(df['time_s'] / 60, df['y_m'], color=color, linewidth=1.5, label=case_name)

        # 质量-时间
        axes[1, 1].plot(df['time_s'] / 60, df['mass_kg'], color=color, linewidth=1.5, label=case_name)

    axes[0, 0].set_xlabel('水平距离 (km)', fontsize=12)
    axes[0, 0].set_ylabel('高度 (m)', fontsize=12)
    axes[0, 0].set_title('弹道曲线', fontsize=13)
    axes[0, 0].legend(fontsize=8, loc='best')
    axes[0, 0].grid(True, alpha=0.3)

    axes[0, 1].set_xlabel('时间 (min)', fontsize=12)
    axes[0, 1].set_ylabel('速度 (m/s)', fontsize=12)
    axes[0, 1].set_title('速度-时间曲线', fontsize=13)
    axes[0, 1].legend(fontsize=8, loc='best')
    axes[0, 1].grid(True, alpha=0.3)

    axes[1, 0].set_xlabel('时间 (min)', fontsize=12)
    axes[1, 0].set_ylabel('高度 (m)', fontsize=12)
    axes[1, 0].set_title('高度-时间曲线', fontsize=13)
    axes[1, 0].legend(fontsize=8, loc='best')
    axes[1, 0].grid(True, alpha=0.3)

    axes[1, 1].set_xlabel('时间 (min)', fontsize=12)
    axes[1, 1].set_ylabel('质量 (kg)', fontsize=12)
    axes[1, 1].set_title('质量-时间曲线', fontsize=13)
    axes[1, 1].legend(fontsize=8, loc='best')
    axes[1, 1].grid(True, alpha=0.3)

    plt.suptitle('BGM-109 巡航导弹弹道仿真结果对比', fontsize=15, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'trajectory_comparison.png'), dpi=150)
    plt.close()
    print("  -> trajectory_comparison.png")

    # 单独绘制基准工况的详细图
    baseline = os.path.join(OUTPUT_DIR, 'trajectory_C1_baseline_sea.csv')
    if os.path.exists(baseline):
        df = pd.read_csv(baseline)
        fig, axes = plt.subplots(3, 2, figsize=(14, 12))

        axes[0, 0].plot(df['x_m'] / 1000, df['y_m'], 'b-', linewidth=2)
        axes[0, 0].set_xlabel('水平距离 (km)'); axes[0, 0].set_ylabel('高度 (m)')
        axes[0, 0].set_title('弹道曲线 (基准工况)'); axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].plot(df['time_s'], df['V_m_s'], 'r-', linewidth=2)
        axes[0, 1].set_xlabel('时间 (s)'); axes[0, 1].set_ylabel('速度 (m/s)')
        axes[0, 1].set_title('速度-时间曲线'); axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].plot(df['time_s'], df['theta_deg'], 'g-', linewidth=2)
        axes[1, 0].set_xlabel('时间 (s)'); axes[1, 0].set_ylabel('弹道倾角 (°)')
        axes[1, 0].set_title('弹道倾角-时间曲线'); axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(df['time_s'], df['T_N'], 'm-', linewidth=2)
        axes[1, 1].set_xlabel('时间 (s)'); axes[1, 1].set_ylabel('推力 (N)')
        axes[1, 1].set_title('推力-时间曲线'); axes[1, 1].grid(True, alpha=0.3)

        axes[2, 0].plot(df['time_s'], df['L_N'], 'c-', linewidth=2, label='升力 L')
        axes[2, 0].plot(df['time_s'], df['D_N'], 'orange', linewidth=2, label='阻力 D')
        axes[2, 0].set_xlabel('时间 (s)'); axes[2, 0].set_ylabel('力 (N)')
        axes[2, 0].set_title('气动力-时间曲线'); axes[2, 0].legend(); axes[2, 0].grid(True, alpha=0.3)

        axes[2, 1].plot(df['time_s'], df['Ma'], 'b-', linewidth=2)
        axes[2, 1].set_xlabel('时间 (s)'); axes[2, 1].set_ylabel('马赫数 Ma')
        axes[2, 1].set_title('马赫数-时间曲线'); axes[2, 1].grid(True, alpha=0.3)

        plt.suptitle('BGM-109 基准工况 (海上巡航 Ma=0.72, H=15m)', fontsize=15, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'trajectory_baseline_detail.png'), dpi=150)
        plt.close()
        print("  -> trajectory_baseline_detail.png")


# ============================================================================
# 4. 参数扫描结果绘图
# ============================================================================
def plot_sweep():
    ensure_dir(FIGURE_DIR)
    print("绘制参数扫描结果...")

    # D-1: 速度扫描
    fpath = os.path.join(OUTPUT_DIR, 'sweep_velocity.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].plot(df['Ma'], df['range_km'], 'b-o', linewidth=2, markersize=8)
        axes[0].set_xlabel('巡航马赫数 Ma', fontsize=12)
        axes[0].set_ylabel('射程 (km)', fontsize=12)
        axes[0].set_title('射程随巡航速度变化', fontsize=13)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(df['Ma'], df['flight_time_min'], 'r-s', linewidth=2, markersize=8)
        axes[1].set_xlabel('巡航马赫数 Ma', fontsize=12)
        axes[1].set_ylabel('飞行时间 (min)', fontsize=12)
        axes[1].set_title('飞行时间随巡航速度变化', fontsize=13)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'sweep_velocity.png'), dpi=150)
        plt.close()
        print("  -> sweep_velocity.png")

    # D-2: 高度扫描
    fpath = os.path.join(OUTPUT_DIR, 'sweep_altitude.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].plot(df['altitude_m'], df['range_km'], 'b-o', linewidth=2, markersize=8)
        axes[0].set_xlabel('巡航高度 (m)', fontsize=12)
        axes[0].set_ylabel('射程 (km)', fontsize=12)
        axes[0].set_title('射程随巡航高度变化', fontsize=13)
        axes[0].grid(True, alpha=0.3)

        axes[1].plot(df['altitude_m'], df['flight_time_min'], 'r-s', linewidth=2, markersize=8)
        axes[1].set_xlabel('巡航高度 (m)', fontsize=12)
        axes[1].set_ylabel('飞行时间 (min)', fontsize=12)
        axes[1].set_title('飞行时间随巡航高度变化', fontsize=13)
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'sweep_altitude.png'), dpi=150)
        plt.close()
        print("  -> sweep_altitude.png")

    # D-3: 质量扫描
    fpath = os.path.join(OUTPUT_DIR, 'sweep_mass.csv')
    if os.path.exists(fpath):
        df = pd.read_csv(fpath)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(df['mass_kg'], df['range_km'], 'g-o', linewidth=2, markersize=8)
        ax.set_xlabel('初始质量 (kg)', fontsize=12)
        ax.set_ylabel('射程 (km)', fontsize=12)
        ax.set_title('射程随初始质量变化', fontsize=13)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(FIGURE_DIR, 'sweep_mass.png'), dpi=150)
        plt.close()
        print("  -> sweep_mass.png")


# ============================================================================
# 主函数
# ============================================================================
def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'

    print("=" * 60)
    print("  BGM-109 巡航导弹弹道仿真 -- 后处理绘图")
    print("=" * 60)
    print(f"数据目录: {OUTPUT_DIR}")
    print(f"图表目录: {FIGURE_DIR}\n")

    if mode in ('all', 'aero'):
        plot_aero()
    if mode in ('all', 'propulsion'):
        plot_propulsion()
    if mode in ('all', 'trajectory'):
        plot_trajectory()
    if mode in ('all', 'sweep'):
        plot_sweep()

    print(f"\n全部绘图完成! 图表保存在: {FIGURE_DIR}/")


if __name__ == '__main__':
    main()
