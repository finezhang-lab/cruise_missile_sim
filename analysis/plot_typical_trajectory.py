#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_typical_trajectory.py -- 典型多阶段战术巡航弹道绘图程序 (精确校准版)
输出:
1. trajectory_typical_altitude.png: 典型弹道全剖面高度轨迹图 (含分段标注与局部放大)
2. trajectory_typical_velocity.png: 典型弹道全剖面速度、马赫数与推力阻力自适应节流特性图
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.font_manager as fm

# 配置中文字体
_CJK_FONT_NAMES = [
    'Microsoft YaHei',
    'SimHei',
    'Source Han Sans SC',
    'Noto Sans CJK SC',
    'WenQuanYi Micro Hei',
    'DejaVu Sans',
]

def _find_cjk_font():
    available = {f.name for f in fm.fontManager.ttflist}
    for name in _CJK_FONT_NAMES:
        if name in available:
            return name
    return 'DejaVu Sans'

cjk_font = _find_cjk_font()
matplotlib.rcParams['font.sans-serif'] = [cjk_font, 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['font.size'] = 11

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.dirname(_SCRIPT_DIR)
OUTPUT_DIR = os.path.join(_PROJECT_ROOT, 'output')
FIGURE_DIR = os.path.join(_SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

csv_path = os.path.join(OUTPUT_DIR, 'trajectory_C0_typical_5k_to_15m.csv')
if not os.path.exists(csv_path):
    print(f"错误: 未找到数据文件 {csv_path}，请先运行仿真生成!")
    exit(1)

df = pd.read_csv(csv_path)

# 转换为常用工程单位
x_km = df['x_m'] / 1000.0
y_m = df['y_m']
t_min = df['time_s'] / 60.0
V_ms = df['V_m_s']
Ma = df['Ma']
T_kN = df['T_N'] / 1000.0
D_kN = df['D_N'] / 1000.0
mass_kg = df['mass_kg']

# 阶段分界点提取与统计
phase_cuts = []
for p in [0, 1, 2, 3, 4]:
    sub = df[df['phase'] == p]
    if not sub.empty:
        phase_cuts.append({
            'phase': p,
            'x_start': sub['x_m'].min() / 1000.0,
            'x_end': sub['x_m'].max() / 1000.0,
            't_start': sub['time_s'].min() / 60.0,
            't_end': sub['time_s'].max() / 60.0,
            'y_max': sub['y_m'].max(),
            'y_min': sub['y_m'].min(),
        })

print(f"加载弹道数据成功: 共 {len(df)} 记录点，总航程 {x_km.iloc[-1]:.3f} km，飞行时间 {t_min.iloc[-1]:.2f} min")

# =============================================================================
# 图 1: 典型弹道全剖面高度轨迹图 (Altitude Profile)
# =============================================================================
fig = plt.figure(figsize=(13, 8), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.8], hspace=0.32, wspace=0.22)

# 主图: 全航程高度曲线
ax_main = fig.add_subplot(gs[0, :])
ax_main.plot(x_km, y_m, color='#b65338', linewidth=2.4, label='典型弹道高度轨迹 (5km爬升平飞 $\\rightarrow$ 俯冲 $\\rightarrow$ 15m掠海平飞)')

# 阶段底色与纵向标注
phase_colors = ['#f5e6cf', '#e2eef7', '#eaf3ee', '#fcf0ed', '#edf5f0']
phase_labels = [
    '① 助推加速段\n(0~0.23km)',
    '② 高空爬升段\n(0.23~93.14km)',
    '③ 5km高空平飞巡航段 (H = 5000 m, 93.35~999.87km)\n[自适应节流巡航, 净增航程 906.52km]',
    '④ 受控俯冲下滑段\n(1000.10~1101.45km, 净增 101.35km)',
    '⑤ 15m低空掠海平飞段 (H = 15 m, 1101.69~2648.39km)\n[海平面突防, 净增 1546.70km, 燃油耗尽]'
]

for i, pc in enumerate(phase_cuts):
    ax_main.axvspan(pc['x_start'], pc['x_end'], color=phase_colors[i], alpha=0.45)
    if i < len(phase_cuts) - 1:
        ax_main.axvline(pc['x_end'], color='#8c3b24', linestyle='--', linewidth=1.0, alpha=0.6)

# 阶段标签文本
y_txt_positions = [1800, 3200, 5250, 3200, 1800]
x_txt_positions = [
    15,
    45,
    (phase_cuts[2]['x_start'] + phase_cuts[2]['x_end']) / 2,
    (phase_cuts[3]['x_start'] + phase_cuts[3]['x_end']) / 2,
    (phase_cuts[4]['x_start'] + phase_cuts[4]['x_end']) / 2
]

for i, txt in enumerate(phase_labels):
    ax_main.text(x_txt_positions[i], y_txt_positions[i], txt,
                 ha='center', va='center', fontsize=9.2, fontweight='bold',
                 color='#1f2724', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#ded9d0', alpha=0.9))

# 标注关键控制点
ax_main.scatter([x_km.iloc[0], phase_cuts[1]['x_end'], phase_cuts[2]['x_end'], phase_cuts[3]['x_end'], x_km.iloc[-1]],
                [y_m.iloc[0], 5000.0, 5000.0, 15.0, y_m.iloc[-1]],
                color='#8c3b24', s=45, zorder=5)

ax_main.annotate('发射出筒 (y=15m, θ=30°)',
                 xy=(x_km.iloc[0], y_m.iloc[0]), xytext=(70, 700),
                 arrowprops=dict(arrowstyle='->', color='#2b3934', lw=1.2),
                 fontsize=9, color='#1f2724')

ax_main.annotate('平稳捕获5000m高空平飞\n(超调量 0.35m, y_max=5000.35m)',
                 xy=(phase_cuts[1]['x_end'], 5000.35), xytext=(phase_cuts[1]['x_end'] + 40, 4400),
                 arrowprops=dict(arrowstyle='->', color='#2b3934', lw=1.2),
                 fontsize=9, color='#1f2724')

ax_main.annotate('转入受控俯冲 (下滑倾角约 -3.0°)',
                 xy=(phase_cuts[2]['x_end'], 4999.0), xytext=(phase_cuts[2]['x_end'] - 200, 4200),
                 arrowprops=dict(arrowstyle='->', color='#2b3934', lw=1.2),
                 fontsize=9, color='#1f2724')

ax_main.annotate('拉平切入海平面 (目标 15.0m)\n(安全下冲余量 0.10m, 最低点 15.10m)',
                 xy=(phase_cuts[3]['x_end'], 15.1), xytext=(phase_cuts[3]['x_end'] + 80, 800),
                 arrowprops=dict(arrowstyle='->', color='#2b3934', lw=1.2),
                 fontsize=9, color='#1f2724')

ax_main.set_xlim(-20, 2700)
ax_main.set_ylim(-150, 5800)
ax_main.set_ylabel('飞行高度 y (m)', fontsize=12, fontweight='bold', color='#1f2724')
ax_main.set_title('BGM-109 典型多阶段全剖面弹道高度轨迹剖面图 (数值仿真计算结果)', fontsize=14, fontweight='bold', pad=12, color='#1f2724')
ax_main.grid(True, linestyle=':', alpha=0.6, color='#c5bfb4')
ax_main.legend(loc='lower left', frameon=True, facecolor='white', edgecolor='#ded9d0')

# 子图 1: 起飞助推与爬升段局部放大 (0 ~ 120 km)
ax_sub1 = fig.add_subplot(gs[1, 0])
mask_climb = x_km <= 120.0
ax_sub1.plot(x_km[mask_climb], y_m[mask_climb], color='#b65338', linewidth=2.0)
ax_sub1.axhline(5000.0, color='#2c597a', linestyle='--', linewidth=1.0, label='5000m 目标高度')
ax_sub1.axvline(phase_cuts[0]['x_end'], color='#8f5c1d', linestyle=':', label=f'助推分离 (t=4.0s, x={phase_cuts[0]["x_end"]:.2f}km)')
ax_sub1.set_xlim(-2, 120)
ax_sub1.set_ylim(0, 5400)
ax_sub1.set_xlabel('水平航程 x (km)', fontsize=11, color='#1f2724')
ax_sub1.set_ylabel('高度 y (m)', fontsize=11, color='#1f2724')
ax_sub1.set_title('局部放大 A: 出架助推、爬升与5km平飞捕获', fontsize=11, fontweight='bold', color='#1f2724')
ax_sub1.grid(True, linestyle=':', alpha=0.6)
ax_sub1.legend(loc='lower right', fontsize=8.5, frameon=True, facecolor='white')

# 子图 2: 俯冲下滑与掠海平飞拉平局部放大 (950 ~ 1200 km)
ax_sub2 = fig.add_subplot(gs[1, 1])
mask_dive = (x_km >= 980.0) & (x_km <= 1180.0)
ax_sub2.plot(x_km[mask_dive], y_m[mask_dive], color='#b65338', linewidth=2.0)
ax_sub2.axhline(15.0, color='#256a4e', linestyle='--', linewidth=1.0, label='15m 掠海目标高度')
ax_sub2.axvline(phase_cuts[2]['x_end'], color='#8c3b24', linestyle=':', label=f'俯冲起点 (x={phase_cuts[2]["x_end"]:.1f}km)')
ax_sub2.axvline(phase_cuts[3]['x_end'], color='#256a4e', linestyle=':', label=f'掠海拉平点 (x={phase_cuts[3]["x_end"]:.1f}km)')
ax_sub2.set_xlim(980, 1180)
ax_sub2.set_ylim(-50, 5400)
ax_sub2.set_xlabel('水平航程 x (km)', fontsize=11, color='#1f2724')
ax_sub2.set_ylabel('高度 y (m)', fontsize=11, color='#1f2724')
ax_sub2.set_title('局部放大 B: 5km高空转受控俯冲至15m掠海拉平', fontsize=11, fontweight='bold', color='#1f2724')
ax_sub2.grid(True, linestyle=':', alpha=0.6)
ax_sub2.legend(loc='upper right', fontsize=8.5, frameon=True, facecolor='white')

alt_fig_path = os.path.join(FIGURE_DIR, 'trajectory_typical_altitude.png')
plt.savefig(alt_fig_path, dpi=300, bbox_inches='tight', facecolor='#ffffff')
plt.close(fig)
print(f"高度图已生成: {alt_fig_path}")

# =============================================================================
# 图 2: 典型弹道全剖面速度、马赫数与推力阻力自适应节流特性图 (Velocity Profile)
# =============================================================================
fig2, (ax_v, ax_t) = plt.subplots(2, 1, figsize=(13, 8), dpi=300, sharex=True)

# 速度与马赫数曲线
ax_v.plot(x_km, V_ms, color='#2c597a', linewidth=2.2, label='飞行速度 V (m/s)')
ax_v.set_ylabel('速度 V (m/s)', fontsize=12, fontweight='bold', color='#2c597a')
ax_v.tick_params(axis='y', labelcolor='#2c597a')
ax_v.set_ylim(0, 280)
ax_v.grid(True, linestyle=':', alpha=0.6)

# 双 Y 轴: 马赫数
ax_ma = ax_v.twinx()
ax_ma.plot(x_km, Ma, color='#b65338', linewidth=1.8, linestyle='-', label='飞行马赫数 Ma')
ax_ma.axhline(0.72, color='#8c3b24', linestyle='--', linewidth=1.2, alpha=0.7, label='目标巡航马赫数 Ma=0.720')
ax_ma.set_ylabel('马赫数 Ma', fontsize=12, fontweight='bold', color='#b65338')
ax_ma.tick_params(axis='y', labelcolor='#b65338')
ax_ma.set_ylim(0.0, 0.90)

# 标注巡航状态点物理参数
ax_ma.annotate('高空巡航闭环: Ma 严格收敛至 0.720\n(V ≈ 230.7 m/s, 当地音速 320.5 m/s)',
               xy=(500, 0.72), xytext=(350, 0.82),
               arrowprops=dict(arrowstyle='->', color='#b65338', lw=1.2),
               fontsize=9, color='#1f2724', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ded9d0'))

ax_ma.annotate('掠海平飞闭环: Ma 严格收敛至 0.720\n(V ≈ 245.0 m/s, 海平面音速 340.2 m/s)',
               xy=(1800, 0.72), xytext=(1600, 0.58),
               arrowprops=dict(arrowstyle='->', color='#b65338', lw=1.2),
               fontsize=9, color='#1f2724', bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#ded9d0'))

ax_v.set_title('BGM-109 典型弹道飞行速度、马赫数与推力阻力自适应节流特性曲线', fontsize=14, fontweight='bold', pad=12, color='#1f2724')

# 合并图例
lines_v, labels_v = ax_v.get_legend_handles_labels()
lines_m, labels_m = ax_ma.get_legend_handles_labels()
ax_v.legend(lines_v + lines_m, labels_v + labels_m, loc='lower right', frameon=True, facecolor='white')

# 下图: 推力与阻力推阻自适应曲线
ax_t.plot(x_km, T_kN * 1000.0, color='#8f5c1d', linewidth=2.0, label='发动机推力 T (N)')
ax_t.plot(x_km, D_kN * 1000.0, color='#256a4e', linewidth=1.8, linestyle='--', label='气动力阻力 D (N)')

for i, pc in enumerate(phase_cuts):
    ax_t.axvspan(pc['x_start'], pc['x_end'], color=phase_colors[i], alpha=0.35)
    if i < len(phase_cuts) - 1:
        ax_t.axvline(pc['x_end'], color='#8c3b24', linestyle='--', linewidth=1.0, alpha=0.5)

# 标注推阻平衡与节流状态
ax_t.annotate('爬升主发全推力加速\n(初始 2692 N, 5km降至 1947 N, 均值 2344 N)',
              xy=(40, 2400), xytext=(80, 2900),
              arrowprops=dict(arrowstyle='->', color='#8f5c1d', lw=1.0),
              fontsize=8.5, color='#1f2724')

ax_t.annotate('高空自适应节流平飞\n(初始加速1946 N, 稳态节流至1486 N, 均值1518 N)',
              xy=(550, 1500), xytext=(450, 2400),
              arrowprops=dict(arrowstyle='->', color='#8f5c1d', lw=1.0),
              fontsize=8.5, color='#1f2724')

ax_t.annotate('俯冲初始怠速节流至 926 N\n(海面拉平加大至 2359 N以克服海面高阻)',
              xy=(1050, 950), xytext=(920, 1800),
              arrowprops=dict(arrowstyle='->', color='#8f5c1d', lw=1.0),
              fontsize=8.5, color='#1f2724')

ax_t.annotate('海平面掠海推阻平衡\n(T ≈ D ≈ 2414 N, 节流自适应 2365~2448 N)',
              xy=(1800, 2400), xytext=(1700, 1600),
              arrowprops=dict(arrowstyle='->', color='#8f5c1d', lw=1.0),
              fontsize=8.5, color='#1f2724')

ax_t.set_xlim(-20, 2700)
ax_t.set_ylim(0, 3500)
ax_t.set_xlabel('水平航程 x (km)', fontsize=12, fontweight='bold', color='#1f2724')
ax_t.set_ylabel('推力 / 阻力 (N)', fontsize=12, fontweight='bold', color='#1f2724')
ax_t.grid(True, linestyle=':', alpha=0.6)
ax_t.legend(loc='upper right', frameon=True, facecolor='white')

vel_fig_path = os.path.join(FIGURE_DIR, 'trajectory_typical_velocity.png')
plt.savefig(vel_fig_path, dpi=300, bbox_inches='tight', facecolor='#ffffff')
plt.close(fig2)
print(f"速度图已生成: {vel_fig_path}")

print("全部绘图成功完成!")
