# -*- coding: utf-8 -*-
# 本脚本由 aerodynamics.cpp 自动生成, 用于绘制升力系数 CL 与阻力系数 CD 随攻角变化关系图
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# 注册中文字体, 避免图中汉字显示为方框
_simhei = r'C:\Windows\Fonts\simhei.ttf'
if os.path.exists(_simhei):
    fm.fontManager.addfont(_simhei)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(HERE, 'aero_CL_CD_zh.csv')
PNG_FILE = os.path.join(HERE, 'aero_CL_CD_zh.png')

alpha, CL, CD, LD = [], [], [], []
with open(CSV_FILE, encoding='utf-8-sig', newline='') as f:
    reader = csv.reader(f)
    in_data = False
    for row in reader:
        if not row:
            continue
        if row[0].strip().startswith('攻角'):
            in_data = True
            continue
        if not in_data:
            continue
        try:
            a = float(row[0]); cl = float(row[1]); cd = float(row[2]); ld = float(row[3])
        except (ValueError, IndexError):
            continue
        alpha.append(a); CL.append(cl); CD.append(cd); LD.append(ld)

fig, ax1 = plt.subplots(figsize=(8.5, 5.5))
ax1.plot(alpha, CL, 'b-o', linewidth=2, markersize=5, label='升力系数 CL')
ax1.set_xlabel('攻角 α (°)', fontsize=12)
ax1.set_ylabel('升力系数 CL', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, linestyle='--', alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(alpha, CD, 'r-s', linewidth=2, markersize=5, label='阻力系数 CD')
ax2.set_ylabel('阻力系数 CD', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=11)

ax1.set_title('升力系数 CL / 阻力系数 CD 随攻角 α 的变化关系', fontsize=13)
plt.savefig(PNG_FILE, dpi=150, bbox_inches='tight')
plt.close(fig)
print('已生成关系图:', PNG_FILE)
print('攻角范围 %.1f° ~ %.1f°, 共 %d 个采样点' % (alpha[0], alpha[-1], len(alpha)))

# ---- 图2: 升阻比 L/D 随攻角变化 ----
LD_PNG = os.path.join(HERE, 'aero_LD_zh.png')
fig2, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(alpha, LD, 'g-o', linewidth=2, markersize=5, label='升阻比 L/D')
imax = LD.index(max(LD))
ax.plot(alpha[imax], LD[imax], 'r*', markersize=14, label='最大升阻比')
ax.annotate('最大 L/D = %.2f @ α = %.1f°' % (LD[imax], alpha[imax]),
            xy=(alpha[imax], LD[imax]),
            xytext=(alpha[imax] - 6.0, LD[imax] - 3.0),
            arrowprops=dict(arrowstyle='->', color='r'),
            fontsize=11, color='r')
ax.axhline(0, color='k', linewidth=0.8)
ax.set_xlabel('攻角 α (°)', fontsize=12)
ax.set_ylabel('升阻比 L/D', fontsize=12)
ax.set_title('升阻比 L/D 随攻角 α 的变化曲线', fontsize=13)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(loc='lower right', fontsize=11)
plt.savefig(LD_PNG, dpi=150, bbox_inches='tight')
plt.close(fig2)
print('已生成升阻比关系图:', LD_PNG)
print('最大升阻比 L/D = %.3f, 出现在攻角 α = %.1f°' % (LD[imax], alpha[imax]))
