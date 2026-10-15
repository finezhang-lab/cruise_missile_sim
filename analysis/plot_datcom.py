# -*- coding: utf-8 -*-
# 本脚本由 aerodynamics.cpp 自动生成: 绘制 DATCOM 数据插值得到的 CL/CD 随攻角曲线
import os
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

_simhei = r'C:\Windows\Fonts\simhei.ttf'
if os.path.exists(_simhei):
    fm.fontManager.addfont(_simhei)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'SimSun']
plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_FILE = os.path.join(HERE, 'aero_datcom.csv')
PNG_FILE = os.path.join(HERE, 'aero_datcom.png')

alpha, CL, CD, LD = [], [], [], []
with open(CSV_FILE, encoding='utf-8', newline='') as f:
    for row in csv.DictReader(f):
        alpha.append(float(row['alpha_deg']))
        CL.append(float(row['CL']))
        CD.append(float(row['CD']))
        LD.append(float(row['LD']))

# DATCOM 数据节点 (攻角每隔 2°)
nodes = [-4, -2, 0, 2, 4, 6, 8, 10, 12, 14, 16]
node_CL = [CL[alpha.index(a)] for a in nodes]
node_CD = [CD[alpha.index(a)] for a in nodes]

fig, ax1 = plt.subplots(figsize=(8.5, 5.5))
ax1.plot(alpha, CL, 'b-', linewidth=2, label='插值曲线 CL')
ax1.plot(nodes, node_CL, 'ro', markersize=7, label='DATCOM 数据节点')
ax1.set_xlabel('攻角 α (°)', fontsize=12)
ax1.set_ylabel('升力系数 CL', color='b', fontsize=12)
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True, linestyle='--', alpha=0.5)

ax2 = ax1.twinx()
ax2.plot(alpha, CD, 'r-', linewidth=2, label='插值曲线 CD')
ax2.plot(nodes, node_CD, 'bs', markersize=7, label='DATCOM 数据节点')
ax2.set_ylabel('阻力系数 CD', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=10)

ax1.set_title('DATCOM 数据插值: 升力系数 CL / 阻力系数 CD 随攻角变化', fontsize=13)
plt.savefig(PNG_FILE, dpi=150, bbox_inches='tight')
plt.close(fig)
print('已生成 DATCOM 插值关系图:', PNG_FILE)
