---
name: oommf-batch-sim
description: 批量 OOMMF 微磁仿真：从基础配置一键生成多个参数组合的 MIF 项目，逐个运行 boxsi 仿真，结果按参数自动命名归档到独立文件夹，并自动完成 .odt 数据处理（轨迹/切换判定/汇总/出图/Jc 提取）。当用户需要对 OOMMF 项目批量修改参数仿真（扫描电流密度 J、材料参数 Ms/Ku/A/alpha、几何尺寸、STT/SOT 参数、脉冲窗口等）、批量运行多个 MIF 项目、自动归档仿真结果并处理输出数据时使用。
---

# OOMMF 批量微磁仿真（oommf-batch-sim）

## 环境

- OOMMF: `C:\oommf\oommf2.1\oommf\oommf\oommf.tcl`
- tclsh: `C:\ActiveTcl\bin\tclsh86t.exe`（找不到时回退 `tclsh.exe` / PATH 上的 `tclsh`）
- Python: 需 `numpy`、`matplotlib`、`PyYAML`
- 不依赖 Ubermag：直接生成 MIF + 调 `boxsi`，与现有 `oommf-sot-jc` 技能互补（该技能做单器件 Jc 精修，本技能做任意参数批量扫描）

## 快速开始

```bash
# 1) 单次仿真（用配置里的默认参数）
python scripts/batch_sim.py --config my_device.yaml --out results

# 2) 参数扫描（笛卡尔积；--sweep 可多次）
python scripts/batch_sim.py --config my_device.yaml \
    --sweep "J=4e11,6e11,8e11" --sweep "Ms=1e6,1.1e6" --out results

# 3) 显式运行列表（每个 run 精确指定参数和名字）
python scripts/batch_sim.py --config my_device.yaml --runs runs.yaml --out results

# 4) 只查看运行计划，不执行
python scripts/batch_sim.py --config my_device.yaml --dry-run

# 5) 只对已有 results 目录做数据处理
python scripts/process_results.py --out results
```

## 配置字段（YAML / JSON）

基础配置即"项目参数"，任意字段都可被扫描覆盖：

| 分类 | 字段 | 说明 |
|---|---|---|
| 几何 | `length` `width` `thickness` | 器件尺寸 (m) |
| | `cell` | 网格尺寸 (m)，须 ≤ 交换长度一半（CoFeB 约 2–3 nm） |
| | `shape` | `box` \| `ellipse` \| `image`（image 用位图定义 MTJ 形状） |
| | `image` `colormap` | image 模式的位图路径 + 颜色→区域映射，如 `{red: MTJ, white: vacuum}` |
| 材料 | `Ms` `A` `Ku` `axis` `alpha` | 饱和磁化、交换、各向异性、易轴、阻尼 |
| 初态 | `m0` | 初始磁化方向，如 `[0.996, 0.087, 0]` |
| 力矩 | `mode` | `none` \| `field` \| `stt` \| `sot` \| `both` |
| | `J` | STT 电流密度或 SOT 电荷电流密度 (A/m²) |
| | `mp` `P` `Lambda` | STT 自旋极化方向、极化率（默认 0.4）、Lambda（**默认 2.0，与 OOMMF 默认一致**；设为 1 会弱化力矩） |
| | `theta_SH` `xi_FL` `current_dir` | SOT 自旋霍尔角、类场比、电流方向 |
| | `t_on` `t_off` | 电流脉冲窗口 (s) |
| | `H_bias` | 恒定外磁场 (A/m) |
| 仿真 | `stage_count` `DeltaStage` | 阶段数 × 每阶段时长（总时长 = 二者乘积） |
| | `basename` | OOMMF 输出文件名前缀 |
| | `schedule_ovf` | 是否同时输出磁化归档 .ovf |
| | `timeout_s` `workers` | 单仿真超时、并行进程数 |
| 扫描 | `sweep` | `{参数: [值列表]}`，笛卡尔积 |
| | `runs` | `[{name: 名字, params: {参数: 值}}]`，显式列表（优先于 sweep） |

## 结果命名与归档

每个参数组合一个独立文件夹，自动命名：`参数名_值` 用 `__` 连接，如：

```
results/
├── J_4e11__Ms_1.1e6/
│   ├── run.mif          # 该组合生成的 MIF 项目
│   ├── params.yaml      # 该组合完整参数记录
│   ├── run.odt          # OOMMF DataTable 输出
│   └── results.json     # 处理后的指标
├── batch_config.yaml    # 本次批次配置快照
├── run_plan.csv         # 运行计划（人可读）
├── summary.csv / summary.json   # 全部 run 汇总
├── trajectory_overlay.png       # 全部轨迹叠加
├── dependence_J.png             # 单参数扫描依赖图
├── heatmap_p1_p2.png            # 双参数扫描热图
└── report.md                    # 自动生成的报告
```

## 数据处理（自动执行）

batch_sim.py 跑完后自动调用 `process_results.py`：

1. **解析 .odt**：读取 `# Columns:` 表头 + 数据行，提取时间、mx/my/mz、Stage、电流、能量
2. **指标**：沿易轴的投影磁化初末值/极值、脉冲窗口极值、切换判定（投影磁化跨零）、能量极值
3. **汇总**：summary.csv / summary.json（每 run 一行，含全部参数与指标）
4. **出图**：轨迹叠加图；单参数扫描→依赖图（J 自动对数轴）；双参数→末态热图 + 切换标记图
5. **Jc 提取**：扫描参数含 `J` 且有切换时，报告 Jc 区间（最高未切换 J 与最低切换 J 的几何均值）

## 力矩模式说明

- `stt`：`Oxs_SpinXferEvolve`，`J_profile` Tcl 脚本控制脉冲窗口（与 STT-MTJ 单层项目同款写法）
- `sot`：`Oxs_ScriptUZeeman` 时变等效场（脚本返回场+时间导数，smoothstep 包络保证可微）；H_DL = ℏθ_SH·J/(2eμ₀Ms·t)·(σ×m̂₀)，H_FL = ξ·H_DL·σ，m 依赖用初始方向 m̂₀ 近似（零温趋势研究足够，高精度需自定义 SOT 扩展）
- `both`：STT + SOT 同时注入（无场翻转的对称破缺方案）
- `field`：仅恒定外磁场 `H_bias`
- `none`：纯弛豫

## 注意事项

1. **网格**：cell 必须 ≤ 交换长度一半，否则丢失空间模式、切换阈值失真（本项目 1nm 网格下 J=4e11 可切换，2nm 网格不切换——网格必须与物理一致）
2. **image 模式**：位图颜色必须与 colormap 精确匹配；JPEG 压缩可能模糊边缘，建议用 PNG 纯色位图
3. **温度**：`T` 仅支持 0（零温 RungeKutta）；有限温需 `xf_thermspinxferevolve` 热扩展，本技能未内置
4. **切换判据**：默认取沿易轴投影磁化是否跨零（PMA 器件看 mz，面内器件看 mx，自动按 `axis` 选择）
5. **失败处理**：单 run 失败会写 `FAILED.txt` 并继续其他 run，汇总中标记状态
