---
name: oommf-wer-sim
description: 自动化 SOT/STT-MRAM 器件写入错误率 WER 仿真。用 OOMMF 热演化器（默认 Xf_ThermHeunEvolve 纯热 Heun 积分 + Oxs_ScriptUZeeman 电流等效场，已验证热统计正确；可选 Xf_ThermSpinXferEvolve 原生 STT）注入 300K 热噪声，实现 REE（Rare Event Enhancement / importance splitting，Fine et al. arXiv:1603.08512）重要性分裂：分段 boxsi 仿真 + 按 m_proj 分裂滞后轨迹 + 权重记账，以 ~10²~10³ 条轨迹把 WER 推到朴素蒙特卡洛无法到达的 10⁻⁹ 量级。当用户需要计算写入错误率、WER(J)/WER(t) 曲线、REE 稀有事件增强仿真、热噪声随机 LLG 统计、SOT-MRAM/STT-MRAM 写入可靠性评估时使用；也支持朴素 MC 交叉验证。
---

# OOMMF WER 仿真（oommf-wer-sim）

## 快速开始

```bash
# 1) 单组合 WER（默认 REE）
python scripts/wer_sim.py --config my_device.yaml --out results

# 2) 电流扫描 → WER(J) 曲线（对数间隔采样）
python scripts/wer_sim.py --config my_device.yaml \
    --sweep "J=4e11,6e11,8e11,1.2e12" --out results

# 3) 朴素 MC 交叉验证（验证 REE 逻辑，建议同时跑同一配置 method=mc）
python scripts/wer_sim.py --config my_device_mc.yaml --out results_mc

# 4) 只查看计划
python scripts/wer_sim.py --config my_device.yaml --dry-run

# 5) 只对已有 results 做数据处理
python scripts/process_wer.py --out results
```

## 环境

- OOMMF: `C:\oommf\oommf2.1\oommf\oommf\oommf.tcl`（需已编译热扩展：
  `Xf_ThermHeunEvolve`（默认，推荐）与 `Xf_ThermSpinXferEvolve`（可选）；
  用下述最小 MIF 冒烟测试验证）
- tclsh: `C:\ActiveTcl\bin\tclsh86t.exe`
- Python: numpy / matplotlib / PyYAML（matplotlib 仅 process_wer.py 出图需要）

### 扩展可用性冒烟测试

跑不通就先用最小 MIF 验证热扩展（注意 boxsi 需 `-tk 0 -exitondone 1`，
否则完成后 Tk 窗口常驻、进程不退出）：

```bash
tclsh oommf.tcl boxsi -tk 0 -exitondone 1 test.mif
# 出现 ODT 且含 "Xf_ThermHeunEvolve:evolve:Temperature" 列 => 扩展可用
```

### 热引擎选择（重要）

| evolver | 热统计 | 说明 |
|---|---|---|
| `thermheun`（T>0 默认） | ✅ 已验证正确 | Xf_ThermHeunEvolve 纯热 Heun 积分。50nm 宏自旋 Δ=9.66 实测平衡 ⟨mz⟩=0.945 vs 理论 0.948（0.3% 误差），无虚假翻转。**SOT/STT 力矩用 Oxs_ScriptUZeeman 等效场承载**（常数 m̂₀ 近似）。 |
| `thermspin`（宏自旋 ❌ / 真实多格 ✅） | 分场景 | Xf_ThermSpinXferEvolve 原生 STT 力矩 + 热。**50nm 宏自旋 T>0 出现虚假翻转（0.1ns 内翻 -z）❌**；但 **SARR 90° 真实多格器件（130nm、4225 格、逐区 mp/J + DMI、T=300）热统计正常（J=0 自检 mz=-0.906 零翻转），且零温动力学与原模型 Oxs_SpinXferEvolve 逐位一致（3ns 末态 +0.132224）✅**。结论：**宏自旋/单格场景用 thermheun 等效场；真实多格原生 STT 器件可直接用 thermspin 跑热 WER**。 |

`evolver: auto`（默认）：T>0 → thermheun；T=0 且 mode 含 stt → thermspin（零温
Jc 扫描用原生 STT）；T=0 其余 → thermheun。真实多格原生 STT 器件（带逐区 mp/J、
DMI 等）热 WER 应显式选 `thermspin`（见上表，宏自旋退化场景除外）。

**每次更换配置后先跑自检**（J=0、T、短时长）：
`⟨mz⟩_eq ≈ 1 - 1/(2Δ)`（Δ=K_eff·V/(k_B T)，误差 <5% 为合格；出现 mz→-1 或 ⟨mz⟩ 显著
偏低即为热统计异常，停止并检查配置）。

## 配置字段（YAML / JSON）

几何/材料/力矩字段与 `oommf-batch-sim` 一致，额外增加热学与 WER 字段：

| 分类 | 字段 | 默认 | 说明 |
|---|---|---|---|
| 几何 | `length width thickness cell` | 40e-9 40e-9 1e-9 5e-9 | 器件尺寸与网格（m）；**网格须 ≤ 交换长度一半** |
| | `shape image colormap` | box | box / ellipse / image（位图定义 MTJ 形状） |
| 材料 | `Ms A Ku axis alpha m0` | 1100e3 13e-12 0.85e4 [1,0,0] 0.08 [0.996,0.087,0] | `axis`=易轴，也是 WER 预测坐标投影方向（PMA 用 [0,0,1]） |
| | `use_exchange use_demag` | true true | 宏自旋验证可关（等效各向异性吸收退磁场） |
| 力矩 | `mode` | stt | stt / sot / both / none（SOT 为 ScriptUZeeman 常数 m̂₀ 近似） |
| | `J mp P Lambda` | 4e11 [-1,0,0] 0.4 2.0 | STT 电流密度 / 极化方向 / 极化率 / Lambda |
| | `theta_SH xi_FL current_dir` | 0.3 0.1 [1,0,0] | SOT 参数 |
| | `t_on t_off` | 1e-9 100e-9 | 电流脉冲窗口；**t_on 之前为无电流热化段**，WER 在 t_off 报告 |
| 热学 | `T` | 300.0 | 温度 K（必须 >0，WER 依赖热噪声） |
| | `evolver` | auto | auto / thermheun（推荐）/ thermspin（仅对照，见上） |
| | `evolver_opts` | {} | 透传 evolver 选项，用户键覆盖默认。**T>0 默认注入作者黄金配置** `{min_timestep: 1e-14, max_timestep: 1e-14, energy_precision: -1}`（固定 10fs 步长 + 关能量精度检查，热统计可靠且比默认 1fs 快 10 倍）；显式给 `min_timestep` 可再提速（噪声按该值标定，>5e-14 需先自检） |
| WER | `method` | ree | ree（默认，文献方法）/ mc（朴素蒙特卡洛验证） |
| | `n_traj` | 200 | ree: 初始轨迹数；mc: 总轨迹数 |
| | `checkpoint_ns` | 10.0 | 检查点间隔 (ns)，建议脉冲时长的 1/5~1/10 |
| | `split_rule` | laggard | **laggard**（分裂 m_proj>mz_th 的滞后轨迹，文献验证）/ frontier（分裂最接近翻转的轨迹，PPT 变体） |
| | `mz_th` | 0.95 | 分裂阈值；电流越大越接近 1（2Ic0→0.95, 3Ic0→0.97） |
| | `split_topk` | 0 | frontier 用：每检查点最多分裂条数，0=全部候选 |
| | `split_factor` | 2 | 每条分裂子轨迹数（文献 2，PPT 2-4） |
| | `max_traj` | 500 | 活跃轨迹上限（超出则按优先级裁剪分裂） |
| | `switched_mz` | 0.0 | m_proj 低于此值判已切换（文献 0） |
| | `seed` | 1 | 随机种子基准（`random_seed: false` 时每段派生唯一种子） |
| | `random_seed` | true | **True=每条轨迹独立随机 seed（默认，推荐）**；False=按 `seed` 确定性派生 |
| 仿真 | `DeltaStage` | 0.2e-9 | ODT 记录步长（段被切成 时长/DeltaStage 个 stage） |
| | `keep_ovf` | true | 是否保留全部段 OMF（false 只留每轨迹末段重启点，省磁盘） |
| | `workers threads_per_boxsi` | 4 2 | 并行 boxsi 进程数 / 单进程线程数（避免 CPU 超订） |
| | `timeout_s basename` | 3600 wer | 单段超时 / 输出前缀 |

### 宏自旋验证（推荐先跑）

```yaml
# 宏自旋单格：1 cell，关交换退磁，等效各向异性吸收退磁场
length: 50e-9, width: 50e-9, thickness: 1e-9, cell: 50e-9   # 或 4nm（Ku 相应放大）
A: 0                    # 关交换（use_exchange: false）
shape: box              # 单格即宏自旋（不要 image）
Ms: 1100e3
Ku: 1.6e4               # 50nm：Δ=K_eff·V/(k_B·T)≈9.66（4nm 用 2.5e6）
axis: [0, 0, 1]         # PMA
m0: [0, 0, 1]           # 写入 0 态（+z）；勿用面内方向——那是势垒顶点
mode: sot, J: 8e11, theta_SH: 0.3, xi_FL: 0.1, current_dir: [1,0,0]
T: 300
evolver: auto           # -> thermheun（热统计可靠）
method: ree, n_traj: 200, checkpoint_ns: 10, mz_th: 0.95, max_traj: 500
t_on: 1e-9, t_off: 3e-9     # 1ns 热化 + 2ns 脉冲（宏自旋验证用小窗口）
```

**先自检**：同一配置 `J: 0, t_on: 0, t_off: 5e-9` 跑 5ns，⟨mz⟩≈0.948 为合格。

### 文献参数组（Δ=60，已验证 ✅）

与 PPT 冒烟组（Δ≈9.66）同框架，改用文献 arXiv:1603.08512 的材料参数：

```yaml
length: 50e-9, width: 50e-9, thickness: 1e-9, cell: 50e-9
Ms: 1.1e6, Ku: 8e5, axis: [0,0,1], m0: [0,0,1]
use_exchange: false
use_demag: true        # 关键！保留退磁，让单格退磁因子给出 K_eff≈1.0e5 → Δ≈60
                       # （K_eff = Ku - ½μ₀Ms²(N_zz-N_xx) ≈ 8e5 - 7.0e5 = 1.0e5）
mode: sot, J: 4e12, theta_SH: 0.3, xi_FL: 0.1
T: 300, evolver: auto  # -> thermheun
method: mc（或 ree）, n_traj: 4-8, checkpoint_ns: 1, mz_th: 0.95
t_on: 1e-9, t_off: 3e-9   # 1ns 热化 + 2ns 脉冲
```

本机实测结果（50nm 宏自旋，全通过）：
- **Δ=60 精确成立**：J=0 自检 5ns 末态 ⟨mz⟩=0.9913 vs 理论 0.9917（0.04% 误差）
- **J 扫描单调**：J=2/4/8e12 → WER=0.75/0.25/0.0（MC n=4），Jc0≈4e12
- **REE vs MC 一致**：J=4e12 → REE 0.375 vs MC 0.25（n=4-8 统计噪声内）
- 注意：Δ=60 时 Jc0 比 Δ≈9.66 组高约一个量级（K_eff 比例 1.0e5/1.6e4≈6.3），
  电流扫描要用 2e12-1e13 量级；t_off 建议 ≤5ns（Δ 大时翻转更快）

## 工作流（REE）

```
配置 → 时间网格 [0, ck, ..., t_off]
  → 初始 n_traj 条轨迹（权重 1/n_traj）
  → 逐检查点：
      并行跑本段 boxsi（热演化器：默认 Xf_ThermHeunEvolve + ScriptUZeeman
                     电流等效场；10fs 固定步长热统计经自检；段末 OMF，
                     RandomSeed 每段唯一）
      读 ODT 末行 m_proj → 判切换（m_proj<switched_mz → 计权入切换桶）
      分裂（laggard: m_proj>mz_th 的滞后轨迹 ×split_factor，权重均分；
           超出 max_traj 按最滞后优先裁剪）
  → 汇总：WER(t) = 1 - 累计切换权重；WER_final 在 t_off
  → process_wer.py：WER(t)/WER(J) 图、轨迹叠加、population、report.md
```

### 结果归档

```
results/
├── <combo>/                    # 每个参数组合独立目录（如 J_4e11）
│   ├── summary.json            # WER_final、WER(t) 轨迹、权重分辨率
│   ├── wer_trace.csv
│   ├── trajectories.json       # 每条轨迹: 权重/父代/切换时刻/mz(t)
│   └── trajectories/traj_xxxxx/# 段 MIF / ODT / OMF / boxsi 日志
├── batch_config.yaml
├── summary_all.csv / summary_all.json
├── wer_vs_time.png             # WER(t)（log-y，主图）
├── wer_vs_<param>.png          # 扫描参数依赖（如 WER vs J）
├── trajectories.png / population.png
└── report.md
```

## 数据处理与判读

- `wer_final`：脉冲结束时刻 WER（主结果）。
- `min_active_weight`：末态最小单条权重 = 当前统计量可分辨的 WER 下限。
- MC 可测下限 ≈ 1/n_traj；REE 与 MC 在 WER ≥ 1/n_mc 区间必须一致。
- WER(t) 尾部噪声增大 / 台阶伪影：调大 mz_th、减小检查点间隔、或加大 n_traj。

## 与 oommf-batch-sim / oommf-sot-jc 的分工

| 需求 | 用哪个 |
|---|---|
| WER / WER(J) / 写入可靠性 | **本技能** |
| 任意参数批量零温扫描、Jc 粗估 | `oommf-batch-sim` |
| 单器件 Jc 二分精修 | `oommf-sot-jc` |

## 注意事项

1. **算力现实**：热演化器 T>0 时强制固定步长 = `min_timestep`（作者设计，误差检查
   被旁路）。默认 10fs → 50nm 宏自旋 ~30-60 s/ns、真实多格器件数分钟/ns；
   100ns 脉冲 × 200 轨迹 × 每 ns 1 分钟 ≈ 数天（workers 并行后仍以小时~天计）。
   先用宏自旋 + 短脉冲（1-10ns）验证方法，再评估真实器件算力预算。
2. **热引擎**：宏自旋/单格默认 thermheun（热统计已验证）；`thermspin` 的 T>0 热噪声
   在 **50nm 宏自旋单格**实测异常（0.1ns 虚假翻转），但 **真实多格原生 STT 器件
   （SARR 90°：4225 格、逐区 mp/J + DMI）热统计正常**（J=0 自检 mz=-0.906 零翻转，
   零温 3ns 末态与原模型逐位一致 +0.132224）。选型见"热引擎选择"表。
9. **WER 时长必须覆盖翻转窗口（截断伪影）**：真实器件临界电流下翻转慢（SARR 90°
   零温 10.5ns 才完成写入，1-3ns 面内徘徊、3ns 后弛豫翻转）。若仿真时长 < 翻转时间，
   WER 恒等于 1.0 是**截断伪影**而非物理结论。判定规则：
   - 先跑零温长时 ODT（原模型常有）确定翻转时间 t_sw；
   - WER 仿真总时长 ≥ t_sw（含弛豫）；6ns 对 SARR Jc 组得 WER=0.75（n=4）；
   - 若只关心脉冲末态，须在报告明确"写入含弛豫、按末态判定"。
3. **SOT 定量局限**：默认 SOT 用常数 m̂₀ 近似（趋势研究足够）；定量 SOT-WER
   需自定义 SOT 扩展（随瞬时 m 更新 σ×m）——本技能的分段重启 + REE 记账
   框架可无缝替换为任意自定义扩展的 evolver。
4. **热化**：t_on 之前为无电流热化段；初始 m0 取 [0,0,1]（平衡附近）即可，
   宏自旋从平衡态出发，脉冲本身主导写入。
5. **阈值**：mz_th 默认 0.95；电流越大越接近 1。出现台阶伪影即调大。
6. **失败处理**：单段失败写 FAILED.txt 并跳过该轨迹，其余继续；汇总含状态。
7. **本机 RNG 多流现象（重要）**：OOMMF 2.1a1 的 `Oc_Random` 全局状态在本机实测
   存在"多流交替"——同一 MIF（SHA256 全同、seed 写死）多次运行可能得到 2~3 种
   固定随机流（如宏自旋冒烟中 -0.167/-0.353/+0.441 三终态交替），串行/并行/独立
   TMP/独立 `-logfile` 均无法消除。其根因在 OOMMF 层（非技能层可修）。**影响与对策**：
   - 单条轨迹的"同 seed 复现"不可保证 → 技能默认 `random_seed: true` 随机化，
     保证各轨迹流独立，**WER 统计无偏**（流本身的噪声分布是物理正确的）；
   - REE 重要性分裂不依赖同 seed 复现（子轨迹本就应独立演化），逻辑不受影响；
   - 如需逐条核对，`trajectories/traj_*/<basename>_*.mif` 记录了 intended seed，
     实际流可能不同。
8. **方法细节**：完整算法、阈值选取、宏自旋验证协议见
   `references/ree_method.md`。
