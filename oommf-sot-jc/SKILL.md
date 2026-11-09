---
name: oommf-sot-jc
description: 自动化提取 SOT-MRAM 临界切换电流密度 Jc。基于 OOMMF 微磁学仿真（通过 Ubermag/oommfc Python DSL 或直接 MIF+boxsi），支持电流密度参数扫描、二分法精修、脉冲宽度依赖曲线、STT 辅助无场翻转等模式。当用户需要计算 SOT-MRAM 的 Jc、扫描写入电流、提取翻转阈值、或做 SOT 器件参数优化时使用。
---

# OOMMF SOT-MRAM Jc 自动提取

## 环境

- OOMMF: `C:\oommf\oommf2.1\oommf\oommf\oommf.tcl`（2.1a1，含 `xf_thermspinxferevolve` 等扩展）
- tclsh: `C:/ActiveTcl/bin/tclsh.exe`
- Ubermag（可选）: `pip install ubermag`，自动生成 MIF + 调用 OOMMF + 读 OVF

## 核心脚本

`scripts/extract_jc.py` — 一站式 Jc 提取。

```bash
# 快速演示（默认器件参数）
python scripts/extract_jc.py --quick --out jc_demo

# 自定义器件
python scripts/extract_jc.py --config my_device.yaml --out results/

# 只扫一段电流范围
python scripts/extract_jc.py --config dev.yaml --J-min 5e11 --J-max 2e12 --n-points 12
```

输出：
- `jc_scan.csv` — J, switched, mz_final
- `jc_scan_results.json` — Jc, 不确定度, 完整扫描数据
- `jc_scan.png` — 切换映射图

## 工作流

### 1. 定义器件

配置文件（JSON 或 YAML），关键字段：

```yaml
length: 60e-9      # x (m)
width: 60e-9       # y (m)
thickness: 1.2e-9  # z (m)
cell: 3e-9         # mesh cell (m)
Ms: 1.1e6          # A/m
A: 20e-12          # exchange (J/m)
Ku: 1.2e6          # anisotropy (J/m3)
axis: [0, 0, 1]    # PMA
alpha: 0.02
theta_SH: 0.3      # spin Hall angle
xi_FL: 0.1         # field-like/damping-like ratio
current_dir: [1, 0, 0]  # charge current direction
t_pulse: 5e-9      # pulse width (s)
J_min: 1e11        # scan range (A/m2)
J_max: 3e12
n_points: 15
```

### 2. SOT 实现方式

脚本自动选择后端：

- **Ubermag 可用** → `oommfc` 生成 MIF，`Slonczewski` 项做 STT，`Zeeman` 时变场做 SOT 等效场
- **Ubermag 不可用** → 直接写 MIF + `tclsh oommf.tcl boxsi`，SOT 用 `Oxs_UZeeman` 等效场

SOT 等效场（每个时间步根据当前 m 更新）：

```
H_DL = (ℏ θ_SH J)/(2 e μ0 Ms t) * (σ × m)    # damping-like
H_FL = ξ_FL * H_DL_mag * σ                    # field-like
```

### 3. 无场翻转

设置 `use_stt_assist: true`，脚本会同时注入：
- SOT 面内电流（主翻转力矩）
- STT 垂直电流（破缺对称，`J_stt = J_stt_ratio * J_sot`）

### 4. Jc 提取

1. **粗扫**：对数间隔 10–20 个 J 值，每个跑 relax → pulse → settle
2. **判定**：`mz_final < 0` 视为翻转
3. **二分精修**：在阈值附近做 5 次几何中点迭代，精度 ~0.1%
4. **输出**：Jc ± 不确定度

### 5. 脉冲宽度依赖（进阶）

循环多个 `t_pulse`，提取 Jc(τ)，拟合热激活模型：

```
Jc(τ) = Jc0 * (1 + τ0/τ)
```

外推到 τ→∞ 得零温 Jc0。在脚本外用循环调用，或参考 `references/sot_physics.md`。

## 关键注意事项

1. **cell 尺寸**：必须 ≤ 交换长度的一半（~2–3 nm for CoFeB），否则空间模式丢失
2. **SOT 等效场近似**：纯 UZeeman 方案中 H_DL 的 m 依赖用平均 mz 近似；高精度结果应使用自定义 SOT 扩展或每子步更新场
3. **热化**：有限温 Jc 需要 `xf_thermspinxferevolve` 原生热噪声；零温 Jc 用 `Oxs_RungeKuttaEvolve`
4. **退磁场**：纳米点器件中退磁场不可忽略，必须包含 `Oxs_Demag`
5. **Jc 公式校验**：跑仿真前先用解析公式估算（见 references），确认扫描范围合理

## 参考

- SOT 物理、Jc 公式、OOMMF 实现选项：`references/sot_physics.md`
- Ubermag 论文：arXiv:2105.08355
- Ubermag 文档：https://ubermag.github.io
