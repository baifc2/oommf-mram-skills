# OOMMF 批量仿真工作流详解

## 1. 总体流程

```
基础配置 (yaml/json)
   │  --sweep / --runs / CLI
   ▼
组合展开（笛卡尔积 或 显式列表）
   │
   ▼
每个组合：build_mif() 生成 run.mif
   │
   ▼
boxsi 运行（独立文件夹，支持并行 workers）
   │
   ▼
process_results.py：ODT 解析 → 指标 → 汇总 → 出图 → Jc
```

## 2. 配置 Schema 全表

### 几何 / 网格

| 字段 | 默认 | 说明 |
|---|---|---|
| `length` | 40e-9 | x 方向尺寸 (m) |
| `width` | 40e-9 | y 方向尺寸 (m) |
| `thickness` | 1e-9 | z 方向厚度 (m) |
| `cell` | 5e-9 | 面内网格 (m)，z 向固定一层 |
| `shape` | box | `box` / `ellipse` / `image` |
| `image` | "" | image 模式的位图绝对路径 |
| `colormap` | {red: MTJ, white: vacuum} | 位图颜色 → 区域名映射 |

### 磁性材料

| 字段 | 默认 | 说明 |
|---|---|---|
| `Ms` | 1100e3 | 饱和磁化 (A/m) |
| `A` | 13e-12 | 交换常数 (J/m) |
| `Ku` | 0.85e4 | 单轴各向异性常数 (J/m³) |
| `axis` | [1,0,0] | 易轴方向（同时决定切换判据投影轴） |
| `alpha` | 0.08 | 阻尼系数 |
| `m0` | [0.996,0.087,0] | 初始磁化方向 |

### 电流 / 力矩

| 字段 | 默认 | 说明 |
|---|---|---|
| `mode` | stt | none / field / stt / sot / both |
| `J` | 4e11 | STT 电流密度或 SOT 电荷电流密度 (A/m²) |
| `mp` | [-1,0,0] | STT 自旋极化方向 |
| `P` | 0.4 | STT 极化率（OOMMF 默认 0.4） |
| `Lambda` | 2.0 | STT Lambda（**OOMMF 默认 2，不是 1**；设 1 会去掉 η 对 m·mp 的依赖、弱化力矩） |
| `t_on` | 1e-9 | 电流开启时刻 (s) |
| `t_off` | 4e-9 | 电流关闭时刻 (s) |
| `theta_SH` | 0.3 | SOT 自旋霍尔角 |
| `xi_FL` | 0.1 | SOT 类场/类阻尼比 |
| `current_dir` | [1,0,0] | SOT 电荷电流方向 |
| `H_bias` | [0,0,0] | 恒定外磁场 (A/m) |

### 仿真控制

| 字段 | 默认 | 说明 |
|---|---|---|
| `stage_count` | 50 | TimeDriver 阶段数 |
| `DeltaStage` | 0.2e-9 | 每阶段时长 (s)，总时长 = 二者乘积 |
| `T` | 0.0 | 温度；仅 0 支持（零温） |
| `basename` | run | OOMMF 输出前缀（产生 run.odt） |
| `schedule_ovf` | false | 是否输出磁化归档 .ovf |
| `timeout_s` | 3600 | 单 run 超时 |
| `workers` | 1 | 并行 boxsi 进程数（>1 时用线程池） |

### 扫描 / 列表

```yaml
# sweep：笛卡尔积
sweep:
  J: [2e11, 4e11, 6e11]
  Ms: [1e6, 1.1e6]
  # 列表型参数也支持：
  # m0: [[0.996,0.087,0], [1,0,0]]
  # mp: [[-1,0,0], [1,0,0]]

# runs：显式列表（优先于 sweep）
runs:
  - name: baseline
    params: {J: 4e11}
  - name: high_ku
    params: {Ku: 2e4, J: 6e11}
```

CLI 等价写法：`--sweep "J=2e11,4e11,6e11" --sweep "Ms=1e6,1.1e6"`。

## 3. MIF 生成逻辑

每个组合的 `run.mif` 由 `build_mif()` 生成，结构：

1. **Atlas**：`Oxs_BoxAtlas` / `Oxs_EllipseAtlas`（axis 半径 = L/2, W/2）/ `Oxs_ImageAtlas`（位图 + colormap）
2. **Mesh**：`Oxs_RectangularMesh`，z 向单层
3. **能量项**：`Oxs_UniaxialAnisotropy`（image 模式用 `Oxs_AtlasScalarField` 按区域给 Ku，vacuum 为 0）、`Oxs_UniformExchange`、`Oxs_Demag`
4. **外场**：`H_bias` 非零时加 `Oxs_UZeeman:bias`（恒场，`Hrange` 的 steps=0）
5. **SOT 场**：`mode` 含 sot 时加 `Oxs_ScriptUZeeman:sot`，脚本按 `t_on/t_off` 窗口返回场+时间导数（smoothstep 包络，避免 evolver 步长崩溃）
6. **Evolver**：STT 用 `Oxs_SpinXferEvolve`（含 `J_profile` 窗口脚本）；否则 `Oxs_RungeKuttaEvolve`
7. **Driver**：`Oxs_TimeDriver`，`stopping_time=DeltaStage`、`stage_count`；image 模式 Ms 用 AtlasScalarField
8. **Tcl 脚本**：`Jprofile_script`（STT 电流窗口）、`SOT_field_script`（SOT 场窗口）写在 MIF 末尾

### SOT 等效场（常数近似）

```
Hmag = ℏ θ_SH J / (2 e μ₀ Ms t)          # 类阻尼场幅度
H_DL = Hmag * (σ̂ × m̂₀)                    # 类阻尼项，m 用初始方向近似
H_FL = xi_FL * Hmag * σ̂                   # 类场项
```

局限：真实 SOT 的 H_DL 随瞬时 m 变化（σ̂ × m(t)），脚本方案用 m̂₀ 恒定方向近似；
做 J 扫描趋势、阈值粗估足够；定量结果需 OOMMF 自定义 SOT 扩展。

## 4. 数据处理细节

### ODT 解析

- `# Columns:` 行：花括号包裹含空格的列名，用正则 `\{[^}]*\}|\S+` 切分后去括号
- 数据行：空白分隔，与列一一对应
- 关键列匹配：`Simulation time`、`::mx/::my/::mz`、`Total energy`、`average J`、`::Stage`

### 指标定义

| 指标 | 定义 |
|---|---|
| `m_proj_initial/final` | 沿易轴投影磁化的初/末值（`m·axiŝ`） |
| `m_proj_min/max` | 全程投影磁化极值 |
| `m_pulse_min/max` | 脉冲窗口 [t_on,t_off] 内的投影磁化极值 |
| `switched` | `m_initial × m_final < 0`（跨零） |
| `E_min_J/E_final_J` | 总能量极值/末值 |
| `J_applied_max` | 脉冲期间最大电流 |

### Jc 粗提取

对 `J` 扫描：按 J 排序，最高未切换 J (J_lo) 与最低切换 J (J_hi) 之间的几何均值
作为 Jc 估计，报告区间 `[J_lo, J_hi]`。精修请用 `oommf-sot-jc` 技能的二分法。

### 出图

- `trajectory_overlay.png`：全部 run 的投影磁化-时间曲线叠加
- `dependence_<param>.png`：单参数扫描（末态投影磁化 vs 参数；J 自动对数轴）
- `heatmap_<p1>_<p2>.png`：双参数末态热图；`_switched` 后缀为切换标记图

## 5. 与 oommf-sot-jc 的分工

| 需求 | 用哪个 |
|---|---|
| 单个器件精修 Jc（二分法，~0.1% 精度） | `oommf-sot-jc` |
| 任意参数批量扫描、结果归档、数据处理 | `oommf-batch-sim`（本技能） |
| SOT 无场翻转（STT 辅助）扫描 | 本技能 `mode: both` |

## 6. 常见问题

1. **数值被当字符串**：PyYAML 把 `2e11`（无符号指数）解析为字符串，脚本已内置 `_coerce`
   递归纠正；手写配置时也可用 `2e+11` 规避
2. **image 模式不切换**：先确认位图颜色与 colormap 精确一致（JPEG 压缩会糊边缘），
   再检查网格是否足够细（用户项目 1nm 网格可切换、2nm 不切换）
3. **boxsi 卡住**：检查是否误用 GUI 版 `tclsh.exe`；脚本默认用 `tclsh86t.exe`
4. **超时**：`timeout_s` 调大，或减小网格/阶段数；`workers` 可并行加速
5. **多 odt 文件**：一个 run 可能产生多个 .odt（多 Destination），脚本取最大者为主表
