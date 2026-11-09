# REE 方法详解：SOT/STT-MRAM 写入错误率（WER）仿真

## 1. 物理背景与 WER 定义

写入错误率（Write Error Rate, WER）= 在给定电流脉冲（幅值 J、时长 t_p）结束后，
器件仍未完成磁化翻转的概率：

```
WER(t_p) = P( 位元在脉冲结束时刻未切换 )
```

热涨落使翻转时间服从宽分布：同一电流下不同轨迹（热噪声实现）的翻转时刻不同，
因此在有限脉冲内存在"未写完"的概率。目标 WER 通常 ≤ 10⁻⁹（有 ECC）甚至 10⁻¹⁹。

朴素蒙特卡洛（MC）要分辨 WER≈10⁻⁹ 需要 ≥10⁹ 条独立随机 LLGS 轨迹，不可行。
REE（Rare Event Enhancement，即重要性分裂 importance splitting）通过**放大稀有事件
的采样密度、同时按权重还原概率**，用 ~10³ 条在行轨迹即可推到 10⁻⁹ 量级。

**主要文献**
- J. Fine, A. Hovorka, P. D'Aquino, C. Serpico, M. d'Aquino, G. Rowlands, ...
  "Write error rate of spin-transfer-torque random access memory including
  micromagnetic effects using rare event enhancement", IEEE Trans. Magn.
  (2016) / arXiv:1603.08512。本技能实现的默认算法。
- 重要性采样变体：arXiv:1901.03890（biased Monte Carlo / IS 求切换与非切换概率）。
- 本技能 OOMMF 后端（默认）：Xf_ThermHeunEvolve（K. Fong, NUS 纯热扩展，
  Heun 积分 + 朗之万热噪声）+ Oxs_ScriptUZeeman 电流等效场；
  可选 Xf_ThermSpinXferEvolve（原生 STT，**50nm 宏自旋 T>0 热统计异常；真实多格原生 STT 器件热统计正常**，见 §4）。

## 2. REE 算法（arXiv:1603.08512）

### 2.1 预测坐标（progress coordinate / predictor）

用沿目标轴（PMA 为 z）的归一化投影磁化：

```
m_proj(t) = ⟨m(t)⟩·axiŝ / |⟨m(t)⟩|
```

m_proj 越大 → 面内分量越小 → STT 越弱 → 翻转越慢 → 越可能成为**稀有未切换事件**。
同时越接近 1 时热噪声相对越重要、轨迹越"随机"，分裂才有意义
（完全确定性的轨迹分裂无增益）。

### 2.2 分裂规则（本技能默认 laggard，与文献一致）

仿真被切成检查点段 [t_l, t_{l+1}]。每个检查点：

| m_proj 区间 | 处理 |
|---|---|
| m_proj < 0（已切换） | 判为写入成功，该轨迹权重计入"切换桶"，不再演化 |
| m_proj > m_z,th（滞后，接近 +1） | **分裂**为 M 条子轨迹，每条权重 = 父权重 / M |
| m_z,th ≥ m_proj ≥ 0 | 原样继续，权重不变 |

- 权重全程守恒：任何分裂只是把一条权重 w 换成 M 条 w/M。
- 子轨迹从父轨迹当前磁化状态（段末 OMF）继续，各带独立热噪声实现
  （独立 RandomSeed），互不相关。
- WER(t) = 1 - 到 t 为止的累计切换权重。

### 2.3 阈值 m_z,th 的选择

文献的宏自旋研究（图 2）：m_z,th = 0.95 在 0.7/0.95/0.99 三者中波动最小、最稳。
经验规律：**电流越大，最优 m_z,th 越接近 1**（2I_c0→0.95, 3I_c0→0.97, 4I_c0→0.99）。

系统化选法（图 3）：对若干阈值计算"跃迁时间"（从 m_proj=m_z,th 落到 0 的时间）
的样本标准差，取标准差-阈值曲线的拐点。拐点左侧路径随机性强、方差大，右侧接近
确定性，拐点即随机→确定的过渡，对应最优分裂阈值。

实操：先用 m_z,th=0.95 起步；若 WER(t) 曲线出现明显台阶状伪影，向 1 调大阈值。

### 2.4 检查点间隔与轨迹预算

- 文献用**自适应间隔**：当"已切换数 ≈ 滞后数"时终止当前段，使在行轨迹数守恒。
- 本技能默认**固定间隔**（checkpoint_ns，PPT 验证用 10 ns），用 max_traj 上限兜底：
  候选按优先级（laggard 最滞后优先 / frontier 最接近翻转优先）排序后，
  只分裂到不超过 max_traj 为止。
- 建议：检查点间隔取脉冲时长的 1/5~1/10；max_traj ≥ 2×n_traj。

### 2.5 与 PPT 变体（frontier 规则）的区别

本技能同时实现 PPT 汇报中描述的变体（`split_rule: frontier`：分裂**最接近翻转**
的 K 条）。注意：**文献验证到 10⁻⁹ 的是 laggard 规则**；frontier 规则放大的是
"接近翻转"的常见分支（适合短时/低电流、WER 偏高的验证区间），推低 WER 的能力
未经文献验证。两者在 WER 不低的区间（如 10⁻¹ 量级）应给出相同结果，可用于交叉验证。

## 3. 宏自旋验证协议（先验证方法，再上真实器件）

文献与 PPT 都用宏自旋单格先验证 REE 与 MC 的一致性：

1. **器件**：1×1×1 单格（如 4×4×1 nm 一个 cell），配置
   `use_exchange: false, use_demag: false`（关闭交换与退磁），
   用等效各向异性吸收退磁场：K_eff = Ku - ½μ₀Ms²(N_zz - N_xx)，
   薄膜退磁因子 N_xx=N_yy≈0.027、N_zz≈0.946（按具体长宽比算）。
   若保留 `use_demag: true`，则 Oxs_Demag 对单格给出该格的自然退磁因子，
   Ku 直接按目标 Δ 设定即可（Δ 公式见下）。
2. **热稳定因子**：Δ = V/(k_B T)·[Ku - ½μ₀Ms²(N_zz-N_xx)]，
   文献 Δ=60；PPT 验证用 Δ=9.7（Ku=2.5e6）。
3. **电流**：取 J = (2~4)×J_c0 的过驱动区，WER 随时间指数下降；
   也取 J≈0.2×J_c0 的低电流区验证重叠（PPT 验证：2-6 ns 区间 REE 与 MC 的
   WER=0.6/0.2/0.2 完全一致）。
4. **热化**：初始分布应热化——从偏离平衡的 m0（文献从 mz=0 落到势垒顶）无电流跑
   t_on（如 5 ns）到脉冲开始。
5. **初始态方向**：PMA 器件写入"0"态必须从 mz≈+1 出发（如 m0=[0,0,1]）。
   **若误用面内方向（mz=0），等于把器件放在势垒顶点**：无电流时热噪声会把它
   随机推向 ±z，造成虚假切换、WER 恒等于 ~0.5——这是最常见的配置错误。
5. **判据**：REE 与 MC 在 WER ∈ [1/n_mc, 0.5] 区间逐点一致；
   MC 的可测下限 ≈ 1/n_mc，低于该值只信 REE。

## 4. OOMMF 后端实现要点（本机已验证）

| 项 | 用法 |
|---|---|
| 热演化器 | 默认 `Specify Xf_ThermHeunEvolve:evolve { alpha; temperature 300; min_timestep 1e-14; max_timestep 1e-14; energy_precision -1 }`（纯热 Heun 积分，热统计已验证正确）；可选 `Xf_ThermSpinXferEvolve`（原生 STT 力矩，**50nm 宏自旋 T>0 热统计异常；真实多格原生 STT 器件（SARR 90°，4225 格）热统计正常且零温与原模型逐位一致，可直接用于热 WER**，见下方验证结论） |
| 电流力矩 | SOT/STT 一律用 `Oxs_ScriptUZeeman` 等效场承载（SOT: H=Hmag(σ̂×m̂₀)+ξHmagσ̂；STT: H=-Hmag(mp̂×m̂₀)，Hmag=ℏ·P·J/(2eμ₀Ms·t)），smoothstep 脉冲包络，与热噪声叠加 |
| 随机种子 | MIF 顶层 `RandomSeed <int>`（同时初始化 Tcl 与 C++ 随机源）；默认 `random_seed: true` 每条轨迹用独立随机 seed（子轨迹独立演化）。**本机实测 OOMMF 2.1a1 的 Oc_Random 全局状态存在"多流交替"**：同一 MIF 多次运行可能得到 2~3 种固定随机流（同 seed 不可复现），串行/并行/独立 TMP/独立 -logfile 均无法消除——根因在 OOMMF 层（Oc_Srand/Oc_RandomState 初始化时序，非技能层可修）。对策：随机化 seed 保证各轨迹流独立 → WER 统计无偏；REE 分裂不依赖同 seed 复现（子轨迹本就应独立演化），逻辑不受影响。详见 SKILL.md 注意事项 7 |
| 分段重启 | 段末 `Schedule Oxs_TimeDriver::Magnetization archive stage N` 写 OMF（命名 `<basename>-Oxs_TimeDriver-Magnetization-<stage>-<iter>.omf`）；下一段 `Specify Oxs_FileVectorField:m0f { atlas :atlas; file <omf> }` + `m0 :m0f` |
| 脉冲窗口锚定 | 每段 MIF 内嵌全局起始时刻 t_start，SOT/STT 场脚本用 `gt = t_start + total_time` 判断是否落在 [t_on, t_off] |
| 进程管理 | boxsi 用 Popen + 文件重定向 + 超时 kill（capture_output 会挂起）；**必须加 `-tk 0 -exitondone 1`**：`-tk 0` 关闭 Tk 窗口（否则 MIF 报错弹对话框永久等待），`-exitondone 1` 让 boxsi 求解完成即退出（否则窗口常驻、进程不结束被误判超时） |
| 线程 | `-threads <n>` 限制单进程线程数，避免 workers×默认 8 线程超订 CPU |
| 并行 | ThreadPoolExecutor 多进程 boxsi；每进程独立工作目录 |

### 热引擎验证结论（关键，本机实测）

两个 Xf 扩展共享 rk2heun 热积分框架，热统计表现分场景：

| 测试 | 结果 |
|---|---|
| `Xf_ThermHeunEvolve`（50nm 宏自旋，Δ=9.66，T=300，J=0，5ns） | **⟨mz⟩=0.945 vs 理论 0.948（0.3% 误差），零翻转** ✅ |
| `Xf_ThermSpinXferEvolve`（同 50nm 宏自旋配置） | 0.1ns 内翻至 -z 并停留（虚假翻转，T=0 对照完美冻结）❌ |
| `Xf_ThermSpinXferEvolve`（**SARR 90° 真实多格器件**：130nm、4225 格、逐区 mp/J + DMI，T=300） | **J=0 自检 mz=-0.906 稳定零翻转（热统计正常）✅**；**零温 3ns 末态 +0.132224 与原始模型 Oxs_SpinXferEvolve 逐位一致 ✅** |

结论：
- **宏自旋/单格场景**（退化）：T>0 用 thermheun 等效场路线；thermspin 只在 T=0 用。
- **真实多格原生 STT 器件**（带逐区 mp/J、DMI、形状）：thermspin 可直接跑热 WER
  （忠实原生 STT，免等效场近似）——先跑 J=0 热自检（≥2ns 无虚假翻转）+ 零温对照
  （与原模型 ODT 同刻逐位一致）验证后使用。
- **每次改配置后跑自检**：J=0、T、短时长（≥2τ），⟨mz⟩ 与 `1-1/(2Δ)` 偏差 <5% 为合格。

### 时间步长与算力

热演化器 **T>0 时强制固定步长 = min_timestep**（作者设计：误差检查头被旁路，
`allowed_absolute_step_error=PI` 等），**自适应/大步长提速在 T>0 不可行**。
- 默认 min_timestep=1e-15（1fs）过慢；**作者 example1 推荐 1e-14（10fs）**
  （`min_timestep 1e-14; max_timestep 1e-14; energy_precision -1`），本技能
  T>0 默认注入该配置。
- 实测 50nm 宏自旋 10fs 步长 ~20-40 s/ns（单进程）；4nm 单格 ~5-10 倍更慢；
  真实多格器件数分钟/ns。
- 100ns 脉冲 × 200 轨迹 ≈ 每轨迹 40-80 分钟 × 200 ≈ 数天（workers 并行后
  以小时~天计）。**先用宏自旋 + 短脉冲（1-10ns）验证方法**，真实器件 WER
  量产前先评估算力预算（或改用 GPU 微磁/解析模型）。
- 显式给 `evolver_opts.min_timestep` 可提速（噪声按该值标定，>5e-14 需先自检）。

### SOT 模式的局限

OOMMF 无原生 SOT 项。本技能按 oommf-batch-sim 的常数 m̂₀ 近似实现
（H_DL = ℏθ_SH J/(2eμ₀Ms·t)·(σ̂×m̂₀)，H_FL = ξ·H_DL·σ̂，smoothstep 包络注入
ScriptUZeeman，与热噪声叠加）。定量 SOT-WER 需要随瞬时 m 更新的自定义 SOT 扩展，
或采用"弯曲电流"空间渐变 STT 等效 SOT（Oxs_SpinXferEvolve 原生支持，
见 PPT 器件实例与 oommf-sot-jc 技能）。本技能的分段重启 + REE 记账框架可
无缝替换为任意自定义扩展的 evolver（改 MIF 生成器即可）。

## 5. 结果解读

- `wer_final`：脉冲结束时刻 WER（主结果）。
- `min_active_weight`：末态最小单条权重，指示当前统计量能分辨的 WER 下限。
- WER(t) 曲线尾部噪声增大是正常现象（样本变稀疏），报告时给出
  min_active_weight 作为置信度说明。
- 若 WER(t) 出现明显台阶/伪影：调大 m_z,th（更接近 1）或减小检查点间隔。
- 与 MC 的交叉验证：在 WER ≥ 1/n_mc 区间两者必须一致，否则检查分裂权重记账。
