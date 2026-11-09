# OOMMF MRAM 仿真技能包

OOMMF 微磁学仿真的三个可复用技能（Skill），面向 SOT-MRAM / STT-MRAM 器件的写入电流阈值、批量参数扫描与写入错误率（WER）评估。

每个技能是一个独立目录：`SKILL.md`（技能说明，agent 读取后触发）、`scripts/`（可独立运行的 Python 脚本）、`references/`（方法与物理参考）。

## 目录结构

```
oommf-batch-sim/                批量零温参数扫描
├── SKILL.md
├── references/batch_workflow.md
└── scripts/
    ├── batch_sim.py            # 生成 MIF + 批量跑 boxsi + 自动归档
    └── process_results.py      # ODT 解析 / 指标 / 汇总 / 出图 / Jc 粗提取

oommf-sot-jc/                   单器件 Jc 二分精修
├── SKILL.md
├── references/sot_physics.md
└── scripts/extract_jc.py

oommf-wer-sim/                  写入错误率 WER（REE 重要性分裂 + 朴素 MC 交叉验证）
├── SKILL.md
├── references/ree_method.md
└── scripts/
    ├── wer_sim.py
    └── process_wer.py
```

## 分工

| 需求 | 用哪个技能 |
|---|---|
| WER / WER(J) / 写入可靠性 | `oommf-wer-sim` |
| 任意参数批量零温扫描、Jc 粗估 | `oommf-batch-sim` |
| 单器件 Jc 二分精修 | `oommf-sot-jc` |

## 环境依赖

- OOMMF 2.1a1：`C:/oommf/oommf2.1/oommf/oommf/oommf.tcl`
  - WER 仿真需已编译热扩展 `Xf_ThermHeunEvolve`（默认）与 `Xf_ThermSpinXferEvolve`（可选）
- tclsh：`C:/ActiveTcl/bin/tclsh86t.exe`
- Python 3 + `numpy` / `matplotlib` / `PyYAML`
  - `oommf-sot-jc/scripts/extract_jc.py` 可选 Ubermag（`oommfc`），未安装时自动回退直接写 MIF + boxsi
- 脚本顶部的 `OOMMF_TCL` / `TCLSH_CANDIDATES` 按本机路径写死，换机器需相应修改

## 安装

把技能目录整体复制到你的 skills 根目录下，例如把 `oommf-wer-sim/` 放到 `<skills_root>/oommf-wer-sim/`，保持目录内部结构不变（`SKILL.md` 相对引用 `scripts/` 与 `references/`）。

## 已验证的关键约定

详见各 `SKILL.md` 的注意事项一节，几条最容易踩的：

1. **热引擎选型**：宏自旋 / 单格场景用 `Xf_ThermHeunEvolve`（SOT/STT 力矩以 `Oxs_ScriptUZeeman` 等效场承载，热统计已验证）；真实多格原生 STT 器件（带逐区 mp/J、DMI）可用 `Xf_ThermSpinXferEvolve`。
2. **boxsi 参数**：必须 `-tk 0 -exitondone 1`，否则 Tk 窗口常驻、进程不退出，被误判为超时。
3. **WER 时长**：必须覆盖完整翻转时间（含弛豫），否则 WER≈1 是截断伪影而非物理结论。
4. **随机种子**：OOMMF 2.1a1 的 `Oc_Random` 存在多流交替，同 seed 不可复现；技能默认 `random_seed: true` 保证轨迹独立、WER 统计无偏。
5. **STT Lambda**：OOMMF 默认为 2.0（不是 1），设 1 会弱化力矩。
6. **网格**：cell 必须 ≤ 交换长度一半，否则空间模式丢失、切换阈值失真。

## 说明

个人科研用脚本，参数默认值取自 HM/CoFeB/MgO 体系与 SOT-MRAM 典型量级，使用前请按自己的器件参数核对。
