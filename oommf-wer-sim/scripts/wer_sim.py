#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oommf-wer-sim 主脚本：SOT/STT-MRAM 器件写入错误率 (WER) 自动化仿真。

方法
----
1) REE（Rare Event Enhancement / importance splitting，默认）
   文献方法：Fine et al., "Write error rate of spin-transfer-torque random access
   memory including micromagnetic effects using rare event enhancement",
   IEEE Trans. Magn. (2016), arXiv:1603.08512。
   - 演化被切成若干检查点段；每段一条/多条独立 boxsi 进程（热演化器注入 300 K
     热噪声），段末 OMF 由下一段通过 Oxs_FileVectorField 重启。
   - 每个检查点按"预测坐标" m_proj = m·axiŝ（PMA 即 mz，归一化到总磁化模长）排序：
       * 已切换（m_proj < switched_mz）：判为写入成功，其权重计入"切换桶"，不再演化；
       * 滞后轨迹（m_proj > mz_th，仍停在初态附近、最可能成为稀有未切换事件）：
         分裂为 M 条子轨迹，每条权重 = 父权重 / M（权重守恒）；
       * 其余轨迹原样继续（权重不变）。
   - WER(t) = 1 - 到 t 为止的累计切换权重；脉冲结束时刻的 WER 即最终写入错误率。
   - 分裂规则可选：
       * laggard（默认，文献方法）：分裂 m_proj > mz_th 的滞后轨迹；
       * frontier（PPT 变体）：分裂最接近翻转（m_proj 最小）的若干条。
2) MC（朴素蒙特卡洛，用于验证 REE）
   - n_traj 条独立轨迹直接跑到脉冲结束，WER = 未切换条数 / 总数；
   - 只适合 WER 不太低（>= 1e-2 量级）的验证区间，REE 与 MC 的重叠区应一致。

热噪声：OOMMF 热扩展。默认 Xf_ThermHeunEvolve（纯热 Heun 积分，热统计本机
已验证正确）+ Oxs_ScriptUZeeman 电流等效场（SOT/STT 均为常数 m̂₀ 近似）；
可选 Xf_ThermSpinXferEvolve（原生 STT 力矩，T>0 热统计本机实测异常，仅零温
对照）。T>0 时扩展强制固定步长 = min_timestep（默认注入作者黄金配置
1e-14 s = 10 fs + energy_precision -1）。
随机种子：MIF 顶层 RandomSeed 命令，每段唯一种子，保证可复现。

用法
----
  python wer_sim.py --config base.yaml --out results
  python wer_sim.py --config base.yaml --sweep "J=4e11,6e11,8e11" --out results
  python wer_sim.py --config base.yaml --sweep "T=0.5e9,1e9" --out results
  python wer_sim.py --config base.yaml --dry-run
"""

import argparse
import json
import random
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml

# ---------------------------------------------------------------------------
# 环境
# ---------------------------------------------------------------------------
OOMMF_TCL = r"C:\oommf\oommf2.1\oommf\oommf\oommf.tcl"
TCLSH_CANDIDATES = [
    r"C:\ActiveTcl\bin\tclsh86t.exe",
    r"C:\ActiveTcl\bin\tclsh.exe",
    "tclsh",
]

# ---------------------------------------------------------------------------
# 默认参数（可被 --config 覆盖）
# ---------------------------------------------------------------------------
DEFAULTS = dict(
    # ---- 器件几何（与 oommf-batch-sim 一致）----
    length=40e-9, width=40e-9, thickness=1e-9, cell=5e-9,
    shape="box",                  # box | ellipse | image
    image="",
    colormap={"red": "MTJ", "white": "vacuum"},
    # ---- 磁性材料 ----
    Ms=1100e3, A=13e-12, Ku=0.85e4,
    axis=[1, 0, 0],               # 易轴 = 切换判据投影方向（PMA 用 [0,0,1]）
    alpha=0.08,
    m0=[0.996, 0.087, 0],         # 初始磁化（脉冲前的热化起点）
    use_exchange=True,            # 宏自旋验证可关（A 无效化）
    use_demag=True,               # 宏自旋验证可关（等效各向异性吸收退磁场）
    # ---- 电流 / 力矩 ----
    mode="stt",                   # stt | sot | both | none
    J=4e11,                       # STT 电流密度或 SOT 电荷电流密度 (A/m^2)
    mp=[-1, 0, 0],
    P=0.4,
    Lambda=2.0,
    t_on=1e-9, t_off=100e-9,      # 电流脉冲窗口（t_on 之前为无电流热化段）
    theta_SH=0.3, xi_FL=0.1, current_dir=[1, 0, 0],
    H_bias=[0, 0, 0],
    # ---- 热学 ----
    T=300.0,                      # 温度 K（WER 必须 > 0）
    evolver="auto",               # auto | thermheun | thermspin
                                  #   thermheun : Xf_ThermHeunEvolve（纯热，热统计可靠，推荐）
                                  #              SOT/STT 均以 ScriptUZeeman 等效场承载
                                  #   thermspin : Xf_ThermSpinXferEvolve（原生 STT 力矩，但
                                  #              T>0 热统计在本机实测异常，仅作对照）
                                  #   auto      : T>0 -> thermheun；T=0 且 mode 含 stt -> thermspin；
                                  #              T=0 其余 -> thermheun
    evolver_opts={},              # 透传 evolver 选项（用户值覆盖热默认配置）。
                                  # T>0 默认注入作者黄金配置 {min_timestep: 1e-14,
                                  #   max_timestep: 1e-14, energy_precision: -1}（固定 10fs 步长），
                                  # 如显式给 min_timestep 可提速（噪声按该值标定）。
    # ---- WER / REE ----
    method="ree",                 # ree | mc
    n_traj=200,                   # 初始轨迹数（ree: 初始批次；mc: 总轨迹数）
    checkpoint_ns=10.0,           # 检查点间隔 (ns)
    split_rule="laggard",         # laggard（文献）| frontier（PPT 变体）
    mz_th=0.95,                   # 分裂阈值（laggard: 分裂 m_proj>mz_th；frontier: 候选 m_proj<mz_th）
    split_topk=0,                 # frontier 用：每检查点最多分裂条数；0=全部候选
    split_factor=2,               # 每条分裂出的子轨迹数（文献 2；PPT 2-4）
    max_traj=500,                 # 活跃轨迹上限（超过则优先保留权重最大的子轨迹）
    switched_mz=0.0,              # m_proj 低于此值判为已切换（文献 0）
    seed=1,                       # 随机种子基准（random_seed: false 时每段派生唯一种子）
    random_seed=True,             # True=每条轨迹独立随机 seed（推荐；规避本机 Oc_Random 多流
                                  # 交替对固定 seed 的干扰，保证轨迹独立）；False=确定性派生
    # ---- 仿真控制 ----
    DeltaStage=0.2e-9,            # ODT 记录步长（每段被切成 时长/DeltaStage 个 stage）
    keep_ovf=True,                # 是否保留全部段 OMF（False 只保留每条轨迹末段）
    timeout_s=3600,
    workers=4,                    # 并行 boxsi 进程数
    threads_per_boxsi=2,          # 单进程 boxsi 线程数（避免 CPU 超订）
    basename="wer",
)

MODES = ("none", "field", "stt", "sot", "both")
SHAPES = ("box", "ellipse", "image")
METHODS = ("ree", "mc")
SPLIT_RULES = ("laggard", "frontier")


# ---------------------------------------------------------------------------
# 数值 / 配置工具
# ---------------------------------------------------------------------------
def fmt_tcl(v):
    """Tcl 安全数值格式：指数不带 '+'，如 4e11 / 1.1e6 / 4e-9。"""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, (list, tuple)):
        return " ".join(fmt_tcl(x) for x in v)
    if v == 0:
        return "0"
    s = f"{v:.12e}"
    mant, exp = s.split("e")
    mant = mant.rstrip("0").rstrip(".")
    return f"{mant}e{int(exp)}"


def fmt_val(v):
    """文件夹名安全数值格式。"""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, (list, tuple)):
        return "_".join(fmt_val(x) for x in v)
    s = f"{v:.3g}"
    if "e" in s:
        m, e = s.split("e")
        e = e.replace("+", "")
        e = str(int(e))
        s = f"{m}e{e}"
    return s


_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_INT_RE = re.compile(r"^[+-]?\d+$")


def _coerce(obj):
    """把形如 2e11 的字符串转成数值（PyYAML 会把无符号指数解析成字符串）。"""
    if isinstance(obj, str):
        s = obj.strip()
        if _INT_RE.match(s):
            try:
                return int(s)
            except ValueError:
                pass
        if _NUM_RE.match(s):
            try:
                return float(s)
            except ValueError:
                pass
        return obj
    if isinstance(obj, list):
        return [_coerce(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _coerce(v) for k, v in obj.items()}
    return obj


def load_config(path):
    cfg = dict(DEFAULTS)
    if path:
        text = Path(path).read_text(encoding="utf-8")
        if str(path).lower().endswith((".yaml", ".yml")):
            user_cfg = _coerce(yaml.safe_load(text) or {})
        else:
            user_cfg = _coerce(json.loads(text))
        if "params" in user_cfg:
            base = dict(user_cfg)
            base.update(user_cfg["params"])
            user_cfg = base
        cfg.update(user_cfg)
    _validate(cfg)
    return cfg


def _validate(cfg):
    if cfg["mode"] not in MODES:
        raise ValueError(f"mode 必须是 {MODES} 之一，收到: {cfg['mode']}")
    if cfg["shape"] not in SHAPES:
        raise ValueError(f"shape 必须是 {SHAPES} 之一，收到: {cfg['shape']}")
    if cfg["method"] not in METHODS:
        raise ValueError(f"method 必须是 {METHODS} 之一，收到: {cfg['method']}")
    if cfg["split_rule"] not in SPLIT_RULES:
        raise ValueError(f"split_rule 必须是 {SPLIT_RULES} 之一，收到: {cfg['split_rule']}")
    if cfg["T"] <= 0:
        raise ValueError("WER 仿真必须有热噪声：T 必须 > 0 K")
    if cfg["t_off"] <= cfg["t_on"]:
        raise ValueError("t_off 必须 > t_on（脉冲窗口非空）")
    if cfg["checkpoint_ns"] <= 0:
        raise ValueError("checkpoint_ns 必须 > 0")
    if cfg["n_traj"] < 1:
        raise ValueError("n_traj 必须 >= 1")


def parse_cli_value(s):
    s = s.strip()
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    low = s.lower()
    if low in ("true", "false"):
        return low == "true"
    if s.startswith("[") and s.endswith("]"):
        return [parse_cli_value(x) for x in s[1:-1].split(",")]
    return s


def parse_sweep_arg(s):
    if "=" not in s:
        raise ValueError(f"--sweep 格式应为 key=v1,v2,...，收到: {s}")
    key, vals = s.split("=", 1)
    key = key.strip()
    items = [parse_cli_value(t.strip()) for t in vals.split(",") if t.strip()]
    if not items:
        raise ValueError(f"--sweep 中 {key} 没有值")
    return key, items


def combo_name(keys, params_slice, tag=""):
    parts = [f"{k}_{fmt_val(params_slice[k])}" for k in keys]
    name = "__".join(parts) if parts else (tag or "run0")
    if tag and parts:
        name = f"{name}__{tag}"
    return name


def _clean_for_yaml(params):
    out = {}
    for k, v in params.items():
        if k in ("sweep", "runs"):
            continue
        out[k] = list(v) if isinstance(v, (list, tuple)) else v
    return out


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(x) for x in obj]
    return obj


# ---------------------------------------------------------------------------
# MIF 生成（分段、热演化、重启）
# ---------------------------------------------------------------------------
def build_segment_mif(p, traj_id, seg_index, t_start, t_end, restart_omf, seed,
                      stages):
    """生成单个轨迹单个时间段的 MIF。

    参数
    ----
    p           : 完整配置
    traj_id     : 轨迹编号（决定 basename）
    seg_index   : 段序号
    t_start/t_end : 本段覆盖的全局时间 [t_start, t_end)
    restart_omf : 非 None 时用 Oxs_FileVectorField 从该 OMF 重启（m0）
    seed        : RandomSeed 值（每段唯一）
    stages      : 本段的 stage 数（每 stage 时长 DeltaStage）
    """
    L, W, t = p["length"], p["width"], p["thickness"]
    shape = p["shape"]
    mode = p["mode"]
    bname = f"{p['basename']}_{traj_id:05d}_s{seg_index:03d}"

    # ---- Evolver 选择 ----
    evolver_cfg = str(p.get("evolver", "auto")).strip().lower()
    if evolver_cfg not in ("auto", "thermheun", "thermspin"):
        raise ValueError(f"evolver 必须是 auto|thermheun|thermspin，收到: {evolver_cfg}")
    if evolver_cfg == "thermspin":
        use_thermspin = True
    elif evolver_cfg == "thermheun":
        use_thermspin = False
    else:  # auto
        use_thermspin = (p["T"] <= 0) and mode in ("stt", "both")

    lines = ["# MIF 2.1",
             f"# Auto-generated by oommf-wer-sim (segment {seg_index}, "
             f"t=[{fmt_tcl(t_start)},{fmt_tcl(t_end)}] s)",
             ""]

    # ---- Atlas ----
    if shape == "image":
        img = p["image"].replace("\\", "/")
        cmap = " ".join(f"{c} {r}" for c, r in p["colormap"].items())
        lines += [
            "Specify Oxs_ImageAtlas:atlas {",
            f"  xrange {{0 {fmt_tcl(L)}}}",
            f"  yrange {{0 {fmt_tcl(W)}}}",
            f"  zrange {{0 {fmt_tcl(t)}}}",
            "  viewplane xy",
            f"  image {img}",
            f"  colormap {{ {cmap} }}",
            "}",
            "",
        ]
    elif shape == "ellipse":
        lines += [
            "Specify Oxs_EllipseAtlas:atlas {",
            f"  xrange {{0 {fmt_tcl(L)}}}",
            f"  yrange {{0 {fmt_tcl(W)}}}",
            f"  zrange {{0 {fmt_tcl(t)}}}",
            "  axis z",
            "}",
            "",
        ]
    else:
        lines += [
            "Specify Oxs_BoxAtlas:atlas {",
            f"  xrange {{0 {fmt_tcl(L)}}}",
            f"  yrange {{0 {fmt_tcl(W)}}}",
            f"  zrange {{0 {fmt_tcl(t)}}}",
            "}",
            "",
        ]

    # ---- Mesh ----
    lines += [
        "Specify Oxs_RectangularMesh:mesh {",
        f"  cellsize {{{fmt_tcl(p['cell'])} {fmt_tcl(p['cell'])} {fmt_tcl(t)}}}",
        "  atlas :atlas",
        "}",
        "",
    ]

    # ---- Anisotropy ----
    axis = p["axis"]
    if shape == "image":
        regions = " ".join(f"{r} {fmt_tcl(p['Ku'])}" for r in
                           set(p["colormap"].values()) - {"vacuum"})
        lines += [
            "Specify Oxs_UniaxialAnisotropy {",
            "  K1 { Oxs_AtlasScalarField {",
            "    atlas :atlas",
            "    default_value 0",
            f"    values {{ {regions} vacuum 0 }}",
            "  } }",
            f"  axis {{{fmt_tcl(axis[0])} {fmt_tcl(axis[1])} {fmt_tcl(axis[2])}}}",
            "}",
            "",
        ]
    else:
        lines += [
            "Specify Oxs_UniaxialAnisotropy {",
            f"  K1 {fmt_tcl(p['Ku'])}",
            f"  axis {{{fmt_tcl(axis[0])} {fmt_tcl(axis[1])} {fmt_tcl(axis[2])}}}",
            "}",
            "",
        ]

    # ---- Exchange / Demag（宏自旋验证可分别关闭）----
    if p.get("use_exchange", True):
        lines += [
            "Specify Oxs_UniformExchange {",
            f"  A {fmt_tcl(p['A'])}",
            "}",
            "",
        ]
    if p.get("use_demag", True):
        lines += ["Specify Oxs_Demag {}", ""]

    # ---- 恒定外磁场 ----
    if any(p["H_bias"]):
        hb = p["H_bias"]
        lines += [
            "Specify Oxs_UZeeman:bias {",
            f"  Hrange {{ {{ {fmt_tcl(hb[0])} {fmt_tcl(hb[1])} {fmt_tcl(hb[2])} "
            f"{fmt_tcl(hb[0])} {fmt_tcl(hb[1])} {fmt_tcl(hb[2])} 0 }} }}",
            "}",
            "",
        ]

    # ---- 电流等效场（ScriptUZeeman，smoothstep 包络，常数 m0 近似）----
    # thermheun 路线：SOT / STT 力矩一律以等效场承载（场不依赖瞬时 m，m0 近似）
    use_sot = mode in ("sot", "both")
    use_stt_eff = mode in ("stt", "both") and not use_thermspin
    if use_sot:
        lines += [
            "Specify Oxs_ScriptUZeeman:sot {",
            "  script_args { total_time }",
            "  script SOT_field_script",
            "}",
            "",
        ]
    if use_stt_eff:
        lines += [
            "Specify Oxs_ScriptUZeeman:stt {",
            "  script_args { total_time }",
            "  script STT_field_script",
            "}",
            "",
        ]

    # ---- Evolver ----
    # thermheun：Xf_ThermHeunEvolve（纯热 Heun 积分，T>0 热统计经本机验证可靠）
    # thermspin：Xf_ThermSpinXferEvolve（原生 STT 力矩 + 热；T>0 热统计实测异常，仅作对照）
    use_stt_native = mode in ("stt", "both") and use_thermspin
    evolver_name = "Xf_ThermSpinXferEvolve" if use_thermspin else "Xf_ThermHeunEvolve"
    lines += [f"Specify {evolver_name}:evolve {{", f"  alpha {fmt_tcl(p['alpha'])}"]
    if use_thermspin:
        lines += [
            f"  mp {{{fmt_tcl(p['mp'][0])} {fmt_tcl(p['mp'][1])} {fmt_tcl(p['mp'][2])}}}",
            f"  J {fmt_tcl(p['J']) if use_stt_native else 0.0}",
        ]
        if use_stt_native:
            lines += [
                "  J_profile Jprofile_script",
                "  J_profile_args { total_time }",
            ]
        lines += [
            f"  P {fmt_tcl(p['P'])}",
            f"  Lambda {fmt_tcl(p['Lambda'])}",
        ]
    lines += [f"  temperature {fmt_tcl(p['T'])}"]
    # T>0 热默认配置（作者 example1 黄金配置：固定 10fs 步长 + 关闭能量精度检查）；
    # 用户 evolver_opts 显式键覆盖默认。
    if p["T"] > 0:
        hot_opts = {"min_timestep": 1e-14, "max_timestep": 1e-14,
                    "energy_precision": -1,
                    "uniform_seed": int(seed)}
        hot_opts.update(p.get("evolver_opts") or {})
        for k, v in hot_opts.items():
            lines.append(f"  {k} {fmt_tcl(v)}")
    else:
        for k, v in (p.get("evolver_opts") or {}).items():
            lines.append(f"  {k} {fmt_tcl(v)}")
    lines += ["}", ""]

    # ---- 重启态（Oxs_FileVectorField 必须在 TimeDriver 之前 Specify，
    #     否则 m0 :m0f 解析时找不到该对象）----
    if restart_omf:
        lines += [
            "Specify Oxs_FileVectorField:m0f {",
            "  atlas :atlas",
            f"  file {restart_omf.replace(chr(92), '/')}",
            "}",
            "",
        ]

    # ---- Driver ----
    lines += [
        "Specify Oxs_TimeDriver {",
        f"  basename {bname}",
        "  evolver :evolve",
        f"  stopping_time {fmt_tcl(p['DeltaStage'])}",
        f"  stage_count {int(stages)}",
        "  mesh :mesh",
    ]
    if shape == "image":
        regions = " ".join(f"{r} {fmt_tcl(p['Ms'])}" for r in
                           set(p["colormap"].values()) - {"vacuum"})
        lines += [
            "  Ms { Oxs_AtlasScalarField {",
            "    default_value 0",
            "    atlas :atlas",
            f"    values {{ {regions} vacuum 0 }}",
            "  } }",
        ]
    else:
        lines += [f"  Ms {{ {fmt_tcl(p['Ms'])} }}"]
    if restart_omf:
        lines += [
            "  m0 :m0f",
            "}",
        ]
    else:
        m0 = p["m0"]
        lines += [
            f"  m0 {{{fmt_tcl(m0[0])} {fmt_tcl(m0[1])} {fmt_tcl(m0[2])}}}",
            "}",
        ]

    lines += [
        "",
        f"RandomSeed {int(seed)}",
        "",
        "Destination archive mmArchive",
        "Schedule DataTable archive stage 1",
        # stage 1 = 每个 stage 都归档磁化；OOMMF 的 Stage 事件 count 从 0 起，
        # 若写 stage N 只会在 count%N==0（含 0）输出 → 非段末。取迭代最大的 OMF 为段末。
        "Schedule Oxs_TimeDriver::Magnetization archive stage 1",
        "",
    ]

    # ---- Tcl 脚本（全局时间锚定的脉冲窗口 / 电流等效场）----
    t_on, t_off = p["t_on"], p["t_off"]
    if use_stt_native:
        lines += [
            "proc Jprofile_script { total_time } {",
            f"  set gt [expr {{{fmt_tcl(t_start)}}} + $total_time]",
            f"  if {{{fmt_tcl(t_on)} <= $gt && $gt <= {fmt_tcl(t_off)}}} {{",
            "    return 1",
            "  } else {",
            "    return 0",
            "  }",
            "}",
            "",
        ]
    if use_sot:
        hx, hy, hz = _sot_field_constants(p)
        tr = min(p["DeltaStage"] * 0.25,
                 max((p["t_off"] - p["t_on"]) * 0.25, 1e-12))
        lines += [
            "proc SOT_field_script { total_time } {",
            f"  set tstart {fmt_tcl(t_start)}",
            f"  set ton {fmt_tcl(t_on)}",
            f"  set toff {fmt_tcl(t_off)}",
            f"  set tr {fmt_tcl(tr)}",
            f"  set hx {fmt_tcl(hx)}",
            f"  set hy {fmt_tcl(hy)}",
            f"  set hz {fmt_tcl(hz)}",
            "  set t [expr {$tstart + $total_time}]",
            "  if {$t <= $ton} {",
            "    set f 0.0",
            "    set df 0.0",
            "  } elseif {$t < $ton + $tr} {",
            "    set x [expr {($t - $ton)/$tr}]",
            "    set f [expr {$x*$x*(3-2*$x)}]",
            "    set df [expr {6*$x*(1-$x)/$tr}]",
            "  } elseif {$t <= $toff - $tr} {",
            "    set f 1.0",
            "    set df 0.0",
            "  } elseif {$t < $toff} {",
            "    set x [expr {($toff - $t)/$tr}]",
            "    set f [expr {$x*$x*(3-2*$x)}]",
            "    set df [expr {-6*$x*(1-$x)/$tr}]",
            "  } else {",
            "    set f 0.0",
            "    set df 0.0",
            "  }",
            "  return [list [expr {$hx*$f}] [expr {$hy*$f}] [expr {$hz*$f}] \\",
            "               [expr {$hx*$df}] [expr {$hy*$df}] [expr {$hz*$df}]]",
            "}",
            "",
        ]
    if use_stt_eff:
        hx, hy, hz = _stt_field_constants(p)
        tr = min(p["DeltaStage"] * 0.25,
                 max((p["t_off"] - p["t_on"]) * 0.25, 1e-12))
        lines += [
            "proc STT_field_script { total_time } {",
            f"  set tstart {fmt_tcl(t_start)}",
            f"  set ton {fmt_tcl(t_on)}",
            f"  set toff {fmt_tcl(t_off)}",
            f"  set tr {fmt_tcl(tr)}",
            f"  set hx {fmt_tcl(hx)}",
            f"  set hy {fmt_tcl(hy)}",
            f"  set hz {fmt_tcl(hz)}",
            "  set t [expr {$tstart + $total_time}]",
            "  if {$t <= $ton} {",
            "    set f 0.0",
            "    set df 0.0",
            "  } elseif {$t < $ton + $tr} {",
            "    set x [expr {($t - $ton)/$tr}]",
            "    set f [expr {$x*$x*(3-2*$x)}]",
            "    set df [expr {6*$x*(1-$x)/$tr}]",
            "  } elseif {$t <= $toff - $tr} {",
            "    set f 1.0",
            "    set df 0.0",
            "  } elseif {$t < $toff} {",
            "    set x [expr {($toff - $t)/$tr}]",
            "    set f [expr {$x*$x*(3-2*$x)}]",
            "    set df [expr {-6*$x*(1-$x)/$tr}]",
            "  } else {",
            "    set f 0.0",
            "    set df 0.0",
            "  }",
            "  return [list [expr {$hx*$f}] [expr {$hy*$f}] [expr {$hz*$f}] \\",
            "               [expr {$hx*$df}] [expr {$hy*$df}] [expr {$hz*$df}]]",
            "}",
            "",
        ]

    return "\n".join(lines)


def _sot_field_constants(p):
    """SOT 等效场（常数 m0 近似，与 oommf-batch-sim 一致）：
    Hmag = hbar*theta_SH*J / (2 e mu0 Ms t)；H_DL = Hmag*(sigma×m0̂)；H_FL = xi_FL*Hmag*sigma。
    定量 SOT-WER 需自定义 SOT 扩展（随瞬时 m 更新），本技能用于趋势研究。
    """
    import numpy as np
    hbar = 1.054571817e-34
    e = 1.602176634e-19
    mu0 = 4 * np.pi * 1e-7
    sigma = np.array(p["current_dir"], dtype=float)
    sigma /= np.linalg.norm(sigma)
    m0 = np.array(p["m0"], dtype=float)
    if np.linalg.norm(m0) > 0:
        m0 = m0 / np.linalg.norm(m0)
    Hmag = hbar * p["theta_SH"] * p["J"] / (2 * e * mu0 * p["Ms"] * p["thickness"])
    H = Hmag * np.cross(sigma, m0) + p["xi_FL"] * Hmag * sigma
    return float(H[0]), float(H[1]), float(H[2])


def _stt_field_constants(p):
    """STT 等效场（常数近似，供 thermheun 路线承载 STT）：
    Slonczewski 力矩 T = (hbar J P/(2 e Ms t)) m×(mp×m)，取 g≈P 常数、m≈m0̂，
    等价场 H_STT = -(hbar P J/(2 e mu0 Ms t))·(mp̂×m0̂)。
    局限：真实 STT 随瞬时 m 变化（m×(mp×m)），定量 STT-WER 请用 evolver=thermspin
    或自定义扩展。
    """
    import numpy as np
    hbar = 1.054571817e-34
    e = 1.602176634e-19
    mu0 = 4 * np.pi * 1e-7
    mp = np.array(p["mp"], dtype=float)
    n = np.linalg.norm(mp)
    if n > 0:
        mp = mp / n
    m0 = np.array(p["m0"], dtype=float)
    if np.linalg.norm(m0) > 0:
        m0 = m0 / np.linalg.norm(m0)
    Hmag = hbar * p["P"] * p["J"] / (2 * e * mu0 * p["Ms"] * p["thickness"])
    H = -Hmag * np.cross(mp, m0)
    return float(H[0]), float(H[1]), float(H[2])


# ---------------------------------------------------------------------------
# boxsi 执行（Popen + 超时 kill，避免进程挂起）
# ---------------------------------------------------------------------------
def find_tclsh():
    import shutil
    for cand in TCLSH_CANDIDATES:
        p = Path(cand)
        if p.exists() or shutil.which(cand):
            return cand
    raise RuntimeError("未找到 tclsh，请检查 ActiveTcl 安装")


def run_boxsi(run_dir, mif_name, timeout_s, threads=2):
    """在 run_dir 中运行 boxsi；返回 (ok, detail)。

    -tk 0：关闭 Tk GUI（无窗口、无错误对话框等待）；
    -exitondone 1：问题求解完成即退出进程（否则窗口常驻导致超时误判）；
    -threads 限制单进程线程数，避免 workers×默认线程数超订 CPU；
    -logfile 独立错误日志：多个 boxsi 并行时避免竞争全局 boxsi.errors。
    """
    tclsh = find_tclsh()
    logf = open(run_dir / "boxsi.log", "w", encoding="utf-8", errors="replace")
    errf = open(run_dir / "boxsi.err", "w", encoding="utf-8", errors="replace")
    try:
        cmd = [tclsh, OOMMF_TCL, "boxsi", "-tk", "0", "-exitondone", "1"]
        cmd += ["-logfile", str(run_dir / "boxsi.errors")]
        if threads and threads >= 1:
            cmd += ["-threads", str(threads)]
        cmd.append(mif_name)
        proc = subprocess.Popen(
            cmd,
            cwd=str(run_dir), stdout=logf, stderr=errf)
        try:
            proc.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            proc.kill()
            return False, f"timeout 超过 {timeout_s}s"
        ok = proc.returncode == 0
        detail = ""
        if not ok:
            errf.flush()
            with open(run_dir / "boxsi.err", encoding="utf-8",
                      errors="replace") as f:
                detail = f.read()[-800:]
        return ok, detail
    finally:
        logf.close()
        errf.close()


def find_latest_omf(run_dir, bname):
    """取本段最后写的磁化归档文件（命名 <basename>-Oxs_TimeDriver-Magnetization-<stage>-<iter>.omf）。"""
    hits = sorted(run_dir.glob(f"{bname}-Oxs_TimeDriver-Magnetization-*.omf"))
    if hits:
        return str(hits[-1])
    hits = sorted(run_dir.glob(f"{bname}-*.ovf"))
    return str(hits[-1]) if hits else None


# ---------------------------------------------------------------------------
# ODT 解析
# ---------------------------------------------------------------------------
def parse_odt(path):
    """返回 (col_names, data_matrix)。"""
    cols = []
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if line.startswith("# Columns:"):
                    cols = _parse_columns(line[len("# Columns:"):])
                continue
            if not line.strip():
                continue
            vals = line.split()
            if len(vals) == len(cols):
                rows.append([float(v) for v in vals])
    return cols, rows


def _parse_columns(header):
    tokens = re.findall(r"\{[^}]*\}|\S+", header)
    return [t.strip("{} ").strip() for t in tokens]


def find_col(cols, *fragments):
    for i, c in enumerate(cols):
        if all(frag in c for frag in fragments):
            return i
    return None


def odt_last_state(cols, rows):
    """取 ODT 末行的 (t, m_vec, proj)。列名兼容
    Xf_ThermSpinXferEvolve:evolve:... / Oxs_TimeDriver::mx 等。
    """
    t_i = find_col(cols, "Simulation time")
    mx_i, my_i, mz_i = find_col(cols, "::mx"), find_col(cols, "::my"), find_col(cols, "::mz")
    if t_i is None or mx_i is None or not rows:
        raise ValueError("ODT 缺少必需列（Simulation time / ::mx）")
    r = rows[-1]
    m = (r[mx_i], r[my_i], r[mz_i])
    norm = (m[0] ** 2 + m[1] ** 2 + m[2] ** 2) ** 0.5
    return r[t_i], m, norm


# ---------------------------------------------------------------------------
# REE / MC 核心
# ---------------------------------------------------------------------------
class Trajectory:
    """一条（可能被多次分裂的）轨迹。

    weight>0 恒成立；切换/未切换都由权重累计进 WER。
    """

    def __init__(self, tid, parent, weight, gen):
        self.id = tid
        self.parent = parent          # 父轨迹 id 或 None
        self.weight = weight
        self.gen = gen                # 分裂代数（0=初始）
        self.cur_omf = None           # 当前段末 OMF（下一段重启用）
        self.mz_trace = []            # [(t, m_proj), ...]
        self.switched_at = None       # 切换的全局时刻
        self.status = "active"        # active | switched
        self.seg_count = 0

    def to_json(self):
        return {"id": self.id, "parent": self.parent, "weight": self.weight,
                "gen": self.gen, "switched_at": self.switched_at,
                "status": self.status, "seg_count": self.seg_count,
                "mz_trace": [[round(t, 12), round(v, 6)] for t, v in self.mz_trace]}


def make_time_grid(t_start, t_end, checkpoint_s):
    """检查点网格：[t_start, t1, ..., t_end]。"""
    n = max(1, int(round((t_end - t_start) / checkpoint_s)))
    pts = [t_start + checkpoint_s * i for i in range(n)]
    pts.append(t_end)
    return sorted(set(pts))


def run_wer_combo(cfg, outdir, dry_run=False):
    """跑一个参数组合的 WER 仿真；outdir 为该组合的目录。"""
    t0 = time.time()
    method = cfg["method"]
    n_traj = int(cfg["n_traj"])
    checkpoint_s = cfg["checkpoint_ns"] * 1e-9
    t_end = cfg["t_off"]           # WER 在脉冲结束时刻报告
    grid = make_time_grid(0.0, t_end, checkpoint_s)

    # 轨迹初始化
    total_w = 1.0
    trajs = {i: Trajectory(i, None, total_w / n_traj, 0) for i in range(n_traj)}
    next_tid = n_traj
    switched_weight = 0.0
    wer_trace = []                # [(t, wer, n_active, n_switched_traj)]
    all_traj_records = []
    failed_ids = []               # 段失败轨迹（权重计入 WER 的保守估计）

    print(f"  [{method.upper()}] 时间网格: {len(grid)-1} 段, "
          f"每段 {fmt_val(checkpoint_s*1e9)} ns, t_end={fmt_val(t_end*1e9)} ns, "
          f"初始 {n_traj} 条轨迹")

    for k in range(len(grid) - 1):
        ts, te = grid[k], grid[k + 1]
        seg_dur = te - ts
        stages = max(1, int(round(seg_dur / max(cfg["DeltaStage"], 1e-15))))
        tasks = []
        for tid, tr in trajs.items():
            if tr.status != "active":
                continue
            tasks.append((tid, tr))
        print(f"    段 {k}: t=[{fmt_val(ts*1e9)},{fmt_val(te*1e9)}] ns, "
              f"活跃 {len(tasks)} 条")

        # ---- 并行跑本段 ----
        results = _run_segments(cfg, tasks, ts, te, stages, outdir,
                                dry_run=dry_run)
        if dry_run:
            continue

        # ---- 读取末态 / 判定 / 分裂 ----
        new_trajs = {}
        for tid, tr in tasks:
            r = results[tid]
            if r is None:
                # 段失败：不再演化，权重按"未切换"计入 WER（保守、权重守恒）
                tr.status = "failed"
                failed_ids.append(tid)
                all_traj_records.append(tr.to_json())
                continue
            t_final = r["t_final"]
            m = r["m"]
            norm = r["norm"]
            proj = _project(m, cfg["axis"], norm)
            # mz_trace 使用全局段末时间 te（ODT 内的时间为段内 t，多段下需对齐）
            tr.mz_trace.append((te, proj))
            tr.seg_count += 1
            tr.cur_omf = r["omf"]
            if proj < cfg["switched_mz"]:
                tr.status = "switched"
                tr.switched_at = te
                switched_weight += tr.weight
                all_traj_records.append(tr.to_json())
            else:
                new_trajs[tid] = tr

        if method == "ree" and not dry_run and k < len(grid) - 2:
            # 仅在有下一段时才分裂（最后一段的分裂不会继续演化，只会污染记账）
            new_trajs = _apply_splitting(cfg, new_trajs, next_tid,
                                         switched_weight, all_traj_records)
            next_tid = _next_free_id(new_trajs, next_tid)

        trajs = new_trajs
        n_active = sum(1 for t in trajs.values() if t.status == "active")
        wer = max(0.0, total_w - switched_weight)
        wer_trace.append((te, wer, n_active,
                          int(round(switched_weight / (total_w / n_traj))))

        # 归档每段结果（ODT 复制到组合目录的 trajectories/ 下）
        if not dry_run:
            _archive_segments(cfg, tasks, results, outdir, k)

    # keep_ovf=false：只保留每轨迹最后一段 OMF
    if not dry_run:
        _cleanup_omf(outdir, cfg["basename"], cfg.get("keep_ovf", True))

    # ---- 汇总 ----
    wer_final = max(0.0, total_w - switched_weight)
    min_w = min((t.weight for t in trajs.values()
                 if t.status == "active"), default=0.0)
    summary = {
        "method": method,
        "params": _jsonable(_clean_for_yaml(cfg)),
        "wer_final": wer_final,
        "wer_at_t_end": wer_final,
        "switched_weight": switched_weight,
        "n_initial": n_traj,
        "n_final_active": sum(1 for t in trajs.values()
                              if t.status == "active"),
        "n_failed": len(failed_ids),
        "failed_ids": failed_ids,
        "n_switched_traj": int(round(switched_weight / (total_w / n_traj))),
        "min_active_weight": min_w,
        "effective_weight_resolution": min_w,
        "wer_trace": [[round(t, 12), round(w, 12), n, s]
                      for t, w, n, s in wer_trace],
        "elapsed_s": round(time.time() - t0, 1),
        "checkpoint_ns": cfg["checkpoint_ns"],
        "mz_th": cfg["mz_th"],
        "split_rule": cfg["split_rule"],
        "split_factor": cfg["split_factor"],
        "T_K": cfg["T"],
    }
    if not dry_run:
        with open(outdir / "summary.json", "w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        (outdir / "wer_trace.csv").write_text(
            "t_s,wer,n_active,n_switched\n" +
            "\n".join(f"{t},{w},{n},{s}" for t, w, n, s in wer_trace),
            encoding="utf-8")
        for tr in trajs.values():
            all_traj_records.append(tr.to_json())
        (outdir / "trajectories.json").write_text(
            json.dumps(all_traj_records, indent=2, ensure_ascii=False),
            encoding="utf-8")
    return summary


def _project(m, axis, norm):
    ax = axis
    ln = (ax[0] ** 2 + ax[1] ** 2 + ax[2] ** 2) ** 0.5 or 1.0
    mnorm = norm if norm > 0 else 1.0
    return (m[0] * ax[0] + m[1] * ax[1] + m[2] * ax[2]) / (ln * mnorm)


def _run_segments(cfg, tasks, ts, te, stages, outdir, dry_run=False):
    """并行执行一批段的 boxsi；返回 {tid: dict(omf, t_final, m, norm)}。"""
    results = {}
    workers = max(1, int(cfg["workers"]))

    def one(task):
        tid, tr = task
        run_dir = outdir / "trajectories" / f"traj_{tid:05d}"
        run_dir.mkdir(parents=True, exist_ok=True)
        # seed：默认随机（本机 OOMMF 2.1a1 的 Oc_Random 存在"多流交替"现象，
        # 同 MIF 可能未尊重固定 seed；随机化保证各轨迹流独立、WER 统计无偏）。
        # 设 random_seed: false 可退回确定性派生（结果仍可能受 boxsi 环境影响）。
        if cfg.get("random_seed", True):
            seed = random.randrange(1, 2 ** 31 - 1)
        else:
            seed = (cfg["seed"] + tid * 7919 + tr.seg_count * 104729) % (2 ** 31 - 1)
        bname = f"{cfg['basename']}_{tid:05d}_s{tr.seg_count:03d}"
        restart = tr.cur_omf  # 上一段末 OMF（None=首段用 m0 常量）
        mif = build_segment_mif(cfg, tid, tr.seg_count, ts, te, restart, seed,
                                stages)
        (run_dir / f"{bname}.mif").write_text(mif, encoding="utf-8")
        if dry_run:
            return tid, {"dry": True}
        ok, detail = run_boxsi(run_dir, f"{bname}.mif", cfg["timeout_s"],
                               int(cfg.get("threads_per_boxsi", 2)))
        if not ok:
            (run_dir / "FAILED.txt").write_text(detail or "boxsi 失败",
                                                encoding="utf-8")
            return tid, None
        omf = find_latest_omf(run_dir, bname)
        odt = run_dir / f"{bname}.odt"
        cols, rows = parse_odt(odt)
        t_final, m, norm = odt_last_state(cols, rows)
        return tid, {"omf": omf, "t_final": t_final, "m": m, "norm": norm,
                     "odt": str(odt)}

    if workers == 1:
        for task in tasks:
            tid, r = one(task)
            results[tid] = r
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(one, t) for t in tasks]
            for fut in as_completed(futs):
                tid, r = fut.result()
                results[tid] = r
    return results


def _apply_splitting(cfg, trajs, next_tid, switched_weight, records=None):
    """重要性分裂：按 split_rule 对活跃轨迹分裂，权重守恒。

    laggard（文献）：m_proj > mz_th 的滞后轨迹最可能成为稀有未切换事件 → 分裂。
    frontier（PPT 变体）：m_proj < mz_th 的接近翻转轨迹 → 分裂（按最接近优先）。
    受 max_traj 限制：候选按优先级排序后，只分裂到不超过上限为止。
    records：若提供（all_traj_records），被分裂的父轨迹会以 status="split" 记录，
            保证 trajectories.json 中父子关系可追溯。
    """
    split_factor = max(1, int(cfg["split_factor"]))
    max_traj = int(cfg["max_traj"])
    rule = cfg["split_rule"]
    mz_th = cfg["mz_th"]

    active = [tr for tr in trajs.values() if tr.status == "active"]
    if len(active) < 2 or split_factor <= 1:
        return trajs

    if rule == "laggard":
        cands = [tr for tr in active if tr.mz_trace[-1][1] > mz_th]
        cands.sort(key=lambda tr: tr.mz_trace[-1][1], reverse=True)  # 最滞后优先
    else:  # frontier
        cands = [tr for tr in active if tr.mz_trace[-1][1] < mz_th]
        cands.sort(key=lambda tr: tr.mz_trace[-1][1])               # 最接近翻转优先
        topk = int(cfg.get("split_topk", 0))
        if topk > 0:
            cands = cands[:topk]

    if not cands:
        return trajs

    # 计算最多能分裂多少条而不超过 max_traj
    # 分裂后活跃数 = (n_active - x) + x*split_factor = n_active + x*(split_factor-1) <= max_traj
    n_active = len(active)
    budget = max_traj - n_active
    max_split = max(0, budget // (split_factor - 1)) if split_factor > 1 else 0
    cands = cands[:max_split]
    if not cands:
        return trajs

    new_trajs = dict(trajs)
    for tr in cands:
        # 文献语义：父轨迹被 M 条等权重子轨迹替换（权重守恒 w -> M×(w/M)）
        del new_trajs[tr.id]
        tr.status = "split"
        if records is not None:
            records.append(tr.to_json())
        child_w = tr.weight / split_factor
        for _ in range(split_factor):
            child = Trajectory(next_tid, tr.id, child_w, tr.gen + 1)
            child.cur_omf = tr.cur_omf     # 子轨迹从父轨迹当前状态继续
            child.mz_trace = list(tr.mz_trace)
            new_trajs[next_tid] = child
            next_tid += 1
    return new_trajs


def _next_free_id(trajs, next_tid):
    while next_tid in trajs:
        next_tid += 1
    return next_tid


def _archive_segments(cfg, tasks, results, outdir, seg_k):
    """把本段 ODT 复制到组合目录（轻量归档），并清理中间 OMF。"""
    traj_dir = outdir / "trajectories"
    for tid, r in results.items():
        if not r or not r.get("odt"):
            continue
        src = Path(r["odt"])
        dst = traj_dir / f"traj_{tid:05d}" / f"seg{seg_k:03d}.odt"
        try:
            if not dst.exists():
                dst.write_text(src.read_text(encoding="utf-8", errors="replace"),
                               encoding="utf-8")
        except OSError:
            pass


def _cleanup_omf(outdir, basename, keep_ovf):
    """keep_ovf=false 时，每个轨迹只保留最后一段 OMF（重启点），其余删除。"""
    if keep_ovf:
        return
    traj_dir = outdir / "trajectories"
    if not traj_dir.exists():
        return
    for d in sorted(traj_dir.glob("traj_*")):
        omfs = sorted(d.glob(f"{basename}_*.omf"))
        for omf in omfs[:-1]:
            try:
                omf.unlink()
            except OSError:
                pass


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="OOMMF SOT/STT-MRAM WER 仿真 (REE/MC)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--out", default="wer_results")
    ap.add_argument("--sweep", action="append", default=[])
    ap.add_argument("--tag", default="")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    if args.sweep:
        cli_sweep = dict(cfg.get("sweep") or {})
        for s in args.sweep:
            key, vals = parse_sweep_arg(s)
            cli_sweep[key] = vals
        cfg["sweep"] = cli_sweep
    if args.tag:
        cfg["tag"] = args.tag
    if args.workers:
        cfg["workers"] = args.workers

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "batch_config.yaml").write_text(
        yaml.safe_dump(_clean_for_yaml(cfg), allow_unicode=True, sort_keys=False),
        encoding="utf-8")

    # 组合展开（sweep 笛卡尔积）
    sweep = cfg.get("sweep") or {}
    if sweep:
        import itertools
        keys = list(sweep.keys())
        value_lists = [sweep[k] if isinstance(sweep[k], list) else [sweep[k]]
                       for k in keys]
        combos = []
        for values in itertools.product(*value_lists):
            over = dict(zip(keys, values))
            params = dict(cfg)
            params.update(over)
            combos.append((combo_name(keys, over, cfg.get("tag", "")), params))
    else:
        combos = [(cfg.get("tag") or "run0", dict(cfg))]

    print("=" * 64)
    print(f"OOMMF WER 仿真  共 {len(combos)} 个组合")
    print(f"方法={cfg['method']}  T={cfg['T']} K  器件={cfg['length']*1e9:g}x"
          f"{cfg['width']*1e9:g}x{cfg['thickness']*1e9:g} nm  "
          f"模式={cfg['mode']}  J={cfg['J']:.2e} A/m²")
    if sweep:
        print("扫描参数: " + ", ".join(f"{k}={sweep[k]}" for k in sweep))
    print("-" * 64)

    results = []
    for name, params in combos:
        print(f"组合: {name}")
        combo_dir = outdir / name
        combo_dir.mkdir(parents=True, exist_ok=True)
        s = run_wer_combo(params, combo_dir, dry_run=args.dry_run)
        s["combo"] = name
        results.append(s)
        if not args.dry_run:
            print(f"  -> WER(t_end) = {s['wer_final']:.3e}  "
                  f"({s['elapsed_s']:.1f}s)")
            if s.get("n_failed"):
                print(f"  [warn] {s['n_failed']} 条轨迹段失败 "
                      f"(ids={s['failed_ids']})：其权重按未切换保守计入 WER，"
                      f"请检查 trajectories/traj_*/FAILED.txt")

    if args.dry_run:
        print("\ndry-run：不执行仿真。")
        return

    # 汇总
    with open(outdir / "summary_all.json", "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    # 数据处理（出图 / 报告）
    sys.path.insert(0, str(Path(__file__).parent))
    try:
        import process_wer
        process_wer.process(outdir)
    except Exception as exc:
        print(f"[warn] 数据处理失败: {exc}")


if __name__ == "__main__":
    main()
