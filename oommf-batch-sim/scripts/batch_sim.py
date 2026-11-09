#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oommf-batch-sim 主脚本：批量 OOMMF 微磁仿真。

功能：
  1. 从基础配置生成 OOMMF MIF 项目（支持 box / ellipse / image 三种器件形状，
     none / field / stt / sot / both 五种力矩模式）
  2. 参数扫描（sweep 笛卡尔积）或显式运行列表（runs）批量修改任意参数
  3. 每个参数组合独立命名文件夹归档（run.mif / params.yaml / OOMMF 输出）
  4. 调用 process_results.py 自动完成数据处理（ODT 解析 / 汇总 / 出图 / Jc 提取）

用法：
  python batch_sim.py --config base.yaml --out results
  python batch_sim.py --config base.yaml --sweep "J=4e11,6e11,8e11" --out results
  python batch_sim.py --config base.yaml --sweep "J=4e11,6e11" --sweep "Ms=1e6,1.1e6" \
                      --tag runA --out results
  python batch_sim.py --config base.yaml --runs runs.yaml --out results
  python batch_sim.py --config base.yaml --dry-run          # 只打印运行计划
"""

import argparse
import json
import re
import shutil
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
    # ---- 器件几何 ----
    length=40e-9, width=40e-9, thickness=1e-9, cell=5e-9,
    shape="box",                  # box | ellipse | image
    image="",                     # shape=image 时的位图绝对路径（png/jpg）
    colormap={"red": "MTJ", "white": "vacuum"},   # 位图颜色 -> 区域名
    # ---- 磁性材料 ----
    Ms=1100e3, A=13e-12, Ku=0.85e4,
    axis=[1, 0, 0],               # 易轴（也决定切换判据投影方向）
    alpha=0.08,
    m0=[0.996, 0.087, 0],         # 初始磁化方向
    # ---- 电流 / 力矩 ----
    mode="stt",                   # none | field | stt | sot | both
    J=4e11,                       # STT 电流密度 或 SOT 电荷电流密度 (A/m²)
    mp=[-1, 0, 0],                # STT 自旋极化方向
    P=0.4,                        # STT 极化率（OOMMF 默认 0.4）
    Lambda=2.0,                   # STT Lambda（OOMMF 默认 2；1 会弱化力矩）
    t_on=1e-9, t_off=4e-9,        # 电流脉冲窗口 (s)
    theta_SH=0.3,                 # SOT 自旋霍尔角
    xi_FL=0.1,                    # SOT 类场 / 类阻尼比
    current_dir=[1, 0, 0],        # SOT 电荷电流方向
    H_bias=[0, 0, 0],             # 恒定外磁场 (A/m)
    # ---- 仿真控制 ----
    stage_count=50,               # 阶段数
    DeltaStage=0.2e-9,            # 每阶段时长 (s)
    T=0.0,                        # 温度（0=零温 RungeKutta；>0 需要热扩展，暂不支持）
    basename="run",               # OOMMF 输出文件名前缀
    schedule_ovf=False,           # 是否输出磁化归档 .ovf
    timeout_s=3600,               # 单个仿真的超时秒数
    workers=1,                    # 并行仿真进程数
)

MODES = ("none", "field", "stt", "sot", "both")
SHAPES = ("box", "ellipse", "image")


# ---------------------------------------------------------------------------
# 数值格式化
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
    """文件夹名安全数值格式：1.1e6 / 4e-9 / -1_0_0（列表下划线连接）。"""
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


# ---------------------------------------------------------------------------
# 配置加载
# ---------------------------------------------------------------------------
_NUM_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_INT_RE = re.compile(r"^[+-]?\d+$")


def _coerce(obj):
    """递归把形如 2e11 / 4e-9 的字符串转成数值。

    PyYAML 只把带符号指数的写法（如 2e+11）识别为 float，
    不带符号的 2e11 会解析成字符串，这里统一纠正。
    """
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
    if isinstance(obj, tuple):
        return tuple(_coerce(x) for x in obj)
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
        if "params" in user_cfg:            # 允许整体包一层 params
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
    if cfg["shape"] == "image" and not cfg["image"]:
        raise ValueError("shape=image 时必须提供 image 位图路径")
    if cfg["stage_count"] <= 0:
        raise ValueError("stage_count 必须 > 0（或设置 total_time 自动换算）")


def parse_cli_value(s):
    """解析 --sweep 中的单个值：数字 / 列表 / 布尔。"""
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
    """解析 --sweep "J=4e11,6e11,8e11" 或 "m0=[0,0,1],[1,0,0]"。"""
    if "=" not in s:
        raise ValueError(f"--sweep 格式应为 key=v1,v2,...，收到: {s}")
    key, vals = s.split("=", 1)
    key = key.strip()
    items = []
    for tok in vals.split(","):
        tok = tok.strip()
        if not tok:
            continue
        items.append(parse_cli_value(tok))
    if not items:
        raise ValueError(f"--sweep 中 {key} 没有值")
    return key, items


# ---------------------------------------------------------------------------
# MIF 生成
# ---------------------------------------------------------------------------
def build_mif(p):
    """根据参数生成完整 OOMMF MIF 2.1 文本。"""
    L, W, t = p["length"], p["width"], p["thickness"]
    shape = p["shape"]
    mode = p["mode"]
    lines = ["# MIF 2.1",
             f"# Auto-generated by oommf-batch-sim",
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
        # Oxs_EllipseAtlas: margin 默认 0 -> 椭圆内切于盒子（半径 L/2, W/2）
        lines += [
            "Specify Oxs_EllipseAtlas:atlas {",
            f"  xrange {{0 {fmt_tcl(L)}}}",
            f"  yrange {{0 {fmt_tcl(W)}}}",
            f"  zrange {{0 {fmt_tcl(t)}}}",
            "  axis z",                       # 圆柱对称轴（薄膜默认 z）
            "}",
            "",
        ]
    else:  # box
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

    # ---- Exchange / Demag ----
    lines += [
        "Specify Oxs_UniformExchange {",
        f"  A {fmt_tcl(p['A'])}",
        "}",
        "",
        "Specify Oxs_Demag {}",
        "",
    ]

    # ---- 恒定外磁场（UZeeman 恒场也需 Hrange，steps=0 表示恒定）----
    if any(p["H_bias"]):
        hb = p["H_bias"]
        lines += [
            "Specify Oxs_UZeeman:bias {",
            f"  Hrange {{ {{ {fmt_tcl(hb[0])} {fmt_tcl(hb[1])} {fmt_tcl(hb[2])} "
            f"{fmt_tcl(hb[0])} {fmt_tcl(hb[1])} {fmt_tcl(hb[2])} 0 }} }}",
            "}",
            "",
        ]

    # ---- SOT 等效场（时变 ScriptUZeeman，返回场+时间导数）----
    use_sot = mode in ("sot", "both")
    if use_sot:
        lines += [
            "Specify Oxs_ScriptUZeeman:sot {",
            "  script_args { total_time }",
            "  script SOT_field_script",
            "}",
            "",
        ]

    # ---- Evolver ----
    use_stt = mode in ("stt", "both")
    if use_stt:
        mp = p["mp"]
        lines += [
            "Specify Oxs_SpinXferEvolve:evolve {",
            f"  alpha {fmt_tcl(p['alpha'])}",
            f"  mp {{{fmt_tcl(mp[0])} {fmt_tcl(mp[1])} {fmt_tcl(mp[2])}}}",
            f"  J {fmt_tcl(p['J'])}",
            "  J_profile Jprofile_script",
            "  J_profile_args { total_time }",
            f"  P {fmt_tcl(p['P'])}",
            f"  Lambda {fmt_tcl(p['Lambda'])}",
            "}",
            "",
        ]
    else:
        lines += [
            "Specify Oxs_RungeKuttaEvolve:evolve {",
            f"  alpha {fmt_tcl(p['alpha'])}",
            "}",
            "",
        ]

    # ---- Driver ----
    m0 = p["m0"]
    lines += [
        "Specify Oxs_TimeDriver {",
        f"  basename {p['basename']}",
        "  evolver :evolve",
        f"  stopping_time {fmt_tcl(p['DeltaStage'])}",
        f"  stage_count {int(p['stage_count'])}",
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
    lines += [
        f"  m0 {{{fmt_tcl(m0[0])} {fmt_tcl(m0[1])} {fmt_tcl(m0[2])}}}",
        "}",
        "",
        "Destination archive mmArchive",
        "Schedule DataTable archive stage 1",
    ]
    if p["schedule_ovf"]:
        lines.append("Schedule Oxs_TimeDriver::Magnetization archive stage 1")
    lines.append("")

    # ---- Tcl 脚本（脉冲窗口 / SOT 场）----
    t_on, t_off = p["t_on"], p["t_off"]
    if use_stt:
        lines += [
            "proc Jprofile_script { total_time } {",
            f"  if {{{fmt_tcl(t_on)} < $total_time && $total_time < {fmt_tcl(t_off)}}} {{",
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
            f"  set ton {fmt_tcl(p['t_on'])}",
            f"  set toff {fmt_tcl(p['t_off'])}",
            f"  set tr {fmt_tcl(tr)}",
            f"  set hx {fmt_tcl(hx)}",
            f"  set hy {fmt_tcl(hy)}",
            f"  set hz {fmt_tcl(hz)}",
            "  # smoothstep 包络（可微，避免 evolver 步长崩溃）",
            "  if {$total_time <= $ton} {",
            "    set f 0.0",
            "    set df 0.0",
            "  } elseif {$total_time < $ton + $tr} {",
            "    set x [expr {($total_time - $ton)/$tr}]",
            "    set f [expr {$x*$x*(3-2*$x)}]",
            "    set df [expr {6*$x*(1-$x)/$tr}]",
            "  } elseif {$total_time <= $toff - $tr} {",
            "    set f 1.0",
            "    set df 0.0",
            "  } elseif {$total_time < $toff} {",
            "    set x [expr {($toff - $total_time)/$tr}]",
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
    """
    SOT 等效场（常数近似）：
      H_DL = Hmag * (sigma × m0̂)，H_FL = xi_FL * Hmag * sigma
      Hmag = hbar * theta_SH * J / (2 e mu0 Ms t)
    说明：用初始磁化方向 m0 近似 σ×m（真实 SOT 的 m 依赖需自定义扩展），
    适合零温 J 扫描的趋势研究。
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


# ---------------------------------------------------------------------------
# 组合构建（sweep 笛卡尔积 / 显式 runs）
# ---------------------------------------------------------------------------
def build_combos(cfg):
    """返回 [(run_name, merged_params), ...]。"""
    runs_spec = cfg.get("runs") or []
    sweep = cfg.get("sweep") or {}
    tag = cfg.get("tag", "")

    combos = []
    if runs_spec:
        for i, run in enumerate(runs_spec):
            params = dict(cfg)
            params.update(run.get("params", {}))
            name = run.get("name") or combo_name(list(run.get("params", {}).keys()),
                                                 run.get("params", {}), tag)
            combos.append((name, params))
    elif sweep:
        keys = list(sweep.keys())
        value_lists = [sweep[k] if isinstance(sweep[k], list) else [sweep[k]]
                       for k in keys]
        for values in itertools_product(*value_lists):
            over = dict(zip(keys, values))
            params = dict(cfg)
            params.update(over)
            name = combo_name(keys, over, tag)
            combos.append((name, params))
    else:
        name = tag or "run0"
        combos.append((name, dict(cfg)))
    return combos


def combo_name(keys, params_slice, tag=""):
    parts = [f"{k}_{fmt_val(params_slice[k])}" for k in keys]
    name = "__".join(parts) if parts else "run0"
    if tag:
        name = f"{name}__{tag}" if parts else tag
    return name


def itertools_product(*args):
    import itertools
    return itertools.product(*args)


# ---------------------------------------------------------------------------
# boxsi 执行
# ---------------------------------------------------------------------------
def find_tclsh():
    for cand in TCLSH_CANDIDATES:
        p = Path(cand)
        if p.exists() or shutil.which(cand):
            return cand
    raise RuntimeError("未找到 tclsh，请检查 ActiveTcl 安装")


def run_one_sim(name, params, outdir, dry_run=False):
    """单个组合：写 MIF/参数 -> boxsi -> 返回状态。"""
    run_dir = outdir / name
    run_dir.mkdir(parents=True, exist_ok=True)
    mif = build_mif(params)
    (run_dir / "run.mif").write_text(mif, encoding="utf-8")

    # 记录该组合完整参数
    (run_dir / "params.yaml").write_text(
        yaml.safe_dump(_clean_for_yaml(params), allow_unicode=True,
                       sort_keys=False),
        encoding="utf-8")

    if dry_run:
        return {"run": name, "status": "planned"}

    t0 = time.time()
    try:
        tclsh = find_tclsh()
        proc = subprocess.run(
            [tclsh, OOMMF_TCL, "boxsi", "run.mif"],
            cwd=str(run_dir), capture_output=True, text=True,
            timeout=float(params.get("timeout_s", 3600)))
        ok = proc.returncode == 0 and list(run_dir.glob("*.odt"))
        status = "ok" if ok else "failed"
        detail = ""
        if not ok:
            detail = (proc.stderr or proc.stdout or "").strip()[-800:]
    except subprocess.TimeoutExpired:
        status, detail = "timeout", f"超过 {params.get('timeout_s', 3600)}s"
    except Exception as exc:
        status, detail = "error", str(exc)[-800:]

    elapsed = time.time() - t0
    print(f"  [{status.upper():7s}] {name}  ({elapsed:.1f}s)")
    if detail:
        print(f"          {detail}")

    # 运行失败记录
    if status != "ok":
        (run_dir / "FAILED.txt").write_text(f"{status}\n{detail}\n",
                                            encoding="utf-8")
    return {"run": name, "status": status, "elapsed_s": round(elapsed, 1)}


def _clean_for_yaml(params):
    """把不可序列化项转成可写 YAML 的形式。"""
    out = {}
    for k, v in params.items():
        if k in ("sweep", "runs"):
            continue
        if isinstance(v, (list, tuple)):
            out[k] = list(v)
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="批量 OOMMF 微磁仿真")
    ap.add_argument("--config", default=None, help="基础配置 YAML/JSON")
    ap.add_argument("--out", default="results", help="结果输出目录")
    ap.add_argument("--sweep", action="append", default=[],
                    help='参数扫描，可多次，如 "J=4e11,6e11,8e11"')
    ap.add_argument("--runs", default=None, help="显式运行列表 YAML/JSON")
    ap.add_argument("--tag", default="", help="运行名后缀")
    ap.add_argument("--dry-run", action="store_true", help="只打印运行计划")
    ap.add_argument("--no-analyze", action="store_true", help="仿真后不自动分析")
    ap.add_argument("--workers", type=int, default=None, help="并行进程数")
    args = ap.parse_args()

    cfg = load_config(args.config)

    # CLI 扫描合并进配置
    if args.sweep:
        cli_sweep = dict(cfg.get("sweep") or {})
        for s in args.sweep:
            key, vals = parse_sweep_arg(s)
            cli_sweep[key] = vals
        cfg["sweep"] = cli_sweep
    if args.runs:
        text = Path(args.runs).read_text(encoding="utf-8")
        runs = yaml.safe_load(text)
        cfg["runs"] = _coerce(runs if runs is not None else json.loads(text))
    if args.tag:
        cfg["tag"] = args.tag
    if args.workers:
        cfg["workers"] = args.workers

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    combos = build_combos(cfg)
    print("=" * 64)
    print(f"OOMMF 批量仿真  共 {len(combos)} 个组合")
    print(f"模式={cfg['mode']}  形状={cfg['shape']}  "
          f"器件={fmt_val(cfg['length']*1e9)}x{fmt_val(cfg['width']*1e9)}"
          f"x{fmt_val(cfg['thickness']*1e9)} nm  "
          f"stage={cfg['stage_count']}x{fmt_val(cfg['DeltaStage']*1e9)}ns")
    if cfg.get("sweep"):
        print("扫描参数: " + ", ".join(f"{k}={cfg['sweep'][k]}" for k in cfg["sweep"]))
    print("-" * 64)
    for name, _ in combos:
        print(f"  -> {name}")
    print("=" * 64)

    # 记录批次配置与运行计划
    (outdir / "batch_config.yaml").write_text(
        yaml.safe_dump(_clean_for_yaml(cfg), allow_unicode=True, sort_keys=False),
        encoding="utf-8")
    with open(outdir / "run_plan.csv", "w", newline="", encoding="utf-8") as f:
        import csv
        w = csv.writer(f)
        w.writerow(["run"] + list(combos[0][1].keys()))
        for name, params in combos:
            w.writerow([name] + [params.get(k) for k in combos[0][1].keys()])

    if args.dry_run:
        print("dry-run：不执行仿真。")
        return

    # 执行
    results = []
    workers = max(1, int(cfg.get("workers", 1)))
    if workers == 1:
        for name, params in combos:
            results.append(run_one_sim(name, params, outdir))
    else:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(run_one_sim, name, params, outdir): name
                    for name, params in combos}
            for fut in as_completed(futs):
                results.append(fut.result())

    n_ok = sum(1 for r in results if r["status"] == "ok")
    print(f"\n完成：{n_ok}/{len(combos)} 个组合成功")

    # 数据处理
    if not args.no_analyze:
        sys.path.insert(0, str(Path(__file__).parent))
        try:
            import process_results
            process_results.process(outdir)
        except Exception as exc:
            print(f"[warn] 数据处理失败: {exc}")


if __name__ == "__main__":
    main()
