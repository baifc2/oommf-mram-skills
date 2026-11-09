#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oommf-batch-sim 数据处理脚本。

对 batch_sim.py 生成的 results 目录做：
  1. 解析每个 run 文件夹中的 .odt（OOMMF DataTable）——提取时间 / mx,my,mz /
     Stage / 电流 / 能量
  2. 计算指标：投影磁化初末值、脉冲窗口极值、切换判定、能量极值
  3. 汇总全部 run -> summary.csv / summary.json
  4. 出图：轨迹叠加、单参数依赖、双参数热图
  5. J 扫描时自动提取临界电流 Jc（粗阈值区间）

独立用法：
  python process_results.py --out results
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import yaml


# ---------------------------------------------------------------------------
# ODT 解析
# ---------------------------------------------------------------------------
def parse_columns(header):
    """解析 '# Columns:' 行，返回列名列表（去掉花括号）。"""
    tokens = re.findall(r"\{[^}]*\}|\S+", header)
    return [t.strip("{} ").strip() for t in tokens]


def parse_odt(path):
    """返回 (col_names, data_matrix)。"""
    cols = []
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if line.startswith("# Columns:"):
                    cols = parse_columns(line[len("# Columns:"):])
                continue
            if not line.strip():
                continue
            vals = line.split()
            if len(vals) == len(cols):
                rows.append([float(v) for v in vals])
    return cols, np.asarray(rows, dtype=float)


def find_col(cols, *fragments):
    """按片段匹配列名；返回列索引，找不到返回 None。"""
    for i, c in enumerate(cols):
        if all(frag in c for frag in fragments):
            return i
    return None


# ---------------------------------------------------------------------------
# 指标计算
# ---------------------------------------------------------------------------
def compute_metrics(cols, data, params):
    """从一个 run 的 ODT 数据计算指标字典。"""
    t_i = find_col(cols, "Simulation time")
    mx_i, my_i, mz_i = (find_col(cols, "::mx"), find_col(cols, "::my"),
                        find_col(cols, "::mz"))
    e_i = find_col(cols, "Total energy")
    j_i = find_col(cols, "average J")
    stage_i = find_col(cols, "::Stage")

    if t_i is None or mx_i is None:
        raise ValueError("ODT 缺少必需列（Simulation time / ::mx）")

    t = data[:, t_i]
    mx, my, mz = data[:, mx_i], data[:, my_i], data[:, mz_i]
    axis = np.asarray(params.get("axis", [1, 0, 0]), dtype=float)
    axis = axis / np.linalg.norm(axis)

    # 投影磁化：沿易轴的分量（切换判据）
    m_proj = mx * axis[0] + my * axis[1] + mz * axis[2]
    m_abs = np.sqrt(mx ** 2 + my ** 2 + mz ** 2)

    t_on, t_off = params.get("t_on", 0), params.get("t_off", 1e99)
    pulse_mask = (t >= t_on) & (t <= t_off)

    m_initial = float(m_proj[0])
    m_final = float(m_proj[-1])
    m_min = float(m_proj.min())
    m_max = float(m_proj.max())
    m_pulse_min = float(m_proj[pulse_mask].min()) if pulse_mask.any() else float("nan")
    m_pulse_max = float(m_proj[pulse_mask].max()) if pulse_mask.any() else float("nan")

    # 切换判定：投影磁化跨越零
    switched = bool(m_initial * m_final < 0)

    metrics = {
        "m_proj_initial": round(m_initial, 6),
        "m_proj_final": round(m_final, 6),
        "m_proj_min": round(m_min, 6),
        "m_proj_max": round(m_max, 6),
        "m_pulse_min": None if np.isnan(m_pulse_min) else round(m_pulse_min, 6),
        "m_pulse_max": None if np.isnan(m_pulse_max) else round(m_pulse_max, 6),
        "m_abs_final": round(float(m_abs[-1]), 6),
        "switched": switched,
    }

    if e_i is not None:
        e = data[:, e_i]
        # 能量量级可到 1e-20，round 会抹零；保留 6 位有效数字
        metrics["E_min_J"] = _sig(e.min(), 6)
        metrics["E_final_J"] = _sig(e[-1], 6)
    if j_i is not None:
        metrics["J_applied_max"] = round(float(data[:, j_i].max()), 4)
    if stage_i is not None:
        metrics["stages"] = int(data[:, stage_i].max())
    return metrics


def load_run_params(run_dir):
    p = run_dir / "params.yaml"
    if p.exists():
        return yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    return {}


def process_run(run_dir):
    """处理单个 run 文件夹 -> (run_name, params, metrics, odt_paths)。"""
    name = run_dir.name
    params = load_run_params(run_dir)
    odts = sorted(run_dir.glob("*.odt"))
    if not odts:
        return name, params, {"status": "no_data"}, []
    try:
        # 若有多个 odt（多 Destination），拼接时间列最长的那个为主表
        best = max(odts, key=lambda p: p.stat().st_size)
        cols, data = parse_odt(best)
        metrics = compute_metrics(cols, data, params)
        metrics["status"] = "ok"
        metrics["odt"] = best.name
    except Exception as exc:
        metrics = {"status": "parse_error", "error": str(exc)[-500:]}
    return name, params, metrics, odts


# ---------------------------------------------------------------------------
# 汇总与推断
# ---------------------------------------------------------------------------
def infer_swept_params(runs_meta):
    """从各 run 的 params 推断哪些参数在批次间变化。"""
    if not runs_meta:
        return []
    keys = [k for k in runs_meta[0]["params"] if k not in
            ("image", "colormap")]
    swept = []
    for k in keys:
        vals = {json.dumps(m["params"].get(k), sort_keys=True, default=str)
                for m in runs_meta}
        if len(vals) > 1:
            swept.append(k)
    return swept


def param_value(params, key):
    v = params.get(key)
    return v


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------
def _setup_plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
    return plt


def plot_trajectories(outdir, runs_meta, proj_label):
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(10, 6))
    for m in runs_meta:
        if m["metrics"].get("status") != "ok":
            continue
        odt = Path(m["odt_path"])
        cols, data = parse_odt(odt)
        t_i = find_col(cols, "Simulation time")
        mx_i, my_i, mz_i = (find_col(cols, "::mx"), find_col(cols, "::my"),
                            find_col(cols, "::mz"))
        axis = np.asarray(m["params"].get("axis", [1, 0, 0]), dtype=float)
        axis = axis / np.linalg.norm(axis)
        t = data[:, t_i] * 1e9
        mp = (data[:, mx_i] * axis[0] + data[:, my_i] * axis[1]
              + data[:, mz_i] * axis[2])
        ax.plot(t, mp, lw=1.2, label=m["name"])
        ax.plot(t[-1], mp[-1], "o", ms=4)
    ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.6)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel(proj_label)
    ax.set_title("Trajectory overlay (all runs)")
    if ax.get_legend_handles_labels()[1]:
        ax.legend(fontsize=6, ncol=2, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = outdir / "trajectory_overlay.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_dependence(outdir, runs_meta, param, proj_label):
    """单参数扫描：末态投影磁化 vs 参数值。"""
    plt = _setup_plot()
    xs, ys, sw = [], [], []
    for m in runs_meta:
        if m["metrics"].get("status") != "ok":
            continue
        v = param_value(m["params"], param)
        if v is None:
            continue
        try:
            x = float(v) if not isinstance(v, (list, tuple)) else float(v[0])
        except (TypeError, ValueError):
            continue
        xs.append(x)
        ys.append(m["metrics"]["m_proj_final"])
        sw.append(m["metrics"]["switched"])

    if not xs:
        plt.close("all")
        return None

    xs, ys, sw = map(np.asarray, (xs, ys, sw))
    order = np.argsort(xs)
    xs, ys, sw = xs[order], ys[order], sw[order]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(xs[sw], ys[sw], "o", color="#dc2626", label="switched")
    ax.plot(xs[~sw], ys[~sw], "s", color="#0e9f6e", label="not switched")
    ax.plot(xs, ys, "-", color="#2563eb", lw=1, alpha=0.5)
    ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.6)
    if max(xs) / min(xs) > 100:          # 跨度大时用对数轴（J 扫描常见）
        ax.set_xscale("log")
    ax.set_xlabel(param)
    ax.set_ylabel(f"final {proj_label}")
    ax.set_title(f"Dependence on {param}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = outdir / f"dependence_{param}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_heatmap(outdir, runs_meta, p1, p2, proj_label):
    """双参数扫描：末态投影磁化热图 + 切换热图。"""
    plt = _setup_plot()
    xs, ys, zs = [], [], []
    for m in runs_meta:
        if m["metrics"].get("status") != "ok":
            continue
        try:
            v1 = float(param_value(m["params"], p1))
            v2 = float(param_value(m["params"], p2))
        except (TypeError, ValueError):
            continue
        xs.append(v1)
        ys.append(v2)
        zs.append(m["metrics"]["m_proj_final"])
    if not xs:
        plt.close("all")
        return None

    def grid(x, y, z):
        ux, uy = np.unique(x), np.unique(y)
        g = np.full((len(uy), len(ux)), np.nan)
        for xi, yi, zi in zip(x, y, z):
            g[np.where(uy == yi)[0][0], np.where(ux == xi)[0][0]] = zi
        return ux, uy, g

    ux, uy, g = grid(xs, ys, zs)
    out1 = outdir / f"heatmap_{p1}_{p2}.png"
    out2 = outdir / f"heatmap_{p1}_{p2}_switched.png"

    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(g, origin="lower", aspect="auto", cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(ux)), [f"{v:.3g}" for v in ux], rotation=45, fontsize=7)
    ax.set_yticks(range(len(uy)), [f"{v:.3g}" for v in uy], fontsize=7)
    ax.set_xlabel(p1)
    ax.set_ylabel(p2)
    ax.set_title(f"Final {proj_label}  ({p1} x {p2})")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout()
    fig.savefig(out1, dpi=150)
    plt.close(fig)

    # 切换热图
    sw_map = np.array([[1 if (m["metrics"]["switched"]) else 0
                        for m in runs_meta
                        if m["metrics"].get("status") == "ok"]], dtype=float)
    fig, ax = plt.subplots(figsize=(7, 5.5))
    im = ax.imshow(g * 0, origin="lower", aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    # 逐点画切换标记
    for m in runs_meta:
        if m["metrics"].get("status") != "ok":
            continue
        try:
            xi = np.where(ux == float(param_value(m["params"], p1)))[0][0]
            yi = np.where(uy == float(param_value(m["params"], p2)))[0][0]
        except (IndexError, TypeError, ValueError):
            continue
        ax.plot(xi, yi,
                "o" if m["metrics"]["switched"] else "x",
                color="#facc15" if m["metrics"]["switched"] else "#1e293b",
                ms=9, mew=1.5)
    ax.set_xticks(range(len(ux)), [f"{v:.3g}" for v in ux], rotation=45, fontsize=7)
    ax.set_yticks(range(len(uy)), [f"{v:.3g}" for v in uy], fontsize=7)
    ax.set_xlabel(p1)
    ax.set_ylabel(p2)
    ax.set_title(f"Switching map (o=switched, x=not)  ({p1} x {p2})")
    fig.tight_layout()
    fig.savefig(out2, dpi=150)
    plt.close(fig)
    return out1, out2


# ---------------------------------------------------------------------------
# Jc 提取（J 扫描时）
# ---------------------------------------------------------------------------
def extract_jc(runs_meta, param="J"):
    """粗阈值 Jc：最高未切换 J 与最低切换 J 之间（几何均值）。"""
    pts = []
    for m in runs_meta:
        if m["metrics"].get("status") != "ok":
            continue
        v = param_value(m["params"], param)
        if v is None:
            continue
        try:
            x = float(v)
        except (TypeError, ValueError):
            continue
        pts.append((x, m["metrics"]["switched"], m["metrics"]["m_proj_final"]))
    if not pts:
        return None
    pts.sort()
    unsw = [p for p in pts if not p[1]]
    sw = [p for p in pts if p[1]]
    if not sw:
        return {"found": False, "note": "扫描范围内未观察到切换"}
    j_lo = unsw[-1][0] if unsw else None
    j_hi = sw[0][0]
    jc = (j_lo * j_hi) ** 0.5 if j_lo else j_hi
    return {
        "found": True,
        "param": param,
        "Jc_estimate": jc,
        "J_lo": j_lo,
        "J_hi": j_hi,
        "n_unswitched": len(unsw),
        "n_switched": len(sw),
        "note": f"Jc 位于 [{j_lo:.3e}, {j_hi:.3e}] A/m²，估计 {jc:.3e}"
                if j_lo else f"Jc <= {j_hi:.3e} A/m²",
    }


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(outdir):
    outdir = Path(outdir)
    run_dirs = sorted([d for d in outdir.iterdir() if d.is_dir()])

    runs_meta = []
    for d in run_dirs:
        name, params, metrics, odts = process_run(d)
        runs_meta.append({
            "name": name,
            "params": params,
            "metrics": metrics,
            "odt_path": str(odts[0]) if odts else "",
        })

    # 逐 run 结果 json
    for m in runs_meta:
        (outdir / m["name"] / "results.json").write_text(
            json.dumps({"metrics": m["metrics"],
                        "params": _jsonable(m["params"])},
                       indent=2, ensure_ascii=False), encoding="utf-8")

    swept = infer_swept_params(runs_meta)
    proj_label = "m (along easy axis)"

    # summary.csv / json
    flat_rows = []
    for m in runs_meta:
        row = {"run": m["name"]}
        row.update({f"param:{k}": _flat(v) for k, v in m["params"].items()
                    if k not in ("image", "colormap", "sweep", "runs")})
        row.update(m["metrics"])
        flat_rows.append(row)

    keys = flat_rows[0].keys() if flat_rows else []
    import csv
    with open(outdir / "summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(keys))
        w.writeheader()
        for r in flat_rows:
            w.writerow(r)

    summary = {
        "swept_params": swept,
        "n_runs": len(flat_rows),
        "n_switched": sum(1 for r in flat_rows if r.get("switched")),
        "runs": flat_rows,
    }
    if "J" in swept:
        summary["Jc"] = extract_jc(runs_meta, "J")
    (outdir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8")

    # 出图
    produced = []
    if len(runs_meta) > 1:
        p = plot_trajectories(outdir, runs_meta, proj_label)
        if p:
            produced.append(p.name)
        if len(swept) == 1:
            p = plot_dependence(outdir, runs_meta, swept[0], proj_label)
            if p:
                produced.append(p.name)
        elif len(swept) >= 2:
            outs = plot_heatmap(outdir, runs_meta, swept[0], swept[1], proj_label)
            if outs:
                produced.extend(o.name for o in outs)

    # 报告
    report = build_report(summary, swept, proj_label)
    (outdir / "report.md").write_text(report, encoding="utf-8")

    print(f"\n[数据处理] 完成 {len(runs_meta)} 个 run")
    print(f"  切换 {summary['n_switched']}/{summary['n_runs']}"
          f"  扫描参数: {swept}")
    if summary.get("Jc") and summary["Jc"].get("found"):
        print(f"  Jc ≈ {summary['Jc']['Jc_estimate']:.3e} A/m² "
              f"({summary['Jc']['Jc_estimate']/1e10:.2f} MA/cm²)")
    print(f"  输出: summary.csv / summary.json / report.md / {produced}")

    return summary


def build_report(summary, swept, proj_label):
    lines = ["# OOMMF 批量仿真报告", ""]
    lines.append(f"- Run 总数: {summary['n_runs']}")
    lines.append(f"- 切换数: {summary['n_switched']}")
    lines.append(f"- 扫描参数: {', '.join(swept) if swept else '无（单次仿真）'}")
    jc = summary.get("Jc")
    if jc:
        lines.append(f"- Jc: {jc.get('note', '')}")
    lines.append("")
    lines.append("| run | 扫描参数 | m_proj_initial | m_proj_final | switched | E_final (J) |")
    lines.append("|---|---|---|---|---|---|")
    for r in summary["runs"]:
        sp = " ".join(f"{k}={_flat(r.get('param:' + k, ''))}" for k in swept)
        lines.append(
            f"| {r['run']} | {sp} | {r.get('m_proj_initial', '')} | "
            f"{r.get('m_proj_final', '')} | {r.get('switched', '')} | "
            f"{r.get('E_final_J', '')} |")
    lines.append("")
    lines.append("> 数据文件: summary.csv（机器可读）、summary.json")
    lines.append("> 图: trajectory_overlay.png, dependence_*.png / heatmap_*.png")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 工具
# ---------------------------------------------------------------------------
def _sig(v, n=6):
    """保留 n 位有效数字（能量等极小量用 round 会抹零）。"""
    return float(f"{v:.{n}g}")


def _flat(v):
    if isinstance(v, (list, tuple)):
        return "[" + ",".join(str(x) for x in v) + "]"
    return v


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(x) for x in obj]
    return obj


def main():
    ap = argparse.ArgumentParser(description="OOMMF 批量仿真数据处理")
    ap.add_argument("--out", default="results", help="batch_sim 的输出目录")
    args = ap.parse_args()
    process(args.out)


if __name__ == "__main__":
    main()
