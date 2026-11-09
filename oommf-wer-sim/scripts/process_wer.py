#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oommf-wer-sim 数据处理脚本。

对 wer_sim.py 的输出目录做：
  1. 读取每个组合的 summary.json（WER(t) 轨迹、最终 WER、权重分辨率）
  2. 出图：
     - wer_vs_time.png      每个组合的 WER(t) 曲线（log-y）
     - wer_vs_param.png     单参数扫描（如 J）时最终 WER vs 参数
     - trajectories.png     活跃轨迹 mz(t) 叠加（抽样）
     - population.png       活跃轨迹数随时间（检验分裂效果）
  3. 汇总 summary_all.csv / summary_all.json
  4. 生成 report.md

独立用法：
  python process_wer.py --out wer_results
"""

import argparse
import csv
import json
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# 读取
# ---------------------------------------------------------------------------
def load_batch(outdir):
    outdir = Path(outdir)
    combos = []
    for d in sorted(outdir.iterdir()):
        if not d.is_dir():
            continue
        sj = d / "summary.json"
        if not sj.exists():
            continue
        combos.append((d.name, json.loads(sj.read_text(encoding="utf-8"))))
    return combos


def swept_params(combos):
    if not combos:
        return []
    keys = set()
    for _, s in combos:
        keys.update(s.get("params", {}).keys())
    swept = []
    for k in keys:
        vals = {json.dumps(s.get("params", {}).get(k), default=str)
                for _, s in combos}
        if len(vals) > 1:
            swept.append(k)
    return swept


def _param_val(s, key):
    return s.get("params", {}).get(key)


# ---------------------------------------------------------------------------
# 出图
# ---------------------------------------------------------------------------
def _setup_plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10})
    return plt


def plot_wer_vs_time(outdir, combos):
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for name, s in combos:
        tr = s.get("wer_trace") or []
        if not tr:
            continue
        t = [r[0] * 1e9 for r in tr]
        w = [max(r[1], 1e-300) for r in tr]
        ax.plot(t, w, "-o", ms=3, lw=1.3,
                label=f"{name} ({s.get('method', '?')})")
        ax.annotate(f"WER={s['wer_final']:.2e}",
                    xy=(t[-1], w[-1]), xytext=(6, 0),
                    textcoords="offset points", fontsize=7)
    ax.set_yscale("log")
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("WER  =  P(not switched)")
    ax.set_title("WER vs time (REE / MC)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    out = outdir / "wer_vs_time.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_wer_vs_param(outdir, combos, param):
    plt = _setup_plot()
    xs, ys = [], []
    for name, s in combos:
        v = _param_val(s, param)
        if v is None:
            continue
        try:
            x = float(v) if not isinstance(v, (list, tuple)) else float(v[0])
        except (TypeError, ValueError):
            continue
        xs.append(x)
        ys.append(max(s["wer_final"], 1e-300))
    if not xs:
        plt.close("all")
        return None
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in order]
    ys = [ys[i] for i in order]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(xs, ys, "-o", color="#2563eb", lw=1.5, ms=5)
    for x, y in zip(xs, ys):
        ax.annotate(f"{y:.2e}", (x, y), textcoords="offset points",
                    xytext=(0, 7), fontsize=7, ha="center")
    ax.set_xlabel(param)
    ax.set_ylabel("WER (final)")
    ax.set_title(f"Write error rate vs {param}")
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    out = outdir / f"wer_vs_{param}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_trajectories(outdir, combos, max_lines=200):
    """抽样画每条活跃轨迹的 mz(t)（mz_trace 已在分裂时继承）。"""
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(8, 5))
    plotted = 0
    for name, s in combos[:1]:
        tj = Path(outdir) / name / "trajectories.json"
        if not tj.exists():
            continue
        trajs = json.loads(tj.read_text(encoding="utf-8"))
        for tr in trajs:
            if plotted >= max_lines:
                break
            mz = tr.get("mz_trace") or []
            if len(mz) < 2:
                continue
            t = [r[0] * 1e9 for r in mz]
            v = [r[1] for r in mz]
            color = "#dc2626" if tr["status"] == "switched" else "#2563eb"
            ax.plot(t, v, lw=0.6, alpha=0.35, color=color)
            plotted += 1
        ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.7)
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("m·axiŝ (projection)")
        ax.set_title(f"Trajectory m(t) overlay ({plotted} shown; "
                     f"red=switched, blue=unswitched)")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = outdir / "trajectories.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_population(outdir, combos):
    plt = _setup_plot()
    fig, ax = plt.subplots(figsize=(8, 4))
    for name, s in combos:
        tr = s.get("wer_trace") or []
        if not tr:
            continue
        t = [r[0] * 1e9 for r in tr]
        n = [r[2] for r in tr]
        ax.plot(t, n, "-o", ms=3, lw=1.2, label=name)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("active trajectories")
    ax.set_title("Active trajectory population (splitting effectiveness)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    out = outdir / "population.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(outdir):
    outdir = Path(outdir)
    combos = load_batch(outdir)
    if not combos:
        print("[process_wer] 没有找到任何组合的 summary.json")
        return

    swept = swept_params(combos)

    # 汇总表
    rows = []
    for name, s in combos:
        row = {"combo": name,
               "method": s.get("method"),
               "T_K": s.get("T_K"),
               "n_initial": s.get("n_initial"),
               "n_final_active": s.get("n_final_active"),
               "n_failed": s.get("n_failed"),
               "n_switched_traj": s.get("n_switched_traj"),
               "wer_final": s.get("wer_final"),
               "min_active_weight": s.get("min_active_weight"),
               "elapsed_s": s.get("elapsed_s")}
        for k, v in (s.get("params") or {}).items():
            if isinstance(v, (list, tuple)):
                v = "[" + ",".join(str(x) for x in v) + "]"
            row[f"param:{k}"] = v
        rows.append(row)

    keys = rows[0].keys() if rows else []
    with open(outdir / "summary_all.csv", "w", newline="",
              encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(keys))
        w.writeheader()
        w.writerows(rows)

    (outdir / "summary_all.json").write_text(
        json.dumps([{"combo": n, **s} for n, s in combos], indent=2,
                   ensure_ascii=False, default=str),
        encoding="utf-8")

    # 出图
    produced = []
    p = plot_wer_vs_time(outdir, combos)
    if p:
        produced.append(p.name)
    for param in swept:
        if param in ("sweep",):
            continue
        p = plot_wer_vs_param(outdir, combos, param)
        if p:
            produced.append(p.name)
    p = plot_trajectories(outdir, combos)
    if p:
        produced.append(p.name)
    p = plot_population(outdir, combos)
    if p:
        produced.append(p.name)

    # 报告
    (outdir / "report.md").write_text(
        build_report(combos, swept), encoding="utf-8")

    print(f"\n[process_wer] 完成 {len(combos)} 个组合")
    print(f"  扫描参数: {swept}")
    for name, s in combos:
        print(f"  {name}: WER = {s['wer_final']:.3e}  "
              f"(active {s.get('n_final_active')}, "
              f"min_weight {s.get('min_active_weight', 0):.2e})")
    print(f"  输出: summary_all.csv / report.md / {produced}")


def build_report(combos, swept):
    lines = ["# OOMMF WER 仿真报告", ""]
    for name, s in combos:
        lines += [
            f"## {name}",
            f"- 方法: {s.get('method')}   温度: {s.get('T_K')} K",
            f"- 初始轨迹: {s.get('n_initial')}   最终活跃: {s.get('n_final_active')}",
            f"- 已切换(按权重计): {s.get('n_switched_traj')}",
            f"- **最终 WER = {s.get('wer_final'):.3e}**  (t = t_off)",
            f"- 最小活跃权重(分辨率指示): {s.get('min_active_weight'):.2e}",
            f"- 耗时: {s.get('elapsed_s')} s",
            "",
        ]
        tr = s.get("wer_trace") or []
        if tr:
            lines.append("| t (ns) | WER | active | 累计已切换(条) |")
            lines.append("|---|---|---|---|")
            for t, w, n, sw in tr:
                lines.append(f"| {t*1e9:.3g} | {w:.3e} | {n} | {sw} |")
            lines.append("")
    lines.append("## 说明")
    lines.append("- WER = P(写入脉冲结束时仍未切换) = 1 - 累计切换权重。")
    lines.append("- REE 中权重守恒：分裂只改变采样密度，不改变 WER 估计。")
    lines.append("- MC 的可测下限 ≈ 1/n_traj；REE 下限由最小活跃权重指示。")
    lines.append("- 扫描参数: " + (", ".join(swept) if swept else "无"))
    lines.append("- 数据文件: summary_all.csv / summary_all.json / wer_trace.csv")
    lines.append("- 图: wer_vs_time.png, wer_vs_*.png, trajectories.png, population.png")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="OOMMF WER 数据处理")
    ap.add_argument("--out", default="wer_results")
    args = ap.parse_args()
    process(args.out)


if __name__ == "__main__":
    main()
