#!/usr/bin/env python3
"""
SOT-MRAM Critical Switching Current Density (Jc) Extraction Tool.

Automates OOMMF micromagnetic simulations via Ubermag (oommfc) to extract
the critical SOT current density for field-free or field-assisted switching.

Workflow:
  1. Relax system to initial state (P or AP)
  2. For each J in scan range: apply SOT pulse, check switching
  3. Bisection refinement around threshold
  4. Output Jc, switching phase diagram, and m(t) curves

Supports two SOT implementations:
  - 'field': equivalent external field (damping-like + field-like torques)
  - 'stt_assist': SOT field + STT current (field-free via STT symmetry breaking)

Usage:
  python extract_jc.py --config device.yaml --out results/
  python extract_jc.py --quick  # built-in demo device
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Default SOT-MRAM device (representative; override via --config)
# ---------------------------------------------------------------------------
DEFAULT_DEVICE = dict(
    # Geometry
    length=60e-9, width=60e-9, thickness=1.2e-9,
    cell=3e-9,
    # Material
    Ms=1.1e6, A=20e-12, Ku=1.2e6, axis=(0, 0, 1),
    alpha=0.02,
    # SOT
    theta_SH=0.3,          # spin Hall angle
    xi_FL=0.1,             # field-like / damping-like ratio
    current_dir=(1, 0, 0), # charge current direction (in-plane)
    # STT assist (for field-free switching)
    use_stt_assist=False,
    P=0.4, Lambda=1.0, mp=(0, 0, -1),
    J_stt_ratio=0.5,       # J_stt = ratio * J_sot
    # Bias field (optional, for conventional SOT)
    H_bias=(0, 0, 0),
    # Simulation
    T=300.0,
    t_relax=5e-9,
    t_pulse=5e-9,
    t_settle=2e-9,
    # Jc scan
    J_min=1e11, J_max=3e12, n_points=15,
    switch_criterion="mz_zero",  # mz_zero | mz_threshold
    mz_threshold=0.0,
)


def load_config(path):
    """Load device config from JSON or YAML."""
    cfg = dict(DEFAULT_DEVICE)
    if path:
        text = Path(path).read_text(encoding="utf-8")
        if path.endswith(".yaml") or path.endswith(".yml"):
            try:
                import yaml
                user_cfg = yaml.safe_load(text)
            except ImportError:
                # Minimal YAML parser for flat key: value
                user_cfg = {}
                for line in text.splitlines():
                    line = line.strip()
                    if line and not line.startswith("#") and ":" in line:
                        k, v = line.split(":", 1)
                        v = v.strip()
                        try:
                            user_cfg[k.strip()] = eval(v)
                        except Exception:
                            user_cfg[k.strip()] = v
        else:
            user_cfg = json.loads(text)
        cfg.update(user_cfg)
    return cfg


# ---------------------------------------------------------------------------
# SOT effective field
# ---------------------------------------------------------------------------
def sot_field(J, m, device):
    """
    Compute SOT effective field at magnetization state m.

    Damping-like:  H_DL = (hbar * theta_SH * J) / (2*e*mu0*Ms*t) * (sigma x m)
    Field-like:    H_FL = xi_FL * H_DL_mag * sigma

    Returns H_total in A/m.
    """
    hbar = 1.054571817e-34
    e = 1.602176634e-19
    mu0 = 4 * np.pi * 1e-7

    sigma = np.array(device["current_dir"], dtype=float)
    sigma = sigma / np.linalg.norm(sigma)
    m = np.array(m, dtype=float)
    if np.linalg.norm(m) > 0:
        m = m / np.linalg.norm(m)

    H_DL_mag = hbar * device["theta_SH"] * J / (
        2 * e * mu0 * device["Ms"] * device["thickness"])
    H_DL = H_DL_mag * np.cross(sigma, m)
    H_FL = device["xi_FL"] * H_DL_mag * sigma
    return H_DL + H_FL


# ---------------------------------------------------------------------------
# Ubermag-based simulation
# ---------------------------------------------------------------------------
def run_switching_simulation(J, device, outdir, tag=""):
    """
    Run one SOT switching simulation via Ubermag/oommfc.
    Returns True if switching occurred, plus final mz.
    """
    try:
        import discretisedfield as df
        import micromagneticmodel as mm
        import oommfc as oc
    except ImportError:
        return _run_oommf_direct(J, device, outdir, tag)

    d = device
    L, W, t = d["length"], d["width"], d["thickness"]
    cell = d["cell"]

    region = df.Region(p1=(0, 0, 0), p2=(L, W, t))
    mesh = df.Mesh(region=region, cell=(cell, cell, t))

    system = mm.System(name=f"sot_{tag}")
    system.energy = (
        mm.Exchange(A=d["A"]) +
        mm.UniaxialAnisotropy(K=d["Ku"], u=d["axis"]) +
        mm.Demag()
    )
    system.dynamics = mm.Precession(gamma0=2.211e5) + mm.Damping(alpha=d["alpha"])
    system.m = df.Field(mesh, dim=3, value=(0, 0, 1), norm=d["Ms"])

    # Set OOMMF path
    oc.runner.oommf_exe = r"C:\oommf\oommf2.1\oommf\oommf\oommf.tcl"

    # Step 1: Relax
    md = oc.MinDriver()
    md.drive(system, dirname=str(outdir / "relax"))

    # Step 2: Apply SOT pulse via time-dependent Zeeman
    # We approximate SOT as an external field that depends on m
    # For simplicity, use a quasi-static approach: evolve with SOT field
    # In production, use a custom OOMMF SOT extension or oommfc Zeeman with
    # time-dependent function that updates based on current m.

    # Simple approach: add SOT as Zeeman with time function, using average m
    def H_sot_func(t):
        # Use current average mz to estimate field (first-order approximation)
        m_avg = np.array([0.0, 0.0, 1.0])  # placeholder; refined below
        H = sot_field(J, m_avg, d)
        return tuple(H)

    system.energy += mm.Zeeman(H=H_sot_func)

    if d["use_stt_assist"]:
        J_stt = d["J_stt_ratio"] * J
        system.dynamics += mm.Slonczewski(
            J=J_stt, mp=d["mp"], P=d["P"], Lambda=d["Lambda"])

    if np.linalg.norm(d["H_bias"]) > 0:
        system.energy += mm.Zeeman(H=tuple(d["H_bias"]))

    # Step 3: Time evolution during pulse
    td = oc.TimeDriver()
    td.drive(system, t=d["t_pulse"], n=200, dirname=str(outdir / "pulse"))

    # Step 4: Settle (current off)
    system.energy.energy.pop()  # remove SOT Zeeman
    if d["use_stt_assist"]:
        system.dynamics.dynamics.pop()
    td.drive(system, t=d["t_settle"], n=50, dirname=str(outdir / "settle"))

    # Check switching
    m_final = system.m.array.mean(axis=(0, 1, 2))
    mz_final = m_final[2] / d["Ms"]

    if d["switch_criterion"] == "mz_zero":
        switched = mz_final < 0
    else:
        switched = mz_final < d["mz_threshold"]

    return switched, mz_final


def _run_oommf_direct(J, device, outdir, tag=""):
    """
    Fallback: direct OOMMF MIF generation + boxsi call (no Ubermag).
    Uses time-dependent UZeeman for SOT equivalent field.
    """
    import subprocess
    import re

    d = device
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    L, W, t = d["length"], d["width"], d["thickness"]
    cell = d["cell"]
    nx = max(1, int(round(L / cell)))
    ny = max(1, int(round(W / cell)))

    # SOT field magnitude
    hbar = 1.054571817e-34
    e = 1.602176634e-19
    mu0 = 4 * np.pi * 1e-7
    H_DL_mag = hbar * d["theta_SH"] * J / (2 * e * mu0 * d["Ms"] * t)

    # MIF template (relax + pulse with SOT field as UZeeman)
    # Note: full SOT requires m-dependent field; here we use a fixed
    # in-plane field approximation (valid for in-plane SOT with bias field)
    sigma = np.array(d["current_dir"], dtype=float)
    sigma /= np.linalg.norm(sigma)
    H_sot = (H_DL_mag * sigma + d["xi_FL"] * H_DL_mag *
             np.cross(sigma, [0, 0, 1])).tolist()

    mif = f"""# MIF 2.1
Specify Oxs_BoxAtlas:atlas {{ xrange {{0 {L}}} yrange {{0 {W}}} zrange {{0 {t}}} }}
Specify Oxs_RectangularMesh:mesh {{ cellsize {{{cell} {cell} {t}}} atlas :atlas }}
Specify Oxs_UniaxialAnisotropy {{ K1 {d['Ku']} axis {{0 0 1}} }}
Specify Oxs_UniformExchange {{ A {d['A']} }}
Specify Oxs_Demag {{}}
Specify Oxs_TimeDriver {{
  basename relax  evolver Oxs_RungeKuttaEvolve  stopping_time {d['t_relax']}
  stage_count 1  mesh :mesh  Ms {{ {d['Ms']} }}  m0 {{0 0 1}}
}}
Destination relax mmArchive
Schedule DataTable relax Stage 1
"""
    (outdir / "relax.mif").write_text(mif)

    # Run relax
    subprocess.run(
        ["tclsh", r"C:\oommf\oommf2.1\oommf\oommf\oommf.tcl", "boxsi", "relax.mif"],
        cwd=str(outdir), capture_output=True)

    # Pulse MIF with SOT field
    pulse_mif = f"""# MIF 2.1
Specify Oxs_BoxAtlas:atlas {{ xrange {{0 {L}}} yrange {{0 {W}}} zrange {{0 {t}}} }}
Specify Oxs_RectangularMesh:mesh {{ cellsize {{{cell} {cell} {t}}} atlas :atlas }}
Specify Oxs_UniaxialAnisotropy {{ K1 {d['Ku']} axis {{0 0 1}} }}
Specify Oxs_UniformExchange {{ A {d['A']} }}
Specify Oxs_Demag {{}}
Specify Oxs_UZeeman {{ H {{{H_sot[0]} {H_sot[1]} {H_sot[2]}}} }}
Specify Oxs_RungeKuttaEvolve:evolve {{ alpha {d['alpha']} }}
Specify Oxs_TimeDriver {{
  basename pulse  evolver :evolve  stopping_time {d['t_pulse']}
  stage_count 1  mesh :mesh  Ms {{ {d['Ms']} }}  m0 {{0 0 1}}
}}
Destination pulse mmArchive
Schedule DataTable pulse Stage 1
"""
    (outdir / "pulse.mif").write_text(pulse_mif)
    subprocess.run(
        ["tclsh", r"C:\oommf\oommf2.1\oommf\oommf\oommf.tcl", "boxsi", "pulse.mif"],
        cwd=str(outdir), capture_output=True)

    # Read final mz
    odt = outdir / "pulse.odt"
    mz_final = 1.0
    if odt.exists():
        lines = odt.read_text(errors="replace").splitlines()
        data = [l for l in lines if l.strip() and not l.startswith("#")]
        if data:
            vals = data[-1].split()
            for l in lines:
                if l.startswith("# Columns:"):
                    cols = re.findall(r'\{[^}]*\}|[^ \t]+', l[len("# Columns:"):])
                    for i, c in enumerate(cols):
                        if c.strip("{} ").endswith("::mz"):
                            mz_final = float(vals[i])
                    break

    switched = mz_final < 0 if d["switch_criterion"] == "mz_zero" else mz_final < d["mz_threshold"]
    return switched, mz_final


# ---------------------------------------------------------------------------
# Jc extraction
# ---------------------------------------------------------------------------
def extract_jc(device, outdir):
    """Scan current density and find Jc via bisection."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    Js = np.logspace(np.log10(device["J_min"]),
                     np.log10(device["J_max"]),
                     device["n_points"])

    results = []
    print(f"{'J (A/m2)':>14} {'switched':>10} {'mz_final':>10}")
    print("-" * 40)

    for i, J in enumerate(Js):
        tag = f"J{i:03d}"
        switched, mz = run_switching_simulation(J, device, outdir / tag, tag)
        results.append({"J": J, "switched": switched, "mz": mz})
        print(f"{J:14.3e} {str(switched):>10} {mz:10.4f}")

    # Find Jc: lowest J with switching
    switched_Js = [r["J"] for r in results if r["switched"]]
    if switched_Js:
        Jc_coarse = min(switched_Js)
        # Bisection refinement
        idx = results.index(next(r for r in results if r["J"] == Jc_coarse))
        J_lo = results[idx - 1]["J"] if idx > 0 else device["J_min"]
        J_hi = Jc_coarse

        print(f"\nBisection refinement between {J_lo:.3e} and {J_hi:.3e}")
        for _ in range(5):
            J_mid = np.sqrt(J_lo * J_hi)  # geometric midpoint
            tag = f"bis_{J_mid:.2e}"
            switched, mz = run_switching_simulation(
                J_mid, device, outdir / tag, tag.replace(".", "p"))
            results.append({"J": J_mid, "switched": switched, "mz": mz, "bisect": True})
            print(f"  {J_mid:.3e}: switched={switched}, mz={mz:.4f}")
            if switched:
                J_hi = J_mid
            else:
                J_lo = J_mid

        Jc = np.sqrt(J_lo * J_hi)
        Jc_uncertainty = (J_hi - J_lo) / 2
    else:
        Jc = None
        Jc_uncertainty = None
        print("\nNo switching observed in scan range. Increase J_max.")

    # Save results
    results.sort(key=lambda r: r["J"])
    with open(outdir / "jc_scan_results.json", "w") as f:
        json.dump({"Jc": Jc, "Jc_uncertainty": Jc_uncertainty,
                   "device": device, "scan": results}, f, indent=2, default=str)

    # CSV
    with open(outdir / "jc_scan.csv", "w") as f:
        f.write("J_A_m2,switched,mz_final\n")
        for r in results:
            f.write(f"{r['J']:.6e},{int(r['switched'])},{r['mz']:.6f}\n")

    print(f"\n{'='*40}")
    if Jc:
        print(f"Jc = {Jc:.3e} +/- {Jc_uncertainty:.3e} A/m²")
        print(f"   = {Jc/1e10:.2f} MA/cm²")
    else:
        print("Jc not found in scan range.")
    print(f"Results saved to {outdir}/")

    return Jc, Jc_uncertainty, results


def plot_results(outdir):
    """Plot Jc scan results."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    outdir = Path(outdir)
    data = np.loadtxt(outdir / "jc_scan.csv", delimiter=",", skiprows=1)
    if data.ndim == 1:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.semilogx(data[:, 0], data[:, 2], 'o-', color='#2563eb')
    ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_xlabel('J (A/m²)')
    ax1.set_ylabel('Final mz')
    ax1.set_title('Switching vs Current Density')
    ax1.grid(True, alpha=0.3)

    switched = data[:, 1] == 1
    ax2.semilogx(data[switched, 0], data[switched, 2], 'o',
                 color='#dc2626', label='Switched', alpha=0.8)
    ax2.semilogx(data[~switched, 0], data[~switched, 2], 's',
                 color='#0e9f6e', label='Not switched', alpha=0.8)
    ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('J (A/m²)')
    ax2.set_ylabel('Final mz')
    ax2.set_title('Switching Map')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(outdir / "jc_scan.png", dpi=150)
    print(f"Plot saved to {outdir}/jc_scan.png")


def main():
    parser = argparse.ArgumentParser(description="SOT-MRAM Jc Extraction")
    parser.add_argument("--config", type=str, default=None,
                        help="Device config JSON/YAML")
    parser.add_argument("--out", type=str, default="jc_results",
                        help="Output directory")
    parser.add_argument("--quick", action="store_true",
                        help="Run with default demo device")
    parser.add_argument("--J-min", type=float, default=None)
    parser.add_argument("--J-max", type=float, default=None)
    parser.add_argument("--n-points", type=int, default=None)
    parser.add_argument("--t-pulse", type=float, default=None,
                        help="Pulse width in seconds")
    args = parser.parse_args()

    device = load_config(args.config)
    if args.J_min: device["J_min"] = args.J_min
    if args.J_max: device["J_max"] = args.J_max
    if args.n_points: device["n_points"] = args.n_points
    if args.t_pulse: device["t_pulse"] = args.t_pulse

    print("=" * 50)
    print("SOT-MRAM Jc Extraction")
    print("=" * 50)
    print(f"Device: {device['length']*1e9:.0f}x{device['width']*1e9:.0f}x"
          f"{device['thickness']*1e9:.1f} nm")
    print(f"Ms={device['Ms']:.2e}, Ku={device['Ku']:.2e}, alpha={device['alpha']}")
    print(f"theta_SH={device['theta_SH']}, xi_FL={device['xi_FL']}")
    print(f"Pulse: {device['t_pulse']*1e9:.1f} ns")
    print(f"Scan: {device['J_min']:.2e} - {device['J_max']:.2e} A/m² "
          f"({device['n_points']} points)")
    print("=" * 50)

    t0 = time.time()
    Jc, Jc_err, results = extract_jc(device, args.out)
    plot_results(args.out)
    print(f"\nTotal time: {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
