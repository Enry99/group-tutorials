#!/usr/bin/env python3
"""Compare the runs of a timestep test (nve.in or layout.in) and say, for
each timestep, whether it is safe.

Usage:
    python check_dt.py nve_dt*.log
    python check_dt.py layout_dt*.log --ref 0.25 --warn-heat 2
    python check_dt.py nve_dt*.log --cols v_nH2 c_ts_nfree -o gas_check

What it does, for each log (dt and mass of H from the "TS_TEST" line):
  1. energy: econserve if printed (layout.in), else etotal (nve.in);
     linear fit vs time -> drift, in meV/atom/ps and as the heating it
     would give in K/ns (drift / 3 kB, the heat capacity of a solid);
     RMS of the fluctuations around the fit, sigma(E);
  2. convergence: with velocity Verlet sigma(E) grows as dt^2 (Merz and
     Shirts, PLOS ONE 13, e0202764, 2018): conv = (sigma / sigma_ref) /
     (dt / dt_ref)^2 is about 1 while dt is in that regime, larger beyond;
  3. temperatures (columns c_ts_T*): mean +- standard error (blocking);
     H and the other atoms should have the same mean;
  4. largest H step (v_dxH), free H atoms and H2 (largest number in the
     run), first peak of the H-X and H-H distributions (rdf_<tag>.txt of
     nve.in), "Dangerous builds", lost atoms and errors, speed (ns/day);
  5. verdict against the smallest dt (or --ref, marked REF). The drift of
     the reference is the floor of the model and of the precision, so the
     heating of the other runs is counted with respect to it, and only when
     larger than 3 standard errors: FAIL (error, lost atoms, dangerous
     builds, heating > --fail-heat), WARN (heating > --warn-heat, conv >
     --max-conv, temperatures or counts different from the reference),
     else PASS; writes <out>.csv and <out>.png.
The thresholds are choices, not rules: --warn-heat 5 K/ns means about 5 K
in a 1 ns run of a region without thermostat. Needs numpy, matplotlib.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

KB = 8.617333262e-5          # eV/K


def thermo_blocks(lines):
    """Thermo tables of a LAMMPS log: list of (header, rows)."""
    blocks, i = [], 0
    while i < len(lines):
        if lines[i].split()[:1] == ["Step"]:
            head, rows = lines[i].split(), []
            i += 1
            while i < len(lines):
                w = lines[i].split()
                if len(w) != len(head):
                    break
                try:
                    rows.append([float(x) for x in w])
                except ValueError:
                    break
                i += 1
            if rows:
                blocks.append((head, np.array(rows)))
        else:
            i += 1
    return blocks


def read_log(path):
    lines = open(path, errors="replace").read().splitlines()
    run = {"path": path, "errors": []}
    for ln in lines:
        m = re.match(r"TS_TEST dt (\S+) fs, mH (\S+), tag (\S+), \D*?(\d+)", ln)
        if m:
            run.update(dt=float(m.group(1)), mH=float(m.group(2)), tag=m.group(3),
                       nmd=int(m.group(4)))
        m = re.search(r"thermostat (\S+) K", ln)
        if ln.startswith("TS_TEST") and m:
            run["T"] = float(m.group(1))
        if ln.startswith("ERROR"):
            run["errors"].append(ln.strip())
        m = re.match(r"Dangerous builds = (\d+)", ln)
        if m:
            run["dangerous"] = max(run.get("dangerous", 0), int(m.group(1)))
        m = re.match(r"Performance: ([\d.eE+-]+) ns/day", ln)
        if m:
            run["ns_day"] = float(m.group(1))
    if "dt" not in run:
        sys.exit(f"{path}: no TS_TEST line (not a log of nve.in or layout.in?)")
    blocks = thermo_blocks(lines)
    if not blocks:
        sys.exit(f"{path}: no thermo output")
    head, rows = max(blocks, key=lambda b: len(b[1]))       # the run (not "run 0")
    run["data"] = {h: rows[:, k] for k, h in enumerate(head)}
    run["errors"] = list(dict.fromkeys(run["errors"]))
    return run


def block_sem(x, min_blocks=8):
    """Standard error of the mean of a correlated series (blocking)."""
    x = np.asarray(x, float)
    best = x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan
    while len(x) >= 2 * min_blocks:
        x = 0.5 * (x[: len(x) // 2 * 2:2] + x[1: len(x) // 2 * 2:2])
        best = max(best, x.std(ddof=1) / np.sqrt(len(x)))
    return best


def rdf_peaks(path, rmax_hx=1.7, rmax_hh=1.0):
    """First peak (position, height) of the H-X (X: not H) and H-H
    distributions. rdf_<tag>.txt has g and coordination of H-all and H-H;
    H-X comes from the difference of the coordinations, dn/dr / r^2 (its
    height is in arbitrary units, compared between runs only)."""
    if not os.path.exists(path):
        return {}
    rows = [ln.split() for ln in open(path) if ln.strip() and not ln.startswith("#")]
    rows = np.array([[float(v) for v in w] for w in rows if len(w) == 6])
    if len(rows) == 0:
        return {}
    r, out = rows[:, 1], {}
    n_hx = rows[:, 3] - rows[:, 5]
    g_hx = np.gradient(n_hx, r) / r ** 2
    for name, g, rmax in (("HX", g_hx, rmax_hx), ("HH", rows[:, 4], rmax_hh)):
        sel = np.where(r < rmax)[0]
        if len(sel) < 3 or g[sel].max() <= (0 if name == "HX" else 1.0):
            continue
        k = sel[np.argmax(g[sel])]
        if 0 < k < len(r) - 1:              # parabola through the three top points
            a, b, c = g[k - 1], g[k], g[k + 1]
            den = a - 2 * b + c
            shift = 0.5 * (a - c) / den if den != 0 else 0.0
        else:
            shift = 0.0
        out[name] = (r[k] + shift * (r[1] - r[0]), g[k])
    return out


def analyse(run, skip, extra_cols):
    d = run["data"]
    t = d["Time"] - d["Time"][0]
    keep = t >= skip
    ecol = "Econserve" if "Econserve" in d else "TotEng"
    t, e = t[keep], d[ecol][keep]
    fit = np.polyfit(t, e, 1)
    res = e - np.polyval(fit, t)
    # error of the slope: least squares, with the number of independent
    # samples from the blocking error of the residuals
    neff = min(len(res), (res.std(ddof=1) / block_sem(res)) ** 2)
    se = res.std(ddof=1) / np.sqrt(((t - t.mean()) ** 2).sum() * neff / len(res))
    n = run["nmd"]
    a = {"ecol": ecol, "t": t, "e": e, "drift": fit[0] / n * 1e3, "drift_err": se / n * 1e3,
         "heat": fit[0] / n * 1e3 / (3 * KB), "heat_err": se / n * 1e3 / (3 * KB),
         "sigma": res.std() * 1e3}
    if "KinEng" in d:
        a["sig_ke"] = res.std() / d["KinEng"][keep].std()
    for c in [c for c in d if c.startswith("c_ts_T")] + extra_cols:
        if c in d:
            a[c] = (d[c][keep].mean(), block_sem(d[c][keep]))
    if "v_dxH" in d:
        a["dxH"] = d["v_dxH"].max()
    for c in ("c_ts_nfree", "v_nH2"):
        if c in d:
            a[c] = d[c].max()
    tag = run.get("tag", "")
    a["rdf"] = rdf_peaks(os.path.join(os.path.dirname(run["path"]), f"rdf_{tag}.txt"))
    return a


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+")
    parser.add_argument("--ref", type=float, help="reference dt (fs), default: the smallest")
    parser.add_argument("--skip", type=float, default=0.0, help="ps left out at the start")
    parser.add_argument("--warn-heat", type=float, default=5.0, help="K/ns")
    parser.add_argument("--fail-heat", type=float, default=50.0, help="K/ns")
    parser.add_argument("--max-conv", type=float, default=2.0)
    parser.add_argument("--cols", nargs="+", default=[], help="other thermo columns to compare")
    parser.add_argument("-o", "--output", default="dt_check")
    args = parser.parse_args()

    paths = sorted({p for pat in args.logs for p in glob.glob(pat)})
    runs = sorted((read_log(p) for p in paths), key=lambda r: (r["mH"], r["dt"]))
    for r in runs:
        r["a"] = analyse(r, args.skip, args.cols)
    ref_dt = args.ref if args.ref is not None else min(r["dt"] for r in runs)
    refs = [r for r in runs if r["dt"] == ref_dt and r["mH"] == 0] or \
           [r for r in runs if r["dt"] == ref_dt]
    if not refs:
        sys.exit(f"no run with dt = {ref_dt} fs")
    ref = refs[0]["a"]
    tcols = [c for c in ref if c.startswith("c_ts_T")] + [c for c in args.cols if c in ref]

    # ---------- verdicts ----------
    for r in runs:
        a, why, fail = r["a"], [], []
        if r["errors"]:
            fail += r["errors"][:1]
        if r.get("dangerous", 0):
            fail.append(f"{r['dangerous']} dangerous builds (neighbour list: lower neigh_every)")
        # heating with respect to the reference (its drift is the floor of the
        # model and of the precision), counted only beyond 3 standard errors
        dh = a["heat"] - ref["heat"]
        real = abs(dh) > 3 * np.hypot(a["heat_err"], ref["heat_err"])
        if r is refs[0]:
            if abs(a["heat"]) > max(3 * a["heat_err"], args.warn_heat):
                why.append(f"drift of the reference itself {a['heat']:.1f} +- "
                           f"{a['heat_err']:.1f} K/ns")
        elif real and abs(dh) > args.fail_heat:
            fail.append(f"heating {dh:+.0f} K/ns more than the reference")
        elif real and abs(dh) > args.warn_heat:
            why.append(f"heating {dh:+.1f} K/ns more than the reference")
        a["conv"] = (a["sigma"] / ref["sigma"]) / (r["dt"] / ref_dt) ** 2
        if r["dt"] != ref_dt and a["conv"] > args.max_conv:
            why.append(f"sigma(E) grows faster than dt^2 (x{a['conv']:.1f})")
        for c in tcols:
            if c in a and c in ref and r is not refs[0]:
                dm = a[c][0] - ref[c][0]
                err = np.hypot(a[c][1], ref[c][1])
                if abs(dm) > max(3 * err, 5.0):
                    why.append(f"{c} {dm:+.0f} K (+-{err:.0f})")
        for c in ("c_ts_nfree", "v_nH2"):
            if c in a and c in ref and a[c] > ref[c]:
                why.append(f"{c} up to {a[c]:g} (ref {ref[c]:g})")
        for k, (pos, h) in a["rdf"].items():
            if k in ref["rdf"] and r is not refs[0]:
                p0, h0 = ref["rdf"][k]
                if abs(h / h0 - 1) > 0.05 or abs(pos - p0) > 0.005:
                    why.append(f"{k} peak {pos - p0:+.3f} A, height x{h / h0:.2f}")
        if r is refs[0]:
            r["verdict"] = "REF"
        else:
            r["verdict"] = "FAIL" if fail else ("WARN" if why else "PASS")
        r["why"] = fail + why

    # ---------- table ----------
    print(f"reference: dt = {ref_dt} fs; energy column {ref['ecol']}; "
          f"{runs[0]['nmd']} atoms that move")
    names = [c.replace("c_ts_", "") for c in tcols]
    print(f"{'dt (fs)':>7} {'mH':>5} {'drift (meV/atom/ps)':>20} {'heat (K/ns)':>14} "
          f"{'sigma(E)':>9} {'conv':>5} {'dxH':>6} {'ns/day':>7}  "
          + " ".join(f"{n:>14}" for n in names) + "  verdict")
    for r in runs:
        a = r["a"]
        temps = " ".join(f"{a[c][0]:8.1f}+-{a[c][1]:<5.1f}" if c in a else f"{'-':>14}" for c in tcols)
        print(f"{r['dt']:7.3g} {r['mH'] or '':>5} {a['drift']:10.4f}+-{a['drift_err']:<8.4f} "
              f"{a['heat']:7.1f}+-{a['heat_err']:<5.1f} {a['sigma']:9.3f} {a['conv']:5.2f} "
              f"{a.get('dxH', np.nan):6.3f} {r.get('ns_day', np.nan):7.3g}  {temps}  {r['verdict']}")
    print("sigma(E) in meV (whole system), dxH in A, temperatures in K (mean +- standard error)")
    for r in runs:
        rdf = ", ".join(f"{k} peak {p:.3f} A" for k, (p, h) in r["a"]["rdf"].items())
        counts = ", ".join(f"{c} max {r['a'][c]:g}" for c in ("c_ts_nfree", "v_nH2") if c in r["a"])
        print(f"dt {r['dt']:g}{' mH ' + str(r['mH']) if r['mH'] else ''}: "
              f"{'; '.join(r['why']) or 'no difference from the reference'}"
              + (f" | {counts}" if counts else "") + (f" | {rdf}" if rdf else ""))

    with open(args.output + ".csv", "w") as fh:
        fh.write("dt_fs,mH,drift_meV_atom_ps,drift_err,heat_K_ns,sigma_meV,conv,dxH_A,ns_day,"
                 + ",".join(f"{c},{c}_sem" for c in tcols) + ",verdict\n")
        for r in runs:
            a = r["a"]
            vals = ",".join(f"{a[c][0]:.2f},{a[c][1]:.2f}" if c in a else "," for c in tcols)
            fh.write(f"{r['dt']},{r['mH']},{a['drift']:.4f},{a['drift_err']:.4f},{a['heat']:.2f},{a['sigma']:.4f},"
                     f"{a['conv']:.3f},{a.get('dxH', np.nan):.4f},{r.get('ns_day', np.nan)},"
                     f"{vals},{r['verdict']}\n")
    plot(runs, ref_dt, refs[0], tcols, args.output + ".png")
    print(f"-> {args.output}.csv, {args.output}.png")


def plot(runs, ref_dt, ref_run, tcols, out):
    """Energy vs time, sigma(E) vs dt, mean temperatures vs dt. Runs that
    stopped with an error are left out (named in the first panel)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    ok = [r for r in runs if not r["errors"]]
    blues = plt.get_cmap("Blues")
    dts = sorted({r["dt"] for r in runs})
    color = {dt: blues(0.35 + 0.65 * k / max(len(dts) - 1, 1)) for k, dt in enumerate(dts)}
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for r in ok:
        a = r["a"]
        lab = f"{r['dt']:g} fs" + (f", mH {r['mH']:g}" if r["mH"] else "") + f" ({r['verdict']})"
        ax[0].plot(a["t"], (a["e"] - a["e"][0]) / r["nmd"] * 1e3, "--" if r["mH"] else "-",
                   color=color[r["dt"]], lw=1.2, label=lab)
    failed = [f"{r['dt']:g} fs ({r['errors'][0].replace('ERROR: ', '').split(':')[0]})"
              for r in runs if r["errors"]]
    ax[0].set(xlabel="time (ps)", ylabel=f"{ref_run['a']['ecol']} - start (meV/atom)",
              title="not shown: " + "; ".join(failed) if failed else "")
    ax[0].title.set_fontsize(8)
    ax[0].legend(fontsize=8, title="timestep", title_fontsize=8)

    plain = [r for r in ok if not r["mH"]]
    x = np.array([r["dt"] for r in plain])
    ax[1].loglog(x, [r["a"]["sigma"] for r in plain], "o", color=blues(0.9), ms=7, label="runs")
    xx = np.linspace(x.min(), x.max(), 50)
    ax[1].loglog(xx, ref_run["a"]["sigma"] * (xx / ref_dt) ** 2, color="0.5", lw=1,
                 label=r"$\propto dt^2$ from the reference")
    for r in ok:
        if r["mH"]:
            ax[1].loglog(r["dt"], r["a"]["sigma"], "s", mfc="none", color=blues(0.9), ms=8,
                         label=f"mH {r['mH']:g}")
    ax[1].set_xticks(x, [f"{v:g}" for v in x])
    ax[1].xaxis.set_minor_formatter(NullFormatter())
    ax[1].set(xlabel="timestep (fs)", ylabel="fluctuation of the energy, sigma (meV)")
    ax[1].legend(fontsize=8)

    # temperatures of the region that is not thermostatted: all, H, other atoms
    show = [c for c in tcols if not c.startswith("c_ts_T_")][:3]
    for k, c in enumerate(show):
        rs = [r for r in plain if c in r["a"]]
        ax[2].errorbar([r["dt"] + 0.012 * (k - 1) for r in rs], [r["a"][c][0] for r in rs],
                       [r["a"][c][1] for r in rs], fmt="o-", ms=5, capsize=3, lw=1.2,
                       color=["#2a78d6", "#eb6834", "#1baf7a"][k], label=c.replace("c_ts_", ""))
    if "T" in ref_run:
        ax[2].axhline(ref_run["T"], color="0.5", lw=1, ls="--", label="thermostat")
    ax[2].set_xticks(x, [f"{v:g}" for v in x])
    ax[2].set(xlabel="timestep (fs)", ylabel="mean temperature (K)")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
