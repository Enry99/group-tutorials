#!/usr/bin/env python3
"""Friction stress (and friction coefficient) from the thermo output of
sliding runs, with error bars from blocking.

Usage:
    python friction.py 02_slide.log --force "-f_load[1]" --data Fe_water_Fe.data
    python friction.py 02_slide.log --force "0.5*(f_hold_bot[1]-f_load[1])" \\
        --normal "f_load[3]" --data Fe_water_Fe.data --tmin 5
    python friction.py load_*/02_slide.log --force "-f_load[1]" --normal "f_load[3]" \\
        --data load_*/Fe_water_Fe.data --labels 1GPa 2GPa 5GPa -o friction_vs_load

What it does, for each log:
  1. F, the friction force (eV/A): an expression of thermo columns, the
     force from the rest of the system on the driven part (05_sliding:
     -f_load[1] on layer_top, f_hold_bot[1] on layer_bot, or their mean),
     signed so that friction is positive; the last stage of the log (the
     sliding), without the first --skip fraction or the rows before --tmin
     (acceleration, run-in);
  2. friction stress tau = F / A (GPa), A = lx * ly from --area or from the
     data file (--data: one for all logs, or one per log);
  3. mean, error of the mean by blocking, integrated correlation time and
     independent samples (as friction_vs_load.py, E. Pedretti); with
     --normal, the normal stress P = N / A and mu = tau / P;
  4. writes <out>.csv (one row per log), <out>_blocking.png (blocking
     curves: a flat right-hand tail means a converged error bar; still
     rising means the run is too short and the error a lower bound) and,
     with several logs, <out>.png: tau vs P (or vs --x).
Needs numpy, pandas, matplotlib.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp_io as L  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def data_area(path):
    """lx * ly from the header of a LAMMPS data file (A^2)."""
    lim = {}
    for line in open(path):
        w = line.split()
        if len(w) >= 4 and w[2:4] in (["xlo", "xhi"], ["ylo", "yhi"]):
            lim[w[2]] = float(w[1]) - float(w[0])
        if "Atoms" in line:
            break
    if len(lim) < 2:
        raise SystemExit(f"{path}: no xlo xhi / ylo yhi lines")
    return lim["xlo"] * lim["ylo"]


def series(log, args):
    runs = L.read_log(log)
    if not runs:
        raise SystemExit(f"{log}: no thermo output")
    d = L.join_runs(runs, args.dt)
    stage = d["stage"].max() if args.stage is None else args.stage
    d = d[d["stage"] == stage]
    d = d.assign(time=d["time"] - d["time"].iloc[0])          # from the start of the stage
    if args.tmin is not None:
        d = d[d["time"] >= args.tmin]
    else:
        d = d.iloc[int(args.skip * len(d)):]
    if len(d) < 32:
        raise SystemExit(f"{log}: only {len(d)} rows left for the averages")
    return d


def stats(y):
    return {"mean": y.mean(), "sem": L.blocking_sem(y), "tau": L.tau_int(y), "n": len(y)}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+")
    parser.add_argument("--force", required=True, help="expression of thermo columns (eV/A)")
    parser.add_argument("--normal", help="expression for the normal force (eV/A)")
    area = parser.add_mutually_exclusive_group(required=True)
    area.add_argument("--area", type=float, help="A^2")
    area.add_argument("--data", nargs="+", help="data file(s) for lx * ly")
    parser.add_argument("--stage", type=int, help="stage of each log (default: the last)")
    parser.add_argument("--tmin", type=float, help="time from the start of the stage")
    parser.add_argument("--skip", type=float, default=0.2,
                        help="fraction of the stage left out (default 0.2; not with --tmin)")
    parser.add_argument("--dt", nargs="+", type=float, help="timestep per stage, if not in the log")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("--x", nargs="+", type=float, help="abscissa of the plot (default: P)")
    parser.add_argument("--xlabel", default="normal stress P (GPa)")
    parser.add_argument("-o", "--output", default="friction")
    args = L.parse_args(parser)

    logs = L.expand(args.logs)
    datas = L.expand(args.data) if args.data else None
    if datas and len(datas) not in (1, len(logs)):
        sys.exit(f"{len(datas)} data files for {len(logs)} logs: give one, or one per log")
    labels = args.labels or [os.path.basename(os.path.dirname(os.path.abspath(f))) for f in logs]
    if len(labels) != len(logs):
        sys.exit("--labels: one per log")

    rows, curves = [], []
    for k, log in enumerate(logs):
        d = series(log, args)
        a = args.area if args.area else data_area(datas[min(k, len(datas) - 1)])
        tau = L.evaluate(args.force, d) / a * L.EV_PER_A3_TO_GPA
        s = stats(tau)
        row = {"label": labels[k], "log": log, "area_A2": a, "rows": s["n"],
               "t_from": d["time"].iloc[0], "t_to": d["time"].iloc[-1],
               "tau_GPa": s["mean"], "tau_sem": s["sem"], "tau_int_rows": s["tau"],
               "N_eff": s["n"] / s["tau"]}
        if args.normal:
            p = stats(L.evaluate(args.normal, d) / a * L.EV_PER_A3_TO_GPA)
            mu = s["mean"] / p["mean"]
            row.update(P_GPa=p["mean"], P_sem=p["sem"], mu=mu,
                       mu_err=abs(mu) * np.hypot(s["sem"] / s["mean"], p["sem"] / p["mean"]))
        rows.append(row)
        curves.append((labels[k], L.blocking_sem(tau, return_curve=True)))

    out = pd.DataFrame(rows)
    out.to_csv(args.output + ".csv", index=False)
    print(f"{'label':>14} {'rows':>5} {'tau (GPa)':>18} {'tau_int':>7} {'N_eff':>6}"
          + (f" {'P (GPa)':>16} {'mu':>15}" if args.normal else ""))
    for r in rows:
        line = (f"{r['label']:>14} {r['rows']:>5} {r['tau_GPa']:>9.4f} +- {r['tau_sem']:<6.4f}"
                f" {r['tau_int_rows']:>7.1f} {r['N_eff']:>6.0f}")
        if args.normal:
            line += (f" {r['P_GPa']:>7.3f} +- {r['P_sem']:<5.3f}"
                     f" {r['mu']:>6.3f} +- {r['mu_err']:<5.3f}")
        if r["N_eff"] < 20:
            line += "   <- few independent samples"
        print(line)
    print("errors: one standard error of the mean (blocking); rows from t = "
          + ", ".join(f"{r['t_from']:g}" for r in rows) + " of the stage")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for k, (lab, (_, sizes, sems)) in enumerate(curves):
        ax.plot(sizes, sems / sems[0], "o-", ms=4, color=COLORS[k % 4] if k < 4 else None,
                label=lab)
    ax.set_xscale("log", base=2)
    ax.set(xlabel="block size (rows)", ylabel=r"SEM / ($\sigma/\sqrt{N}$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output + "_blocking.png", dpi=200)
    outs = [args.output + ".csv", args.output + "_blocking.png"]

    if len(rows) > 1 and (args.x or args.normal):
        x = np.array(args.x) if args.x else out["P_GPa"].to_numpy()
        xerr = None if args.x else out["P_sem"].to_numpy()
        fig, ax = plt.subplots(figsize=(5.5, 4))
        ax.errorbar(x, out["tau_GPa"], yerr=out["tau_sem"], xerr=xerr, fmt="o-",
                    color=COLORS[0], capsize=3)
        for xi, yi, lab in zip(x, out["tau_GPa"], out["label"]):
            ax.annotate(lab, (xi, yi), textcoords="offset points", xytext=(5, 5), fontsize=8)
        ax.set(xlabel=args.xlabel, ylabel=r"friction stress $\tau$ (GPa)")
        fig.tight_layout()
        fig.savefig(args.output + ".png", dpi=200)
        outs.append(args.output + ".png")
    print("-> " + ", ".join(outs))


if __name__ == "__main__":
    main()
