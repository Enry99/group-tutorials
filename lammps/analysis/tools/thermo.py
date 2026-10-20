#!/usr/bin/env python3
"""Thermo output of LAMMPS logs: columns vs time, averages with error bars.

Usage:
    python thermo.py 01_load.log 02_slide.log --list
    python thermo.py 01_load.log 02_slide.log -c c_temp_int v_gap "f_load[3]" -o load_slide
    python thermo.py 02_slide.log -c "-f_load[1]" c_temp_int --stats --skip 0.2 --rolling 21

What it does:
  1. reads the thermo tables of the logs, in the order given,
     minimizations left out (--minimize keeps them); a reset_timestep
     starts a new stage, a log that continues a run from a restart joins
     it (rows repeated after the restart counted once, the later copy
     kept); stages are laid end to end on one time axis (the Time column
     if printed, else Step x timestep of the log; --dt, one value per
     stage, replaces both);
  2. --list prints the columns of every run;
  3. -c: columns or expressions of columns ("-f_load[1]",
     "0.5*(f_hold_bot[1]-f_load[1])"), one panel each vs time, stage
     boundaries dashed, --rolling N adds a centred rolling mean over N
     rows; writes <out>.csv and <out>.png;
  4. --stats: for every column and stage, mean, standard deviation, error
     of the mean by blocking (Flyvbjerg-Petersen), integrated correlation
     time and number of independent samples, after leaving out the first
     --skip fraction of the stage (run-in) or before --tmin; also in
     <out>_stats.csv.
Needs numpy, pandas, matplotlib.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp_io as L  # noqa: E402
import pandas as pd  # noqa: E402

COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def stats_row(name, stage, y):
    tau = L.tau_int(y)
    return {"column": name, "stage": stage, "rows": len(y), "mean": y.mean(),
            "std": y.std(ddof=1), "sem": L.blocking_sem(y), "tau_int": tau, "N_eff": len(y) / tau}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logs", nargs="+")
    parser.add_argument("--list", action="store_true", help="print the columns and stop")
    parser.add_argument("-c", "--columns", nargs="+", default=[], help="columns or expressions")
    parser.add_argument("--dt", nargs="+", type=float, help="timestep per stage (log time units)")
    parser.add_argument("--stage", nargs="+", type=int, help="stages to keep (0, 1, ...)")
    parser.add_argument("--minimize", action="store_true",
                        help="keep the minimizations (left out by default)")
    parser.add_argument("--tmin", type=float, help="leave out rows before this time")
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--stats", action="store_true")
    parser.add_argument("--skip", type=float, default=0.2,
                        help="fraction of each stage left out of --stats (default 0.2)")
    parser.add_argument("--rolling", type=int, help="rows of the centred rolling mean")
    parser.add_argument("--units", help="for the time label (default: from the log)")
    parser.add_argument("-o", "--output", default="thermo")
    args = L.parse_args(parser)

    logs = L.expand(args.logs)
    runs = [r for log in logs for r in L.read_log(log)]
    if not runs:
        sys.exit("no thermo output in the logs")
    if args.list:
        for k, r in enumerate(runs):
            kind = "minimization" if r["minimize"] else f"timestep {r['dt']}"
            print(f"run {k} ({os.path.basename(r['log'])}, {len(r['columns']['Step'])} rows, "
                  f"{kind}): {' '.join(r['columns'])}")
        return

    d = L.join_runs(runs, args.dt, args.minimize)
    if args.stage is not None:
        d = d[d["stage"].isin(args.stage)]
    if args.tmin is not None:
        d = d[d["time"] >= args.tmin]
    if args.tmax is not None:
        d = d[d["time"] <= args.tmax]
    if d.empty:
        sys.exit("no rows left after --stage/--tmin/--tmax")
    tu = L.TIME_UNIT.get(L.log_units(runs, args.units), "time units")
    print(f"{len(d)} rows, {d['stage'].nunique()} stage(s), time {d['time'].iloc[0]:g}-"
          f"{d['time'].iloc[-1]:g} {tu}")
    if not args.columns:
        print("columns:", " ".join(c for c in d.columns if c not in ("time", "stage")))
        return

    values = {c.strip(): L.evaluate(c, d) for c in args.columns}
    out = d[["time", "stage", "Step"]].copy()
    for c, v in values.items():
        out[c] = v
    out.to_csv(args.output + ".csv", index=False)

    if args.stats:
        print(f"\n{'column':>22} {'stage':>5} {'rows':>6} {'mean':>13} {'std':>11} "
              f"{'sem (block)':>11} {'tau_int':>8} {'N_eff':>7}")
        table = []
        for c, v in values.items():
            for s in sorted(out["stage"].unique()):
                m = (out["stage"] == s).to_numpy()
                y = v[m]
                if args.tmin is None:
                    y = y[int(args.skip * len(y)):]
                if len(y) < 16:
                    print(f"{c:>22} {s:>5}  only {len(y)} rows: skipped")
                    continue
                r = stats_row(c, s, y)
                table.append(r)
                print(f"{c:>22} {s:>5} {r['rows']:>6} {r['mean']:>13.5g} {r['std']:>11.4g} "
                      f"{r['sem']:>11.4g} {r['tau_int']:>8.1f} {r['N_eff']:>7.1f}"
                      + ("   <- few independent samples" if r["N_eff"] < 20 else ""))
        print("sem: error of the mean from blocking; tau_int in rows; "
              "N_eff = rows / tau_int (below 20: run longer)")
        pd.DataFrame(table).to_csv(args.output + "_stats.csv", index=False)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    n = len(values)
    fig, axes = plt.subplots(n, 1, figsize=(7, 2.3 * n + 0.6), sharex=True, squeeze=False)
    t = out["time"].to_numpy()
    bounds = out["time"][out["stage"].diff() > 0].to_numpy()
    for ax, (c, v) in zip(axes[:, 0], values.items()):
        if args.rolling:
            ax.plot(t, v, color=COLORS[0], lw=0.8, alpha=0.3)
            ax.plot(t, L.rolling(v, args.rolling), color=COLORS[0], lw=1.6)
        else:
            ax.plot(t, v, color=COLORS[0], lw=1.0)
        for b in bounds:
            ax.axvline(b, color="0.5", ls="--", lw=0.8)
        ax.set_ylabel(c)
    axes[-1, 0].set_xlabel(f"time ({tu})")
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    print(f"-> {args.output}.csv, {args.output}.png" + (f", {args.output}_stats.csv" if args.stats else ""))


if __name__ == "__main__":
    main()
