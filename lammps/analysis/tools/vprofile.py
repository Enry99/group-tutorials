#!/usr/bin/env python3
"""Velocity profile v_x(z) of a sliding run from the unwrapped positions of
a LAMMPS text dump (no velocities needed), with the shear-plane analysis;
overlay of several profiles.

Usage:
    python vprofile.py 02_slide.lammpstrj --dt 0.0005 --tmin 10
    python vprofile.py 02_slide.lammpstrj --dt 0.0005 --tmin 10 --elements O H --bin 0.5 --lag 2
    python vprofile.py --overlay 1GPa.csv 5GPa.csv --labels "1 GPa" "5 GPa" -o overlay_load

What it does (based on compute_v_profile.py and overlay_profiles.py, E. Pedretti):
  1. v_x of each atom between two frames lag apart = displacement of xu /
     time between them: the exact mean velocity over that interval, so
     the thermal noise shrinks as the lag grows; z of the sample is the
     mean of the two z (folded into the box if z is periodic). Atoms are
     matched by id ('id' column), or by row order if the dump has none
     (valid only with dump_modify sort id and a constant atom count);
  2. bins of --bin A (default 0.2), mean over all pairs of frames in the
     window (--tmin/--tmax with --dt, stages laid end to end after a
     reset_timestep; --last N: the last N frames); error per bin from the
     spread between pairs (a lower bound if the pairs are correlated);
     bins with fewer than --min-samples samples dropped;
  3. shear analysis on the bins with data (--zmin/--zmax to restrict it;
     empty bins between crystal planes are left out): tanh fit (shear
     plane z0, half-width w), Couette fit (slope over the sheared region,
     where dv/dz > 25 % of its maximum), model-free width w10-90 (z from
     10 % to 90 % of the velocity drop); LOCALISED if 2w < 25 % of the
     profile, else DISTRIBUTED, and only the matching fit is drawn; a
     warning if the tanh fit has R^2 < 0.9 (profile too noisy);
  4. writes <out>.csv (z, vx_mean, vx_stderr, n_samples, in A/time unit
     of the dump) and <out>.png (v in m/s for metal units, dv/dz in
     10^10 s^-1); --overlay: profiles of several CSVs, v / v_top, z = 0
     where v reaches 95 % of the drop, and <out>_widths.csv.
Widths below 3 bins are not resolved: halve --bin and check. Needs numpy,
pandas, matplotlib, scipy.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lmp_io as L  # noqa: E402

V_TO_MS = {"metal": 100.0, "real": 1.0e5}     # A/ps, A/fs -> m/s


# ---------------------------------------------------------------- profile analysis
def boxcar(y, w):
    if w <= 1:
        return np.asarray(y, float).copy()
    pad = w // 2
    yp = np.pad(np.asarray(y, float), pad, mode="edge")
    return np.convolve(yp, np.ones(w) / w, mode="same")[pad:pad + len(y)]


def longest_run(mask):
    """(start, stop) of the longest run of True."""
    best, i, n = (0, 0), 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            if j - i > best[1] - best[0]:
                best = (i, j)
            i = j
        else:
            i += 1
    return best


def plateaus(v):
    k = max(3, len(v) // 10)
    return float(np.mean(v[:k])), float(np.mean(v[-k:]))


def tanh_model(z, v0, dv, z0, w):
    return v0 + 0.5 * dv * (1.0 + np.tanh((z - z0) / w))


def fit_tanh(z, v):
    from scipy.optimize import curve_fit
    v_bot, v_top = plateaus(v)
    dv = v_top - v_bot
    if abs(dv) < 1e-12:
        return None
    zmid = float(np.interp(0.5 * (v_bot + v_top), *((v, z) if dv > 0 else (v[::-1], z[::-1]))))
    try:
        p, cov = curve_fit(tanh_model, z, v, p0=[v_bot, dv, zmid, 0.05 * (z[-1] - z[0])],
                           bounds=([-np.inf, -np.inf, z[0], 1e-3],
                                   [np.inf, np.inf, z[-1], z[-1] - z[0]]), maxfev=20000)
    except Exception:
        return None
    err = np.sqrt(np.diag(cov)) if np.all(np.isfinite(cov)) else np.full(4, np.nan)
    r2 = 1 - np.sum((v - tanh_model(z, *p)) ** 2) / np.sum((v - v.mean()) ** 2)
    return {"v0": p[0], "dv": p[1], "z0": p[2], "w": abs(p[3]), "z0_err": err[2], "w_err": err[3],
            "r2": r2}


def fit_couette(z, v, grad, frac=0.25):
    """Straight line over the sheared region (dv/dz, signed as the drop
    from bottom to top, above frac of its maximum)."""
    v_bot, v_top = plateaus(v)
    a = np.clip(grad * np.sign(v_top - v_bot), 0, None)
    if a.max() <= 0:
        return None
    j0, j1 = longest_run(a >= frac * a.max())
    if j1 - j0 < 4:
        return None
    slope, icpt = np.polyfit(z[j0:j1], v[j0:j1], 1)
    h = float(z[j1 - 1] - z[j0])
    return {"slope": slope, "icpt": icpt, "z_lo": z[j0], "z_hi": z[j1 - 1], "h": h,
            "frac_in_film": slope * h / (v_top - v_bot) if v_top != v_bot else np.nan}


def crossing(z, s, level):
    i = np.where((s[:-1] < level) & (s[1:] >= level))[0]
    if len(i) == 0:
        return np.nan
    i = i[0]
    return z[i] + (z[i + 1] - z[i]) * (level - s[i]) / (s[i + 1] - s[i])


def widths(z, v, top_level=0.95):
    """Model-free 10-90 % width and the z where v reaches top_level of the
    drop (tanh of half-width w: w10-90 = 2.20 w; ramp of thickness h: 0.80 h)."""
    v_bot, v_top = plateaus(v)
    if v_top == v_bot:
        return None
    s = np.maximum.accumulate((v - v_bot) / (v_top - v_bot))
    z10, z90 = crossing(z, s, 0.10), crossing(z, s, 0.90)
    return {"w10_90": z90 - z10, "z10": z10, "z90": z90, "z_top": crossing(z, s, top_level),
            "v_bot": v_bot, "v_top": v_top}


# ---------------------------------------------------------------- profile from the dump
def compute(args):
    idx = L.index_dump(L.expand(args.dumps))
    steps = np.array([s for _, _, s in idx])
    sel, t = L.window(steps, args.dt, args.tmin, args.tmax, args.last)
    pairs = [(sel[k], sel[k + args.lag]) for k in range(len(sel) - args.lag)]
    if not pairs:
        sys.exit("not enough frames in the window for the lag")
    print(f"{len(idx)} frames; window {len(sel)} frames (t = {t[sel[0]]:g}-{t[sel[-1]]:g}), "
          f"{len(pairs)} pairs with lag {args.lag}")

    cache = {}

    def frame(k):
        if k not in cache:
            fn, pos, _ = idx[k]
            cache[k] = L.read_frame(fn, pos)
            for old in [j for j in cache if j < k - args.lag]:
                del cache[old]
        return cache[k]

    f0 = frame(pairs[0][0])
    lo_z, lz, zper = f0["lo"][2], f0["hi"][2] - f0["lo"][2], f0["pbc"][2]
    nb = int(np.ceil(lz / args.bin)) + 1
    edges = lo_z + args.bin * np.arange(nb + 1)
    zc = 0.5 * (edges[:-1] + edges[1:])
    sw, swv = np.zeros(nb), np.zeros(nb)
    per_pair = np.full((len(pairs), nb), np.nan)
    has_id = "id" in f0["cols"]
    if not has_id:
        print("no 'id' column: atoms matched by row order (needs dump_modify sort id)")
    dzs = []
    for p, (ka, kb) in enumerate(pairs):
        A, B = frame(ka)["cols"], frame(kb)["cols"]
        if has_id:
            common, ia, ib = np.intersect1d(A["id"].to_numpy(), B["id"].to_numpy(),
                                            assume_unique=True, return_indices=True)
        else:
            if len(A) != len(B):
                sys.exit("atom count changes and the dump has no 'id': dump the ids")
            ia = ib = np.arange(len(A))
        dt = (t[kb] - t[ka])
        xa, xb = A["xu"].to_numpy()[ia], B["xu"].to_numpy()[ib]
        za, zb = A["zu"].to_numpy()[ia], B["zu"].to_numpy()[ib]
        vx = (xb - xa) / dt
        zm = 0.5 * (za + zb)
        if zper:
            zm = lo_z + np.mod(zm - lo_z, lz)
        m = np.ones(len(vx), bool)
        if args.elements:
            m = np.isin(A["element"].to_numpy()[ia], args.elements)
        dzs.append(np.mean(np.abs(zb - za)[m]))
        b = np.clip(np.floor((zm - lo_z) / args.bin).astype(int), 0, nb - 1)
        cnt = np.bincount(b[m], minlength=nb).astype(float)
        vs = np.bincount(b[m], weights=vx[m], minlength=nb)
        sw += cnt
        swv += vs
        with np.errstate(invalid="ignore", divide="ignore"):
            per_pair[p] = np.where(cnt > 0, vs / cnt, np.nan)
        if p == 0 and np.abs(xb - xa).max() > 0.5 * (f0["hi"][0] - f0["lo"][0]):
            print("WARNING: jumps larger than lx/2: wrapped x? dump xu yu zu")
    good = sw >= args.min_samples
    import warnings
    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)          # bins with < 2 pairs
        v = np.where(good, swv / np.where(sw > 0, sw, 1), np.nan)
        n_ok = np.sum(np.isfinite(per_pair), axis=0)
        err = np.where(n_ok > 1, np.nanstd(per_pair, axis=0, ddof=1) / np.sqrt(np.maximum(n_ok, 1)),
                       np.nan)
    err[~good] = np.nan
    print(f"mean |dz| over one interval: {np.mean(dzs):.3f} A (z blurring; keep it below the "
          "width you measure)")
    return pd.DataFrame({"z": zc, "vx_mean": v, "vx_stderr": err, "n_samples": sw})


def analyse_and_plot(df, args):
    vs = V_TO_MS.get(args.units, 1.0)
    vl = "m/s" if args.units in V_TO_MS else "A/time"
    ok = np.isfinite(df["vx_mean"].to_numpy())        # empty bins (between crystal planes) dropped
    if args.zmin is not None:
        ok &= df["z"].to_numpy() >= args.zmin
    if args.zmax is not None:
        ok &= df["z"].to_numpy() <= args.zmax
    if ok.sum() < 10:
        sys.exit("profile too sparse for the analysis: more frames, larger --bin")
    z, v = df["z"].to_numpy()[ok], df["vx_mean"].to_numpy()[ok]
    e = df["vx_stderr"].to_numpy()[ok]
    grad = np.gradient(boxcar(v, 5), z)
    tf = fit_tanh(z, v)
    cf = fit_couette(z, v, grad)
    w = widths(z, v)
    regime = "UNKNOWN"
    print(f"\nprofile z = {z[0]:.2f}-{z[-1]:.2f} A; plateaus {w['v_bot'] * vs:+.2f} and "
          f"{w['v_top'] * vs:+.2f} {vl}")
    if tf:
        regime = "LOCALISED" if 2 * tf["w"] < 0.25 * (z[-1] - z[0]) else "DISTRIBUTED"
        print(f"tanh fit: z0 = {tf['z0']:.2f} +- {tf['z0_err']:.2f} A, w = {tf['w']:.2f} +- "
              f"{tf['w_err']:.2f} A (2w = {2 * tf['w']:.2f} A), R^2 = {tf['r2']:.2f} -> {regime}")
        if tf["r2"] < 0.9:
            print("   R^2 below 0.9: the profile is too noisy for the fits and the regime; "
                  "run longer, or check with a larger --bin and --lag")
    if cf and regime != "LOCALISED":                    # a line through a step means nothing
        print(f"Couette fit: sheared region z = {cf['z_lo']:.2f}-{cf['z_hi']:.2f} A "
              f"(h = {cf['h']:.2f} A), shear rate {cf['slope'] * vs:.3g} x 10^10 s^-1"
              if args.units in V_TO_MS else f"Couette fit: slope {cf['slope']:.3g}",
              f", {100 * cf['frac_in_film']:.0f} % of the drop inside it")
    note = "   (under 3 bins: not resolved, halve --bin)" if w["w10_90"] < 3 * args.bin else ""
    if not np.isfinite(w["w10_90"]):
        note = "   (no 10 % or 90 % crossing: profile too noisy; more frames, larger --bin or --lag)"
    print(f"model-free width w10-90 = {w['w10_90']:.2f} A{note}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ax = plt.subplots(1, 2, figsize=(7.5, 5), sharey=True,
                           gridspec_kw={"width_ratios": [2.2, 1]})
    ax[0].fill_betweenx(z, (v - e) * vs, (v + e) * vs, color="#2a78d6", alpha=0.25, lw=0)
    ax[0].plot(v * vs, z, color="#2a78d6", lw=1.8)
    ax[1].plot(grad * vs, z, color="#eb6834", lw=1.5)
    if regime == "LOCALISED":
        for a in ax:
            a.axhspan(tf["z0"] - tf["w"], tf["z0"] + tf["w"], color="0.6", alpha=0.2, lw=0)
        ax[0].plot(tanh_model(z, tf["v0"], tf["dv"], tf["z0"], tf["w"]) * vs, z, "k-", lw=1.2,
                   label=f"tanh fit, $w_{{10-90}}$ = {w['w10_90']:.1f} Å")
    elif regime == "DISTRIBUTED" and cf:
        zz = np.linspace(cf["z_lo"], cf["z_hi"], 20)
        for a in ax:
            a.axhspan(cf["z_lo"], cf["z_hi"], color="#1baf7a", alpha=0.15, lw=0)
        ax[0].plot((cf["slope"] * zz + cf["icpt"]) * vs, zz, "--", color="#1baf7a", lw=1.5,
                   label=f"Couette fit, $w_{{10-90}}$ = {w['w10_90']:.1f} Å")
    ax[0].set(xlabel=f"$v_x$ ({vl})", ylabel="z (Å)")
    ax[1].set_xlabel(r"$dv_x/dz$ ($10^{10}$ s$^{-1}$)" if args.units in V_TO_MS else "dv/dz")
    if ax[0].get_legend_handles_labels()[0]:
        ax[0].legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)


def overlay(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    files = L.expand(args.overlay)
    labels = args.labels or [os.path.splitext(os.path.basename(f))[0] for f in files]
    cmap = plt.get_cmap("Blues")
    fig, ax = plt.subplots(figsize=(5, 5.5))
    rows = []
    for k, (f, lab) in enumerate(zip(files, labels)):
        d = pd.read_csv(f).dropna(subset=["vx_mean"])
        z, v = d["z"].to_numpy(), d["vx_mean"].to_numpy()
        w = widths(z, v)
        if w is None:
            print(f"{f}: no velocity drop, skipped")
            continue
        s = (v - w["v_bot"]) / (w["v_top"] - w["v_bot"])
        ax.plot(s, z - w["z_top"], lw=1.8, label=f"{lab} ({w['w10_90']:.1f} Å)",
                color=cmap(0.35 + 0.65 * k / max(len(files) - 1, 1)))
        rows.append({"label": lab, "csv": f, "w10_90_A": w["w10_90"], "z10": w["z10"],
                     "z90": w["z90"], "z_top": w["z_top"]})
    ax.set(xlabel=r"$(v_x - v_{bot}) / (v_{top} - v_{bot})$",
           ylabel="z − z(95 % of the drop) (Å)")
    ax.legend(fontsize=9, title="$w_{10-90}$", title_fontsize=9)
    fig.tight_layout()
    fig.savefig(args.output + ".png", dpi=200)
    pd.DataFrame(rows).to_csv(args.output + "_widths.csv", index=False)
    for r in rows:
        print(f"{r['label']:>16}: w10-90 = {r['w10_90_A']:.2f} A")
    print(f"-> {args.output}.png, {args.output}_widths.csv")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("dumps", nargs="*")
    parser.add_argument("--dt", nargs="+", type=float, help="MD timestep per stage (time unit of the dump)")
    parser.add_argument("--units", default="metal")
    parser.add_argument("--tmin", type=float)
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--last", type=int, help="last N frames of the window")
    parser.add_argument("--lag", type=int, default=1, help="frames between the two positions")
    parser.add_argument("--bin", type=float, default=0.2, help="A")
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--elements", nargs="+", help="only these elements (default all)")
    parser.add_argument("--zmin", type=float, help="analysis window in z (A)")
    parser.add_argument("--zmax", type=float)
    parser.add_argument("--overlay", nargs="+", metavar="CSV", help="overlay profiles instead")
    parser.add_argument("--labels", nargs="+")
    parser.add_argument("-o", "--output", default="vprofile")
    args = parser.parse_args()

    if args.overlay:
        overlay(args)
        return
    if not args.dumps or not args.dt:
        sys.exit("give the dump(s) and --dt (MD timestep, e.g. 0.0005 for 0.5 fs in metal units)")
    df = compute(args)
    df.to_csv(args.output + ".csv", index=False)
    analyse_and_plot(df, args)
    print(f"-> {args.output}.csv, {args.output}.png")


if __name__ == "__main__":
    main()
