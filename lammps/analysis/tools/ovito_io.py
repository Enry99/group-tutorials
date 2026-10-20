"""OVITO helpers shared by structure.py, bonds.py and detached.py (not
meant to be run): loading the dumps as an explicit file list with a
frame-count check, frame window and time axis, pairwise bonds from
cutoffs.py, per-frame arrays (elements, ids, bonds).

Loading follows traj_io.py (E. Pedretti): OVITO 3.16 skips, without a
message, files with no digits at the '*' of a wildcard, so the list is
expanded here and the frames OVITO finds are compared with the
'ITEM: TIMESTEP' records.
"""
import numpy as np

import cutoffs as C
import lmp_io as L


def add_common_args(parser, last_default=None):
    parser.add_argument("dumps", nargs="+", help="text dumps, in chronological order")
    parser.add_argument("--type-names", nargs="+", help="element of type 1, 2, ... (no 'element' column)")
    parser.add_argument("--dt", nargs="+", type=float, help="timestep per stage (time unit of the dump)")
    parser.add_argument("--units", default="metal", help="for the time label (metal: ps)")
    parser.add_argument("--tmin", type=float)
    parser.add_argument("--tmax", type=float)
    parser.add_argument("--last", type=int, default=last_default, help="last N frames of the window")
    parser.add_argument("--every", type=int, default=1, help="every N-th frame")
    parser.add_argument("--cutoff", nargs="+", metavar="A-B:R", help="bond cutoff overrides (A)")
    parser.add_argument("--alternative", action="store_true",
                        help="cutoffs of ovito_cho.py for C-C, C-H, C-O, H-H (cutoffs.py)")


def load(dumps, type_names=None):
    """Pipeline on the explicit file list; returns (pipeline, steps of every
    OVITO frame, elements present in the first or the last frame)."""
    from ovito.io import import_file
    files = L.expand(dumps)
    idx = L.index_dump(files, dedup=False)
    pipeline = import_file(files if len(files) > 1 else files[0])
    n = pipeline.source.num_frames
    print(f"{len(files)} file(s), {n} frames: " + ", ".join(files))
    if n != len(idx):
        raise SystemExit(f"OVITO found {n} frames, the files hold {len(idx)} "
                         "'ITEM: TIMESTEP' records: check the files and the OVITO version")
    ptypes = pipeline.source.data.particles_.particle_types_
    if type_names:
        if all(t.name for t in ptypes.types):           # 'element' column: OVITO ids follow
            print("--type-names ignored: the dump has element names")  # first appearance
        else:
            for t in ptypes.types:
                if t.id <= len(type_names):
                    ptypes.type_by_id_(t.id).name = type_names[t.id - 1]
    els = [t.name for t in ptypes.types]
    if n > 1:                                           # elements that appear later (deposition)
        els += [t.name for t in pipeline.source.compute(n - 1).particles.particle_types.types
                if t.name not in els]
    if not all(els):
        raise SystemExit("types without element names: add 'element' to the dump "
                         "(dump_modify element ...) or give --type-names")
    return pipeline, np.array([s for _, _, s in idx]), els


def window(steps, args):
    """OVITO frames to analyse (restart overlaps dropped, the later copy
    kept, as lmp_io.clean_steps), with their time and stage."""
    keep = L.clean_steps(steps)
    sel, t = L.window(steps[keep], args.dt, args.tmin, args.tmax, args.last, args.every)
    _, stage = L.stage_time(steps[keep], args.dt or [1.0])
    return keep[sel], t[sel], stage[sel]


def time_label(args):
    return f"time ({L.TIME_UNIT.get(args.units, 'time units')})" if args.dt else "step"


def add_bonds(pipeline, elements, args):
    """Wrap + pairwise CreateBondsModifier; returns the cutoff table used."""
    from ovito.modifiers import CreateBondsModifier, WrapPeriodicImagesModifier
    table = C.table_for(elements, C.parse(args.cutoff), args.alternative)
    mod = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
    for (a, b), r in table.items():
        mod.set_pairwise_cutoff(a, b, r)
    pipeline.modifiers.append(WrapPeriodicImagesModifier())
    pipeline.modifiers.append(mod)
    return table


def elements_of(data):
    """Element name of each particle (array of str)."""
    tp = data.particles["Particle Type"]
    lut = np.full(max(t.id for t in tp.types) + 1, "", dtype=object)
    for t in tp.types:
        lut[t.id] = t.name
    return lut[np.asarray(tp)].astype(str)


_warned = []


def ids_of_frame(data):
    """Atom ids; without an 'id' column, 1..N in row order (valid only with
    dump_modify sort id and a constant atom count)."""
    ids = data.particles.identifiers
    if ids is not None:
        return np.asarray(ids)
    if not _warned:
        print("WARNING: no 'id' column: atoms matched by row order (needs dump_modify "
              "sort id and a constant atom count); add 'id' to the dump")
        _warned.append(data.particles.count)
    if data.particles.count != _warned[0]:
        raise SystemExit("the atom count changes and the dump has no 'id' column")
    return np.arange(1, data.particles.count + 1)


def bond_list(data):
    """(i, j) particle indices of the bonds, and the bond vectors r_j - r_i
    with the periodic images (A)."""
    b = data.particles.bonds
    if b is None or b.count == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, np.zeros((0, 3))
    topo = np.asarray(b.topology, dtype=np.int64)
    pos = np.asarray(data.particles.positions)
    h = np.asarray(data.cell)[:, :3]
    vec = pos[topo[:, 1]] - pos[topo[:, 0]] + np.asarray(b.pbc_vectors) @ h.T
    return topo[:, 0], topo[:, 1], vec


def area(data):
    """|a x b| of the cell (A^2), for densities per interface area."""
    h = np.asarray(data.cell)
    return float(np.linalg.norm(np.cross(h[:, 0], h[:, 1])))


def selector(spec, elements, groups=None):
    """Mask function for a '+'-joined spec of elements ('C+H'), id ranges
    ('1:500'), id lists ('1,5,9') or groups of a groups file: f(el, ids)."""
    toks = spec.split("+")
    el_toks = [t for t in toks if t in elements]
    others = [t for t in toks if t not in elements]
    id_set = L.ids_of("+".join(others), groups) if others else np.zeros(0, dtype=int)

    def f(el, ids):
        return np.isin(el, el_toks) | np.isin(ids, id_set)
    desc = " + ".join(([f"elements {', '.join(el_toks)}"] if el_toks else [])
                      + ([f"{len(id_set)} ids ({', '.join(others)})"] if others else []))
    return f, desc


def clusters(n, i, j, mask):
    """Clusters (connected by bonds) of the atoms in mask: label per atom,
    0 = the largest, -1 outside mask."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    keep = mask[i] & mask[j]
    g = coo_matrix((np.ones(keep.sum()), (i[keep], j[keep])), shape=(n, n))
    _, lab = connected_components(g, directed=False)
    out = np.full(n, -1)
    if mask.any():
        _, inv, cnt = np.unique(lab[mask], return_inverse=True, return_counts=True)
        rank = np.empty(len(cnt), dtype=int)
        rank[np.argsort(-cnt, kind="stable")] = np.arange(len(cnt))
        out[mask] = rank[inv]
    return out


def formulas(labels, el):
    """Hill formula of each cluster: {label: formula}."""
    from ase.formula import Formula
    m = labels >= 0
    out = {}
    for lab, e in zip(labels[m], el[m]):
        out.setdefault(lab, []).append(e)
    return {k: Formula.from_list(v).format("hill") for k, v in out.items()}


def plot_series(ax, t, y, color, label, rolling=None, **kw):
    """Raw series, or raw faded (no legend entry) with a centred rolling mean."""
    if rolling:
        ax.plot(t, y, color=color, lw=0.8, alpha=0.25, **kw)
        ax.plot(t, L.rolling(np.asarray(y, dtype=float), rolling), color=color, lw=1.6,
                label=label, **kw)
    else:
        ax.plot(t, y, color=color, lw=1.3, label=label, **kw)


def stage_lines(ax, t, stage):
    for b in np.asarray(t)[np.r_[False, np.diff(stage) > 0]]:
        ax.axvline(b, color="0.5", ls="--", lw=0.8)
