# 03 · Equilibrium lattice parameter (Birch-Murnaghan fit and variable-cell relaxation)

Two ways to the same number: Birch-Murnaghan fit of an E(V) scan, or
variable-cell relaxation. The latter is more general (non-cubic cells,
internal coordinates) and faster for large systems, but it requires a
cutoff converged for the stress (Pulay stress).

`01_eos` does the E(V) scan and fit for fcc Al, `02_vc_relax` the
variable-cell relaxation of the same crystal, and `03_pulay_fe` compares
the two for bcc Fe at three cutoffs, to show what a too-low cutoff does to
a vc-relax. Reference results are in `reference/`. Cutoff and k-grid come
from 02. The commands of each section start from this folder.

## 01_eos: E(V) scan

Compute the energy at a few lattice parameters around the expected one and
fit an equation of state E(V): its minimum gives V0 (hence a0), its
curvature the bulk modulus B0. Cubic crystals only (one parameter). The fit
is done by ASE (`tools/eos_fit.py`), which reads the `.pwo` files directly.

```bash
cd 01_eos
mkdir scan
cd scan
python ../../../tools/make_scan.py ../Al_fcc.pwi --scale 0.97 1.03 7
```

Seven inputs, `scale_0.970.pwi` ... `scale_1.030.pwi`, with the lattice
parameter `A` multiplied by the factor in the name (±3% in a, ±9% in V).

```bash
sbatch ../job.sh
python ../../../tools/eos_fit.py scale_*.pwo --bravais cF   # when the job has finished
```

`--bravais cF` (face-centered cubic) converts the volume of the primitive
cell into the cubic a0; `cI` for bcc, `cP` for simple cubic. Look at
`eos.png`: the minimum must be inside the scanned range and the points on a
smooth curve.

Reference: a0 = **4.012 Å**, B0 = 81 GPa.

## 02_vc_relax

pw.x moves atoms and cell together, driven by forces and stress, until
the pressure is below `press_conv_thr`. The input starts from a = 4.10 Å
on purpose.

```bash
cd 02_vc_relax
sbatch job.sh Al_fcc_vcrelax.pwi
```

Then read the result:

```bash
grep "bfgs converged" Al_fcc_vcrelax.pwo        # converged, how many steps
grep "new unit-cell volume" Al_fcc_vcrelax.pwo  # volume at each step, A^3 in brackets
grep "P=" Al_fcc_vcrelax.pwo | tail -1          # pressure of the final scf
```

For fcc the primitive cell holds 1/4 of the cube: a = (4V)^(1/3). Reference:
V = 16.142 Å³, a = **4.012 Å**, in 4 BFGS steps.

After the relaxation pw.x prints `Final scf calculation at the relaxed
structure.`: a new SCF with the plane waves of the final cell. Its
pressure should be close to zero.

## The two methods must agree: Pulay stress

The stress from pw.x is computed with a fixed set of plane waves; E(V) at
fixed cutoff includes the change of the basis with the volume. At a
converged cutoff the two agree; at a low cutoff they do not (Pulay
stress). For Al the effect is small:

| ecutwfc (Ry) | a0 from E(V) | a0 from vc-relax | difference |
|---|---|---|---|
| 15 | 4.0145 Å | 4.0115 Å | 0.0030 Å |
| 20 | 4.0145 Å | 4.0105 Å | 0.0040 Å |
| 40 | 4.0119 Å | 4.0118 Å | 0.0001 Å |

The final-scf pressure stayed below 0.1 kbar in all three runs: here it
does not reveal the problem. With a hard pseudopotential the effect is
much larger: next section.

## 03_pulay_fe: a hard pseudopotential (bcc Fe)

The pslibrary ultrasoft Fe (`Fe.pbe-spn-rrkjus_psl.1.0.0.UPF`, with 3s3p
semicore states) suggests 71 Ry (`00_basics/03_pseudopotentials.md`).
Both methods at 40, 60 and 90 Ry. First the low cutoff, with the inputs as
they are:

```bash
cd 03_pulay_fe
mkdir ecut40
cd ecut40
python ../../../tools/make_scan.py ../Fe_bcc.pwi --scale 0.97 1.03 7
cp ../Fe_bcc_vcrelax.pwi .
sbatch ../job.sh                 # the 7 scf and the vc-relax, one after the other
```

Then the converged cutoff in another folder (and 60 Ry if you like):

```bash
mkdir ../ecut90
cd ../ecut90
python ../../../tools/make_scan.py ../Fe_bcc.pwi --scale 0.97 1.03 7 --set ecutwfc=90 ecutrho=720
cp ../Fe_bcc_vcrelax.pwi .       # then set ecutwfc = 90, ecutrho = 720 in this copy
sbatch ../job.sh
```

`--set` changes the two cutoffs in all seven inputs. In each folder, when
the job has finished:

```bash
python ../../../tools/eos_fit.py scale_*.pwo --bravais cI
grep "new unit-cell volume" Fe_bcc_vcrelax.pwo | tail -1   # bcc: a = (2V)^(1/3)
grep "P=" Fe_bcc_vcrelax.pwo | tail -1                     # final scf
grep "total magnetization" Fe_bcc_vcrelax.pwo | tail -1
```

Reference (`reference/pulay_fe.txt`):

| ecutwfc | a0 from E(V) | a0 from vc-relax | P of the final scf | M after vc-relax |
|---|---|---|---|---|
| 40 Ry | 2.818 Å | **2.636 Å** | −2930 kbar | 0 μB |
| 60 Ry | 2.823 Å | 2.811 Å | +12 kbar | 2.25 μB |
| 90 Ry | 2.824 Å | 2.824 Å | −6 kbar | 2.27 μB |

What it shows:

- The energy converges much faster than the stress. At 40 Ry E(V) gives
  a0 within 0.006 Å (B0 20% low), but the stress at a = 2.83 Å is
  −2690 kbar instead of about −20 kbar. The vc-relax follows that stress:
  the cell shrinks by 7% and the moment disappears.
- At 60 Ry, still below the suggested 71 Ry, the vc-relax is 0.012 Å short
  and the final-scf pressure (+12 kbar) flags it.
- At 90 Ry the two methods agree.

So: vc-relax, and stress labels for MLIPs, need the cutoff converged for
the pressure (02), often above the value suggested by the pseudopotential.
Check every vc-relax with the final-scf pressure and, once per system,
against E(V).

## D3 and the lattice parameter

Same vc-relax at 40 Ry without `vdw_corr`: a0 = 4.039 Å (the usual PBE
value; experiment ~4.03 Å at low temperature, 4.05 Å at room temperature).
D3 shortens a0 of Al by 0.027 Å (0.7%). Not a reason to drop it: every
calculation you compare must use the same setting.

## Non-cubic crystals

vc-relax with `cell_dofree = 'ibrav'` keeps the Bravais lattice (e.g. hcp
Ti: a and c relax, the shape stays hexagonal). `cell_dofree = 'all'` also
lets the angles change: use it only if the symmetry itself may change.

## Common problems

- `bad line in namelist &cell` pointing at the `K_POINTS` line: the `&CELL`
  namelist is missing (required by vc-relax; `&IONS` by relax and vc-relax).
- `The maximum number of steps has been reached.`: `nstep` missing (default
  50) or too small.
- E(V) points not on a smooth curve: one SCF did not converge or converged
  to another state (magnetic systems); `tools/pwo_summary.py` on each file.
