# Pseudopotentials: what they are and which to use

A pseudopotential replaces the nucleus and the core electrons with a
smooth potential, so that only valence electrons are computed and plane
waves converge at a reasonable cutoff. The choice fixes the cutoff you
need, the accuracy, and which energies you can compare.

## Three types

| Type | Name in files | Cutoff | `ecutrho` | Notes |
|---|---|---|---|---|
| Norm-conserving (NC) | `mt`, `rrkj`, `oncvpsp`, SG15, PseudoDojo | high `ecutwfc` | 4 × `ecutwfc` (the default) | simplest; required by some post-processing codes |
| Ultrasoft (US) | `van`, `rrkjus`, GBRV `uspp` | low `ecutwfc` | 8-12 × | augmentation charges need the dense grid |
| PAW | `kjpaw`, `bpaw` | low `ecutwfc` | 8-12 × | close to all-electron accuracy; can reconstruct the all-electron density |

With US and PAW, sometimes it's better to converge `ecutrho` separately from `ecutwfc`: SSSP
uses ecutrho = 8 × ecutwfc, 12 × for Fe and Mn.

## Reading a file name (QE convention)

`Element.[rel-]functional-[semicore][n]-type_author.version.UPF`

`Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`:

| Part | Meaning |
|---|---|
| `pbe` | functional the pseudo was generated with: use it with PBE |
| `s`, `p`, `d` | semicore s / p / d states in valence (here 3s and 3p: more electrons, more accurate, more expensive) |
| `n` | nonlinear core correction |
| `kjpaw` | PAW (Kresse-Joubert); `rrkjus` = ultrasoft; `van` = Vanderbilt ultrasoft |
| `_psl.1.0.0` | author/library suffix: pslibrary version 1.0.0 |
| `rel-` (prefix) | fully relativistic, only for spin-orbit runs |

Read the header of the `.UPF` file before using it: number of valence
electrons, functional and, for pslibrary, the suggested minimum cutoffs:

```bash
grep -i -m3 "z_valence\|suggested" Al.pbe-n-kjpaw_psl.1.0.0.UPF
```

The suggested cutoffs are minima, not converged values (Al: 29 Ry
suggested, 30 Ry converged in 02).

## Families you will meet

- **pslibrary** (A. Dal Corso): US and PAW for the whole periodic table,
  generated from published inputs with QE's `ld1.x`. Versions 0.x (older,
  `psl.0.1`, `psl.0.2.1`, `psl.0.3.1`) and 1.0.0. Most of the group's
  pseudos for light elements are pslibrary 1.0.0.
- **GBRV** (Garrity, Bennett, Rabe, Vanderbilt): ultrasoft, designed for
  high-throughput work at low cutoff (~40 Ry). Names like `fe_pbe_v1.5.uspp.F.UPF`.
- **Original QE library**: old Vanderbilt ultrasoft files, e.g.
  `Fe.pbe-sp-van.UPF` (Fe with 3s3p semicore, by P. Giannozzi) and the
  `*_van_ak.UPF` files. Still in use in the group for some elements: keep them where a project or dataset already uses
  them, do not introduce them in new projects without testing
  (exception: `Fe.pbe-sp-van.UPF` is the group default for Fe,
  `tools/group_pseudos.json`).
- **SG15 (ONCV)** and **PseudoDojo**: norm-conserving, well tested; needed
  when a code or method supports NC only.
- **SSSP** (Standard Solid-State Pseudopotentials): not a new family but a
  *selection*, per element, of the best tested file among the above, in
  two versions: *efficiency* (lower cutoffs) and *precision*.

## SSSP: how to use the tables

Site: <https://sssp.materialscloud.org> (PBE efficiency table:
<https://sssp.materialscloud.org/pseudopotentials/PBE/efficiency>).

For each element the site shows the chosen file, the recommended
cutoffs, and the verification tests, all computed with QE on elemental
crystals:

- **Δ-factor**: equation of state compared with all-electron (WIEN2k)
  calculations: the accuracy of the pseudo itself;
- convergence with the cutoff of **phonon frequencies** (within 1%),
  **cohesive/formation energy**, **pressure** and **band structure**;
- the recommended cutoff is where all of these are converged.

Use it to: pick a pseudo for an element the group has never used; get a
starting cutoff; check that a pseudo already in use is not an outlier.

<!-- SSSP 1.3.0 PBE efficiency, elements relevant to the group (copied from
the SSSP JSON on 2026-10-16; check the site for updates):

| El. | File | Family | ecutwfc / ecutrho (Ry) |
|---|---|---|---|
| H | `H.pbe-rrkjus_psl.1.0.0.UPF` | pslibrary 1.0.0 US | 60 / 480 |
| C | `C.pbe-n-kjpaw_psl.1.0.0.UPF` | pslibrary 1.0.0 PAW | 45 / 360 |
| N | `N.pbe-n-radius_5.UPF` | THEOS | 60 / 480 |
| O | `O.pbe-n-kjpaw_psl.0.1.UPF` | pslibrary 0.1 PAW | 50 / 400 |
| F | `f_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 | 45 / 360 |
| Al | `Al.pbe-n-kjpaw_psl.1.0.0.UPF` | pslibrary 1.0.0 PAW | 30 / 240 |
| Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` | pslibrary 1.0.0 US | 30 / 240 |
| P | `P.pbe-n-rrkjus_psl.1.0.0.UPF` | pslibrary 1.0.0 US | 30 / 240 |
| S | `s_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 | 35 / 280 |
| Ti | `ti_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 | 35 / 280 |
| Cr | `cr_pbe_v1.5.uspp.F.UPF` | GBRV 1.5 | 40 / 320 |
| Fe | `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` | pslibrary 0.2.1 PAW | 90 / 1080 |
| Ni | `ni_pbe_v1.4.uspp.F.UPF` | GBRV 1.4 | 45 / 360 |
| Cu | `Cu.paw.z_11.ld1.psl.v1.0.0-low.upf` | pslibrary 1.0.0 PAW (low) | 90 / 720 |
| Zn | `Zn_pbe_v1.uspp.F.UPF` | GBRV 1.2 | 40 / 320 |
| Mo | `Mo_ONCV_PBE-1.0.oncvpsp.upf` | SG15 (NC) | 35 / 140 |

Note how different the cutoffs are: in a multi-element system the
hardest element sets the cutoff of the whole calculation (07: H raises
the Al slab from 30-40 to 60 Ry). For Al, the convergence test of 02 gives
the same 30 Ry as SSSP. -->

<!-- ## Group rules

1. **One set per project**, listed in the project notes (and in
   `tools/group_pseudos.json` for the default set). Never mix two pseudos
   for the same element in calculations you compare.
2. Same pseudos, cutoffs, smearing, D3 (and U) in everything that enters
   one energy difference or one MLIP dataset.
3. The cutoff is set by the hardest element of the system and converged
   for the property you need.
4. New element: start from SSSP efficiency, check the UPF header, run a
   cutoff test on a simple crystal of that element, then decide with
   Enrico.
5. Never rename a `.UPF` file: the name is how others know what you used. -->
