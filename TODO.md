# Open points (draft, October 2026)

## Next parts of the tutorial

- [ ] MLIP part: DFT single points for datasets (removed from this draft,
      to be moved there), parameter presets per dataset, training, active
      learning.
- [x] Structure generation (`structures/`, draft): slabs, coherent
      interfaces (Zur-McGill), amorphous/doped starting structures,
      slab/lubricant/slab with LAMMPS groups, tips. Still open:
      reconstructions (pymatgen `ReconstructionGenerator` gave wrong diamond
      slabs), vacancies, stepped surfaces (`make_steps.py`), adsorbates
      (xsorb).
- [x] LAMMPS templates (`lammps/`, draft), checked locally with classical
      potentials, plus a timestep test for H (`09_timestep`). Still open:
      the items below.
- [ ] Enrico: run `lammps/09_timestep` with MACE on Leonardo on the edge
      cases (interface at the highest load and T, confined water/OH at
      high T, CVD gas species, mass of D at 1 fs); then set the default
      timestep of the templates from the results.
- [ ] `lammps/analysis/`: friction with block errors, velocity and density
      profiles, bond and fragment counts (from the tribology-nemd and
      ovito-analysis scripts), RDF, angles, coordination (SiN and Coimbra
      notebooks), thermo-log reader; OVITO from the `.venv/mace` env.
- [ ] Test `jobscripts/lmp_node.sh` (several GPUs) with a slab template:
      `processors * * 1` and the x, y load balance only act there.
- [ ] Enrico: upload the indent + slide templates with torque removal
      (rigid tip with `torque * off off off`, constant force instead of
      constant z), to replace or extend `lammps/06_indentation`.
- [ ] Possibly: band structure and DOS, AIMD, cohesive energy.
