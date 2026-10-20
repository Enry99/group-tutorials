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
      potentials. Still open: runs with a MLIP on Leonardo, `lmp_node.sh`
      not tested, analysis of the trajectories (`lammps/analysis/`).
- [ ] Enrico: upload the indent + slide templates with torque removal
      (rigid tip with `torque * off off off`, constant force instead of
      constant z), to replace or extend `lammps/06_indentation`.
- [ ] Possibly: band structure and DOS, AIMD, cohesive energy.
