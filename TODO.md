# Open points (draft, October 2026)

## Next parts of the tutorial

- [ ] MLIP part: DFT single points for datasets (removed from this draft,
      to be moved there), parameter presets per dataset, training, active
      learning.
- [ ] Structure-generation tutorial (code-independent, ASE/pymatgen), from
      Enrico's scripts: slabs of any crystal and termination (pymatgen
      `SlabGenerator`, removal of the extra layers, reconstructions;
      `create_slab_new.py`), supercell matrices and their conventions,
      heterogeneous interfaces with the Zur-McGill algorithm (pymatgen
      `ZSLGenerator` + `CoherentInterfaceBuilder`; `ZurAlgo.ipynb`),
      amorphous starting structures at the target density for melt-quench,
      substitutional doping and vacancies, adsorbates (xsorb).
- [ ] LAMMPS part (NEMD sliding, MLIP pair styles).
- [ ] Possibly: band structure and DOS, AIMD, cohesive energy.
