# Protocol Description

Remove: 
- FinalMinimizer 
- Transform mover (low resolution docker)
 
Add:
- Do pack rotamers first as recommended by https://meilerlab.org/wp-content/uploads/2022/02/rosetta3_ligand_docking.pdf (step 3)
- MinMover after HighResDocker 
  - do BB minimization for 6 Å radius around the ligand
  - do SC minimization for 18 Å radius around the ligand , and 10 Å radius around the mutated residue