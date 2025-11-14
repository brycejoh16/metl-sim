# Protocol Description

Remove: 
- FinalMinimizer 
- Transform mover (low resolution docker)
 
Add:
- Do FastRelax with constraints 
- MinMover after HighResDocker 
  - do BB minimization for 6 Å radius around the ligand
  - do SC minimization for 18 Å radius around the ligand , and 10 Å radius around the mutated residue