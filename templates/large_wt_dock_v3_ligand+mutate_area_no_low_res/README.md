# Protocol Description
This template directory serves as the home for the a huge dock run just on wildtype 
to see if that can speed up convergence by starting from a more optimal starting point. 

Remove: 
- FinalMinimizer 
 
Add:
- change cycles in high res docker to 300 from 100, this is because current runs need 3 variants per run.
This will switch that to 1 variant per job submission. 
- keep nstructs at 200, but submit 5000 runs of 200 nstructs , for a total of, 100,000 nstructs 
- MinMover after HighResDocker 
  - do BB minimization for 6 Å radius around the ligand
  - do SC minimization for 18 Å radius around the ligand , and 10 Å radius around the mutated residue