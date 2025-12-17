# Protocol Description

### Remove

- **FinalMinimizer**: Removed the generic final minimization step.

### Add

- **Coordinate Constraints**: Added an `AddConstraints` step before docking to rigidly fix metal coordinating residues (H155, G157, H246). This ensures they remain stationary during `HighResDocker`, which does not natively support specific residue exclusions in its MoveMap.
- **MinMover**: Added a tailored minimization step (named `final`) immediately after `HighResDocker`.
  - **Backbone Minimization**: Applied to residues within a **6 Å** radius around the ligand.
  - **Side-chain Minimization**: Applied to residues within an **18 Å** radius around the ligand **OR** a **10 Å** radius around the mutated residue.
  - _Note: Metal coordinating residues are explicitly excluded from all minimization and repacking steps via MoveMap exclusions (in MinMover) and Coordinate Constraints (in HighResDocker)._
