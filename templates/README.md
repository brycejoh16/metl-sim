# Rosetta templates

This directory contains files needed to run Rosetta, including RosettaScripts XML files and simple text files with command line arguments.
These files are copied over to the working directory during runtime.
Some Rosetta hyperparameters are defined in these files.

## SadX FastRelax templates 

Fast Relax - energize_wd_template 


## SadX docking templates

| Name | Template dir | Rep1 vs Rep2 | Pre-dock | NStruct | Low-Res | HR Cycles | HR SC Rad (Å) | Min Mut SC (Å) | Min BB Rad (Å) | Minimizer | Metal | Lig χ |
|---|---|:---:|---|---:|:---:|---:|---:|---:|---:|---|:---:|:---:|
| Tutorial | [`vanilla/`](vanilla) | <a href="../output/htcondor_runs/vanilla/condor_energize_2025-10-17_15-26-06_sadA_rosetta_ligand_vanilla_example/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/vanilla/condor_energize_2025-10-17_15-26-06_sadA_rosetta_ligand_vanilla_example/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 200 | No | 100 | 6 | 6 | 7 | FinalMinimizer | No | Yes |
| Tutorial-Metals | [`v0_protocol_nstructs_100/`](v0_protocol_nstructs_100) | <a href="../output/htcondor_runs/v0_protocol/condor_energize_2025-10-15_13-10-15_sadA_rosetta_ligand_v0_protocol_nstructs_100/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v0_protocol/condor_energize_2025-10-15_13-10-15_sadA_rosetta_ligand_v0_protocol_nstructs_100/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 100 | No | 6 | 6 | 6 | 7 | FinalMinimizer | Yes | No |
| Tutorial-AllAtomMin | [`v0_protocol_nstructs_100_FinalMin_all_atom/`](v0_protocol_nstructs_100_FinalMin_all_atom) | <a href="../output/htcondor_runs/v0_protocol_nstructs_100_FinalMin_all_atom/condor_energize_2025-10-29_10-41-50_sadA_rosetta_ligand_v0_protocol_nstructs_100_FinalMin_all_atom/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v0_protocol_nstructs_100_FinalMin_all_atom/condor_energize_2025-10-29_10-41-50_sadA_rosetta_ligand_v0_protocol_nstructs_100_FinalMin_all_atom/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 100 | No | 100 | 6 | 12 | 15 | FinalMinimizer | Yes | No |
| Tutorial-AllAtomMin-NoMetals | [`v0_protocol_nstructs_100_FinalMin_all_atom_metals_setup_false/`](v0_protocol_nstructs_100_FinalMin_all_atom_metals_setup_false) | <a href="../output/htcondor_runs/v0_protocol_nstructs_100_FinalMin_all_atom_setup_metals_false/condor_energize_2025-10-31_12-43-53_sadA_rosetta_ligand_v0_protocol_nstructs_100_FinalMin_all_atom_setup_metals_false/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v0_protocol_nstructs_100_FinalMin_all_atom_setup_metals_false/condor_energize_2025-10-31_12-43-53_sadA_rosetta_ligand_v0_protocol_nstructs_100_FinalMin_all_atom_setup_metals_false/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 100 | No | 100 | 6 | 12 | 15 | FinalMinimizer | No | No |
| LowRes | [`v1_ligand+mutate_area/`](v1_ligand+mutate_area) | <a href="../output/htcondor_runs/v1_ligand+mutate_area/condor_energize_2025-10-28_16-51-53_sadA_rosetta_ligand_v1_ligand+mutate_area/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v1_ligand+mutate_area/condor_energize_2025-10-28_16-51-53_sadA_rosetta_ligand_v1_ligand+mutate_area/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 200 | Yes | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| LowRes-HighSampling | [`v1_ligand+mutate_area_nstruct_1000/`](v1_ligand+mutate_area_nstruct_1000) | <a href="../output/htcondor_runs/v1_ligand+mutate_area_nstruct_1000/condor_energize_2025-10-30_12-06-17_sadA_rosetta_ligand_v1_ligand+mutate_area_nstruct_1000/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v1_ligand+mutate_area_nstruct_1000/condor_energize_2025-10-30_12-06-17_sadA_rosetta_ligand_v1_ligand+mutate_area_nstruct_1000/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 1000 | Yes | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| **Final** | [`v1_ligand+mutate_area_no_low_res/`](v1_ligand+mutate_area_no_low_res) | <a href="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res/condor_energize_2025-10-30_12-30-57_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res/condor_energize_2025-10-30_12-30-57_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 200 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| Final-GlobalMin | [`v1_ligand+mutate_area_no_low_res_MinMover_all/`](v1_ligand+mutate_area_no_low_res_MinMover_all) | <a href="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res_MinMover_all/condor_energize_2025-11-04_18-27-25_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res_MinMover_all/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res_MinMover_all/condor_energize_2025-11-04_18-27-25_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res_MinMover_all/processed_run/replicate1_vs_replicate2.png" width="200"></a> | None | 200 | No | 100 | 18 | all | all | MinMover+Cα | No | Yes |
| Final-Pack | [`pack_rotamers_v1+v1_ligand+mutate_area_no_low_res/`](pack_rotamers_v1+v1_ligand+mutate_area_no_low_res) | <a href="../output/htcondor_runs/pack_rotamers_v1+v1_ligand+mutate_area_no_low_res/condor_energize_2025-11-14_16-48-48_sadA_rosetta_ligand_pack_rotamers_v1+v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/pack_rotamers_v1+v1_ligand+mutate_area_no_low_res/condor_energize_2025-11-14_16-48-48_sadA_rosetta_ligand_pack_rotamers_v1+v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png" width="200"></a> | Packing after Mut | 200 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| Final-FastRelax | [`fast_relaxv1+v1_ligand+mutate_area_no_low_res/`](fast_relaxv1+v1_ligand+mutate_area_no_low_res) | <a href="../output/htcondor_runs/fast_relaxv1+v1_ligand+mutate_area_no_low_res/condor_energize_2025-11-14_16-05-44_sadA_rosetta_ligand_fast_relaxv1+v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/fast_relaxv1+v1_ligand+mutate_area_no_low_res/condor_energize_2025-11-14_16-05-44_sadA_rosetta_ligand_fast_relaxv1+v1_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png" width="200"></a> | FastRelax after Mut | 200 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| Final-RelaxInline | [`v1_ligand+mutate_area_no_low_res_FastRelax/`](v1_ligand+mutate_area_no_low_res_FastRelax) | <a href="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res_FastRelax/condor_energize_2025-11-04_18-12-47_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res_FastRelax/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/v1_ligand+mutate_area_no_low_res_FastRelax/condor_energize_2025-11-04_18-12-47_sadA_rosetta_ligand_v1_ligand+mutate_area_no_low_res_FastRelax/processed_run/replicate1_vs_replicate2.png" width="200"></a> | FastRelax in `dock.xml` ‡ | 25 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| Final-HighCompute | [`large_wt_dock_v3_…_no_low_res/`](large_wt_dock_v3_ligand+mutate_area_no_low_res) + [`v1_…_no_low_res/`](v1_ligand+mutate_area_no_low_res)† | <a href="../output/htcondor_runs/large_dock_mutants_v3_ligand+mutate_area_no_low_res/condor_energize_2025-11-24_16-22-07_sadA_rosetta_ligand_large_dock_mutants_v3_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png"><img src="../output/htcondor_runs/large_dock_mutants_v3_ligand+mutate_area_no_low_res/condor_energize_2025-11-24_16-22-07_sadA_rosetta_ligand_large_dock_mutants_v3_ligand+mutate_area_no_low_res/processed_run/replicate1_vs_replicate2.png" width="200"></a> | Lowest WT – 500k | 200 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
| Final-Cluster | [`v1_ligand+mutate_area_no_low_res_nstruct_10_cluster/`](v1_ligand+mutate_area_no_low_res_nstruct_10_cluster) | — | None | 10 | No | 100 | 18 | 10 | 6 | MinMover+Cα | No | Yes |
† `large_wt_dock_v3_ligand+mutate_area_no_low_res` was used only to generate the
starting wild-type structure (lowest-energy of 5×10⁵ replicates of SadA wildtype); all subsequent
docking used `v1_ligand+mutate_area_no_low_res`.

‡ FastRelax is applied as a mover inside the docking XML (`repeats=5`, `hard_rep`,
design disabled), after mutation and before high-resolution docking, under Cα
coordinate restraints (sd 0.5) spanning the whole protein; restraints are cleared
before docking. In `Final-RelaxFirst` the same relax step is instead run as a
separate preceding protocol with the proper flags for a fast relax. 


§ Docking parameters are identical to **Final**; this directory differs only in
downstream analysis, where output structures were clustered locally. All
structures collapsed into a single cluster, so clustering was not used for
structure selection.



**Columns:** NStruct = docking replicates · Low-Res = low resolution docking
HR = high-resolution docking · SC Rad = side-chain repacking radius for high res docking·
Min Mut SC / Min BB Rad = mutant side-chain and backbone minimization radii ·
Metal = metal ions retained in the input structure ·
Lig χ = ligand χ-angle sampling enabled. "all" = unrestricted (whole-protein) minimization.



## SadX Final docking runs — NStruct sweep

All runs use the **Final** protocol above, with the number of docking
trajectories (NStruct) as the only varied parameter, to establish the minimum
sampling that still gives converged scores and so reduce runtime.

|           NStruct | Template dir |
|------------------:|---|
|                20 | [`v1_ligand+mutate_area_no_low_res_nstruct_20/`](v1_ligand+mutate_area_no_low_res_nstruct_20) |
|                25 | [`v1_ligand+mutate_area_no_low_res_nstruct_25/`](v1_ligand+mutate_area_no_low_res_nstruct_25) |
|                30 | [`v1_ligand+mutate_area_no_low_res_nstruct_30/`](v1_ligand+mutate_area_no_low_res_nstruct_30) |
|                35 | [`v1_ligand+mutate_area_no_low_res_nstruct_35/`](v1_ligand+mutate_area_no_low_res_nstruct_35) |
|                50 | [`v1_ligand+mutate_area_no_low_res_nstruct_50/`](v1_ligand+mutate_area_no_low_res_nstruct_50) |
|                75 | [`v1_ligand+mutate_area_no_low_res_nstruct_75/`](v1_ligand+mutate_area_no_low_res_nstruct_75) |
|  100 _(selected)_ | [`v1_ligand+mutate_area_no_low_res_nstruct_100/`](v1_ligand+mutate_area_no_low_res_nstruct_100) |
| 200 *(reference)* | [`v1_ligand+mutate_area_no_low_res/`](v1_ligand+mutate_area_no_low_res) |
