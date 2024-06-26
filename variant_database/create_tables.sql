CREATE TABLE IF NOT EXISTS  `variant` (
    `pdb_fn` TEXT NOT NULL,
    `mutations` TEXT NOT NULL,
    `job_uuid` TEXT NOT NULL,
    `start_time` TEXT,
    `run_time` INTEGER,
    `mutate_run_time` INTEGER,
    `relax_run_time` INTEGER,
    `filter_run_time` INTEGER,
    `centroid_run_time` INTEGER,

    `total_score` REAL,
    `dslf_fa13` REAL,
    `fa_atr` REAL,
    `fa_dun` REAL,
    `fa_elec` REAL,
    `fa_intra_rep` REAL,
    `fa_intra_sol_xover4` REAL,
    `fa_rep` REAL,
    `fa_sol` REAL,
    `hbond_bb_sc` REAL,
    `hbond_lr_bb` REAL,
    `hbond_sc` REAL,
    `hbond_sr_bb` REAL,
    `lk_ball_wtd` REAL,
    `omega` REAL,
    `p_aa_pp` REAL,
    `pro_close` REAL,
    `rama_prepro` REAL,
    `ref` REAL,
    `yhh_planarity` REAL,

    `filter_total_score` REAL,
    `buried_all` REAL,
    `buried_np` REAL,
    `contact_all` REAL,
    `contact_buried_core` REAL,
    `contact_buried_core_boundary` REAL,
    `degree` REAL,
    `degree_core` REAL,
    `degree_core_boundary` REAL,
    `exposed_hydrophobics` REAL,
    `exposed_np_AFIMLWVY` REAL,
    `exposed_polars` REAL,
    `exposed_total` REAL,
    `one_core_each` REAL,
    `pack` REAL,
    `res_count_all` REAL,
    `res_count_buried_core` REAL,
    `res_count_buried_core_boundary` REAL,
    `res_count_buried_np_core` REAL,
    `res_count_buried_np_core_boundary` REAL,
    `ss_contributes_core` REAL,
    `ss_mis` REAL,
    `total_hydrophobic` REAL,
    `total_hydrophobic_AFILMVWY` REAL,
    `total_sasa` REAL,
    `two_core_each` REAL,
    `unsat_hbond` REAL,

    `centroid_total_score` REAL,
    `cbeta` REAL,
    `cenpack` REAL,
    `env` REAL,
    `hs_pair` REAL,
    `linear_chainbreak` REAL,
    `overlap_chainbreak` REAL,
    `pair` REAL,
    `rg` REAL,
    `rsigma` REAL,
    `sheet` REAL,
    `ss_pair` REAL,
    `vdw` REAL,

    PRIMARY KEY (`pdb_fn`,`mutations`,`job_uuid`),
    FOREIGN KEY (`pdb_fn`) REFERENCES pdb_file(`pdb_fn`),
    FOREIGN KEY (`job_uuid`) REFERENCES job(`uuid`));

CREATE INDEX mutations_index ON variant(mutations);
CREATE INDEX pdb_fn_index ON variant(pdb_fn);
CREATE INDEX job_uuid_index ON variant(job_uuid);

CREATE TABLE IF NOT EXISTS  `pdb_file` (
    `pdb_fn` TEXT,
    `aa_sequence` TEXT,
    `seq_len` INTEGER,
    PRIMARY KEY (`pdb_fn`));

CREATE INDEX  pdb_file_pdb_fn_index ON pdb_file(pdb_fn);

CREATE TABLE IF NOT EXISTS  `job` (
    `uuid` TEXT,
    `cluster` TEXT,
    `process` TEXT,
    `hostname` TEXT,
    `github_tag` TEXT,
    `script_start_time` TEXT,
    `hp_mutate_default_max_cycles` INTEGER,
    `hp_relax_repeats` INTEGER,
    `hp_relax_nstruct` INTEGER,
    `hp_relax_distance` REAL,
    PRIMARY KEY (`uuid`));


CREATE INDEX job_job_uuid_index ON job(uuid);
CREATE INDEX job_cluster_index ON job(cluster);
CREATE INDEX job_process_index ON job(process);
CREATE INDEX job_cluster_process_index ON job(cluster, process)
