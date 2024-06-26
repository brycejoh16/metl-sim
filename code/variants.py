""" generate variants from pdb files """
import argparse
import hashlib
import itertools
import math
import time
from os.path import join, basename, isfile
from collections import Counter
import random
from typing import Optional, Sequence, Union

import pandas as pd
from Bio.SeqIO.PdbIO import AtomIterator
from Bio.PDB import PDBParser
import numpy as np
import sqlalchemy as sqla


# silence warnings when reading PDB files generated from Rosetta (which have comments which aren't parsed by my
# approach for getting sequences from PDB files w/ Bio.SeqIO...
import warnings

import utils

warnings.filterwarnings("ignore", message="Ignoring unrecognized record ")


def gen_all_variants(base_seq, num_subs, chars, seq_idxs):
    """ generates all possible variants of base_seq with the given number of substitutions
        using the given available chars and valid sequence idxs for substitution"""
    # positions is a tuple of (pos(1), pos(2), ... pos(num_subs))
    for positions in itertools.combinations(seq_idxs, num_subs):
        # new_aas is a tuple of (aa(1), aa(2), ... aa(num_subs))
        for new_aas in itertools.product(chars, repeat=num_subs):
            if np.all([base_seq[pos] != new_aa for pos, new_aa in zip(positions, new_aas)]):
                # note the pos+1 for 1-based indexing
                variant = ",".join(["{}{}{}".format(base_seq[pos], pos+1, new_aa) for pos, new_aa in zip(positions, new_aas)])
                # should be in sorted order already, but just in case, sort it here again
                variant = utils.sort_variant_mutations(variant)
                yield variant


def gen_sample(base_seq, num_mutants, num_subs, chars, seq_idxs, rng):
    """ generates a random sample of variants with the given number of substitutions """

    # using a set and a list to maintain the order
    # this is slower and uses 2x the memory, but the final variant list will be orderd
    mutants = set()
    mutant_list = []

    for mut_num in range(num_mutants):

        found_valid_mutant = False
        while not found_valid_mutant:

            # choose the positions to mutate and sort them in ascending order
            positions = rng.choice(seq_idxs, num_subs, replace=False)
            positions.sort()

            # choose new amino acids for each of the selected positions
            subs = []
            for pos in positions:
                base_aa = base_seq[pos]
                new_aa = rng.choice([c for c in chars if c != base_aa])
                # note the pos+1 for 1-based indexing
                sub = "{}{}{}".format(base_aa, pos+1, new_aa)
                subs.append(sub)

            # generate the mutant string
            mutant = ",".join(subs)

            if mutant not in mutants:
                found_valid_mutant = True
                mutants.add(mutant)
                mutant_list.append(mutant)

    # this variant list should already be in sorted order (we sort the positions above)
    # but just in case, sort it again here. we need variants in sorted order to avoid accidental dupes.
    mutant_list = utils.sort_variant_mutations(mutant_list)

    return mutant_list


def max_possible_variants(seq_len, num_subs, num_chars):
    return math.comb(seq_len, num_subs) * ((num_chars - 1) ** num_subs)


def distribute_into_buckets(n, num_buckets, bucket_sizes):
    """ distributes n items evenly into num_buckets with sizes bucket_sizes
        i'm sure there's a better way to do this, but this works and was quick to figure out """

    # first check if the buckets have enough capacity to store all n items
    if n > sum(bucket_sizes):
        raise ValueError("buckets do not have enough capacity to store all n items")

    # all buckets start with zero items
    buckets = [0] * num_buckets

    # keep track of which buckets have free space
    free_bucket_idxs = list(range(num_buckets))

    while True:
        # remaining items left to distribute
        remaining = n - sum(buckets)

        # if we evenly split remaining number of items, how much goes into each bucket?
        amount_per_bucket = remaining // len(free_bucket_idxs)
        if amount_per_bucket == 0:
            # no more to distribute (either distributed all of them, or an odd number left)
            break

        # can the free buckets accommodate extra items?
        has_capacity = [(buckets[bi] + amount_per_bucket) <= bucket_sizes[bi] for bi in free_bucket_idxs]
        if all(has_capacity):
            # there is capacity, so distribute those items
            for bi in free_bucket_idxs:
                buckets[bi] += amount_per_bucket
        else:
            # one of the buckets doesn't have the capacity
            fbi = has_capacity.index(False)  # index into free_bucket_idxs (not buckets)
            bi = free_bucket_idxs[fbi]  # index into buckets
            buckets[bi] = bucket_sizes[bi]
            free_bucket_idxs.remove(bi)

    # distribute any stragglers
    num_stragglers = n - sum(buckets)
    free_bucket_idxs = [i for i in range(num_buckets) if buckets[i] < bucket_sizes[i]]
    # at most there will be 1 extra item to place in each free bucket so no need to loop this
    for bi in free_bucket_idxs:
        if num_stragglers > 0:
            buckets[bi] += 1
            num_stragglers -= 1

    # final sanity check
    over_capacity = [buckets[bi] > bucket_sizes[bi] for bi in range(num_buckets)]
    if any(over_capacity):
        raise AssertionError("There are buckets over capacity but that should be impossible. This function"
                             " needs to be fixed")

    return buckets


def single_pdb_local_variants(seq, target_num, num_subs_list, chars, seq_idxs, rng):
    """ generate local variants for a single PDB file.
        given the target number of variants, and the max number of substitutions,
        this function tries to generate an equal number of variants for each possible number of substitutions """

    # print out some info
    print("aa sequence: {}".format(seq))
    print("aa seq len: {}".format(len(seq)))

    # want to distribute number of target seqs evenly across range(max_subs)
    # single mutants probably not have enough possible variants
    max_variants = [max_possible_variants(len(seq_idxs), num_subs, len(chars)) for num_subs in num_subs_list]
    # print("max variants: {}".format(max_variants))

    # distribute the target_num variants to the range of substitutions
    variants_per_num_subs = distribute_into_buckets(target_num, len(num_subs_list), max_variants)

    # now generate the actual variants
    variants = []
    for num_subs, num_v, max_v in zip(num_subs_list, variants_per_num_subs, max_variants):
        # print("getting sample: {} subs, {} variants".format(num_subs, num_v))
        if num_v == max_v:
            print("num_subs: {} num_v: {} max_v: {} approach: gen all".format(num_subs, num_v, max_v))
            variants += list(gen_all_variants(seq, num_subs, chars, seq_idxs))
        elif num_v / max_v > 0.4:
            print("num_subs: {} num_v: {} max_v: {} approach: gen all, then sample".format(num_subs, num_v, max_v))
            # gen_sample could be slow if we are generating a sample approaching the max number of variants
            # in that case it would be much faster to just generate all and select a sample from the pre-generated ones
            all_variants = list(gen_all_variants(seq, num_subs, chars, seq_idxs))
            # variants += random.sample(all_variants, num_v)
            variants += rng.choice(all_variants, num_v, replace=False).tolist()
        else:
            print("num_subs: {} num_v: {} max_v: {} approach: sample".format(num_subs, num_v, max_v))
            variants += gen_sample(seq, num_v, num_subs, chars, seq_idxs, rng)

    return variants


def get_subvariants(variant, num_subs):
    # num_subs must be less than number of substitutions in the given variant
    if num_subs >= len(variant.split(",")):
        raise ValueError("num_subs must be less than the number of substitutions in the given variant ({})".format(
            len(variant.split(","))))

    sv = [",".join(muts) for muts in list(itertools.combinations(variant.split(","), num_subs))]
    # should be in sorted order if the given main variant is in sorted order, but sort here just in case
    sv = utils.sort_variant_mutations(sv)
    return sv


def load_db_variants(db_fn: str, pdb_fn: str) -> set:
    # todo: can this be sped up using connectorX?
    print("Loading existing database variants for pdb file: {}...".format(basename(pdb_fn)))
    start = time.time()
    engine = sqla.create_engine('sqlite:///{}'.format(db_fn))
    conn = engine.connect().execution_options(stream_results=True)
    query = "SELECT mutations FROM variant WHERE `pdb_fn` == \"{}\"".format(basename(pdb_fn))
    db = set(pd.read_sql_query(query, conn, coerce_float=False)["mutations"])
    conn.close()
    engine.dispose()
    print("Loaded existing database variants in {}".format(time.time() - start))
    return db


def gen_subvariants_vlist(seq: str,
                          target_num: int,
                          min_num_subs: int,
                          max_num_subs: int,
                          chars: Union[list[str], tuple[str, ...]],
                          seq_idxs: Sequence[int],
                          rng: np.random.Generator,
                          db_pdb_fn: str,
                          db_fn: Optional[str] = None):

    # max_num_subs determines the maximum number of substitutions for the main variants
    # min_num_subs determines the minimum number of substitutions for subvariants
    #  so for example, if min_num_subs is 2, then this function won't generate subvariants with 1 substitution
    # target_num is the number of variants to generate (approximate)

    # If db_fn is specified, this function will check to see if the generated variants exists in the DB already,
    # and if so, it won't return them from this function. note it only some of the subvariants are in the db,
    # then this will still return the ones that aren't in the DB.
    db = None
    if db_fn is not None:
        db = load_db_variants(db_fn, db_pdb_fn)

    # using a set and a list to maintain the order
    # this is slower and uses 2x the memory, but the final variant list will be ordered
    variants_set = set()
    variants_list = []

    while len(variants_list) < target_num:
        # generate a variant with max_num_subs substitutions
        main_v = gen_sample(seq, num_mutants=1, num_subs=max_num_subs, chars=chars, seq_idxs=seq_idxs, rng=rng)[0]

        # if the variant has already been generated, continue to next one
        if main_v in variants_set:
            continue

        # now generate all subvariants for this variant
        # generating subvariants for all number of substitutions down to single variants
        av = [main_v]
        for i in reversed(range(min_num_subs, max_num_subs)):
            av += get_subvariants(main_v, i)

        # now add this variant and all subvariants to the main list (as long as they are not already there)
        for v in av:
            # check if the variant is already in the set or already in the DB
            variant_in_set = v in variants_set
            if variant_in_set:
                print("Generated variant already in set: {}".format(v))

            variant_in_db = False
            if db_fn is not None:
                if v in db:
                    variant_in_db = True
                    print("Generated variant already in database: {}".format(v))

            if not variant_in_set and not variant_in_db:
                # only add variant to master list if it's not already in the set and it's not in the db
                variants_set.add(v)
                variants_list.append(v)

    return variants_list


def gen_subvariants_sample(db_fn: str,
                           db_pdb_fn: str,
                           target_num: int,
                           min_num_subs: int,
                           max_num_subs: int,
                           rng: np.random.Generator):

    """
    Generate a subvariants sample of an existing database...
    Will only include variants that exist in the given database,
    but will sample those variants using a subvariants approach
    """

    # load all the variants for the given pdb_fn from the database
    db_variants_set = load_db_variants(db_fn, db_pdb_fn)
    db_variants = list(db_variants_set)

    df = pd.DataFrame({"variant": db_variants, "num_mutations": [len(v.split(",")) for v in db_variants]})

    # iteratively sample variants with max_num_subs
    df_max_subs = df[df["num_mutations"] == max_num_subs]

    # using a set and a list to maintain the order
    variants_set = set()
    variants_list = []

    # create the list of max_subs variants to sample from, basically just shuffle df_max_subs
    df_max_subs = df_max_subs.sample(frac=1, random_state=rng.bit_generator).reset_index(drop=True)
    main_v_index = 0

    while len(variants_list) < target_num:

        # sample the next variant from df_max_subs
        if main_v_index > len(df_max_subs) - 1:
            # there aren't enough df_max_subs variants to sample from to put together the target number of variants
            raise ValueError("Not enough {}-variants to sample from".format(max_num_subs))

        main_v = df_max_subs.iloc[main_v_index]["variant"]
        main_v_index += 1

        # generate all subvariants for this variant
        av = [main_v]
        for i in reversed(range(min_num_subs, max_num_subs)):
            av += get_subvariants(main_v, i)

        # add this variant and subvariants to the main list
        # but first check to ensure they are not already in the list
        # and also check to make sure they ARE in the database (because we are sampling the existing db)
        # note the main_v is guaranteed to be in the database because that's where we sampled it from,
        # but it's just easier to include it with all the subvariants
        for v in av:
            variant_in_set = v in variants_set
            if variant_in_set:
                print("Generated variant already in set: {}".format(v))

            variant_in_db = True
            if v not in db_variants_set:
                variant_in_db = False
                print("Generated subvariant does NOT exist in database, skipping: {}".format(v))

            if not variant_in_set and variant_in_db:
                variants_set.add(v)
                variants_list.append(v)

    return variants_list


def human_format(num):
    """https://stackoverflow.com/questions/579310/formatting-long-numbers-as-strings-in-python"""
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])


def gen_random_main(pdb_fn, seq, seq_idxs, chars, target_num, num_subs_list, num_replicates, seed, out_dir):

    out_fn = "{}_random_TN-{}_NR-{}_NS-{}_RS-{}.txt".format(basename(pdb_fn)[:-4],
                                                            human_format(target_num), num_replicates,
                                                            ",".join(map(str, num_subs_list)), seed)
    out_fn = join(out_dir, out_fn)
    if isfile(out_fn):
        raise FileExistsError("Output file already exists: {}".format(out_fn))
    print("Output file will be {}".format(out_fn))

    # create a random number generator for this call
    rng = np.random.default_rng(seed=seed)

    # generate the variants
    variants = single_pdb_local_variants(seq, target_num, num_subs_list, chars, seq_idxs, rng)
    # multiply number of variants for variance testing
    variants *= num_replicates
    print_variant_info(variants)

    with open(out_fn, "w") as f:
        for v in variants:
            f.write("{} {}\n".format(basename(pdb_fn), v))


def gen_all_main(pdb_fn: str,
                 seq: str,
                 seq_idxs: Sequence[int],
                 chars: list[str],
                 num_subs_list: list[int],
                 out_dir: Optional[str],
                 db_fn: Optional[str] = None,
                 db_mode: Optional[str] = None,
                 db_pdb_fn: Optional[str] = None,
                 ignore_existing_out_file: bool = False):
    """
    Generate all variants for a single PDB file
    db_mode: if 'filter', exclude variants that are in the given database
             if 'sample', only include variants that are in the given database
    """

    if (db_mode is None) ^ (db_fn is None):
        raise ValueError("Both db_fn and db_mode should be specified or left as None")

    # error checking -- db_mode must be "filter" or "sample"
    if db_mode not in [None, "filter", "sample"]:
        raise ValueError("db_mode must be None, 'filter' or 'sample'")

    # if db_pdb_fn is None, set it equal to pdb_fn
    # note db_pdb_fn will only be used if db_mode is "filter" or "sample"
    if db_pdb_fn is None:
        db_pdb_fn = pdb_fn

    # if db_fn is specified, we need to have a hash of the database in the filename
    db_hash = hash_db(db_fn)

    # determine the output filename
    if db_mode == "sample":
        # only sampling variants from the given database
        out_fn = "{}_all_NS-{}_sampled-DB-{}-{}.txt".format(basename(pdb_fn)[:-4],
                                                            ",".join(map(str, num_subs_list)),
                                                            db_hash,
                                                            basename(db_pdb_fn)[:-4])
    elif db_mode == "filter":
        # excluding variants that are in the database
        out_fn = "{}_all_NS-{}_filtered-DB-{}-{}.txt".format(basename(pdb_fn)[:-4],
                                                             ",".join(map(str, num_subs_list)),
                                                             db_hash,
                                                             basename(db_pdb_fn)[:-4])
    else:
        # no database specified, just generate all variants
        out_fn = "{}_all_NS-{}.txt".format(basename(pdb_fn)[:-4], ",".join(map(str, num_subs_list)))

    out_fn = join(out_dir, out_fn)
    # output file already exists
    if isfile(out_fn):
        if ignore_existing_out_file:
            new_out_fn = out_fn
            i = 1
            while isfile(new_out_fn):
                new_out_fn = "{}_{}.txt".format(out_fn[:-4], i)
                i += 1
            out_fn = new_out_fn
        else:
            raise FileExistsError("Output file already exists: {}".format(out_fn))
    print("Output file will be {}".format(out_fn))

    for i in num_subs_list:
        mp = max_possible_variants(len(seq_idxs), i, len(chars))
        print("Generating {} {}-mutation variants".format(mp, i))

    variants = []
    for i in num_subs_list:
        variants += list(gen_all_variants(seq, i, chars, seq_idxs))

    if db_mode == "sample":
        # database sample mode, only include variants that are in the database
        db_variants = load_db_variants(db_fn, db_pdb_fn)
        variants = [v for v in variants if v in db_variants]

    # database filter mode, exclude any variants that are in the database
    elif db_mode == "filter":
        # filter out variants already in the database if one is provided
        db_variants = load_db_variants(db_fn, db_pdb_fn)
        variants = [v for v in variants if v not in db_variants]

    print_variant_info(variants)

    with open(out_fn, "w") as f:
        for v in variants:
            f.write("{} {}\n".format(basename(pdb_fn), v))


def hash_db(db_fn):
    db_hash = 0
    if db_fn is not None:
        print("Hashing database...")
        start = time.time()
        hash_obj = hashlib.shake_128()
        with open(db_fn, "rb") as f:
            for byte_block in iter(lambda: f.read(4096), b""):
                hash_obj.update(byte_block)
        db_hash = hash_obj.hexdigest(4)
        print("Hashing database finished in {}".format(time.time() - start))
    return db_hash


def gen_subvariants_main(pdb_fn: str,
                         seq: str,
                         seq_idxs: Sequence[int],
                         chars: list[str],
                         target_num: int,
                         max_num_subs: int,
                         min_num_subs: int,
                         seed: int,
                         out_dir: str,
                         db_fn: Optional[str] = None,
                         db_mode: Optional[str] = None,
                         db_pdb_fn: Optional[str] = None):

    if (db_mode is None) ^ (db_fn is None):
        raise ValueError("Both db_fn and db_mode should be specified or left as None")

    if db_mode is not None and db_mode not in ["filter", "sample"]:
        raise ValueError("db_mode must be None, 'filter' or 'sample'")

    # db_pdb_fn is used to query the database for database modes 'filter' and 'sample'
    # if None, then use the same PDB file for which we are generating variants
    if db_pdb_fn is None:
        db_pdb_fn = pdb_fn

    # if db_fn is specified, we need to have a hash of the database in the filename
    db_hash = hash_db(db_fn)

    # determine the output filename
    # todo: hard for the filename can't communicate all the provenance...maybe have additional metadata file?
    out_fn_template = "{}_subvariants_TN-{}_MAXS-{}_MINS-{}_{}-DB-{}-{}_RS-{}.txt"
    out_fn_template_args = [
        basename(pdb_fn).rsplit('.', 1)[0],
        human_format(target_num),
        max_num_subs,
        min_num_subs,
        "sampled" if db_mode == "sample" else "filtered",
        db_hash,
        basename(db_pdb_fn).rsplit('.', 1)[0],
        seed
    ]
    out_fn = out_fn_template.format(*out_fn_template_args)

    out_fn = join(out_dir, out_fn)
    if isfile(out_fn):
        raise FileExistsError("Output file already exists: {}".format(out_fn))
    print("Output file will be {}".format(out_fn))

    # create a random number generator for this call
    rng = np.random.default_rng(seed=seed)

    # generate the variants
    if db_mode == "sample":
        # just a type hint because if db_mode is "sample" then the error checking ensures db_fn is str
        db_fn: str
        # sampling needs a special function that selects the main variant from the database
        variants = gen_subvariants_sample(db_fn, db_pdb_fn, target_num, min_num_subs, max_num_subs, rng)
    elif db_mode == "filter" or db_mode is None:
        # this can handle db_fn being None or db_mode being "filter"
        variants = gen_subvariants_vlist(seq, target_num, min_num_subs, max_num_subs, chars, seq_idxs, rng, db_pdb_fn, db_fn)
    else:
        raise ValueError("db_mode must be None, 'filter' or 'sample'")

    print_variant_info(variants)

    # save output to file
    with open(out_fn, "w") as f:
        for v in variants:
            f.write("{} {}\n".format(basename(pdb_fn), v))


def print_variant_info(variants):
    # print out info about the generated variants
    print("Generated {} variants".format(len(variants)))
    count = Counter([len(v.split(",")) for v in variants])
    for k, v in count.items():
        print("{}-mutants: {}".format(k, v))


def get_seq_idxs(seq: str,
                 seq_idxs_range_start: Optional[int],
                 seq_idxs_range_end: Optional[int]) -> np.ndarray:

    range_start = seq_idxs_range_start
    if range_start is None:
        range_start = 0

    range_end = seq_idxs_range_end
    if range_end is None:
        range_end = len(seq)

    seq_idxs = np.arange(range_start, range_end)

    return seq_idxs


def main(args):

    chars = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

    for pdb_fn in args.pdb_fn:
        print("Generating variant list for {}".format(pdb_fn))
        seq = utils.extract_seq_from_pdb(pdb_fn, chain_id=args.chain_id, error_on_multiple_chains=True)
        seq_idxs = get_seq_idxs(seq, args.seq_idxs_range_start, args.seq_idxs_range_end)

        # grab a random, random seed
        seed = args.seed
        if seed is None:
            seed = random.randint(100000000, 999999999)

        if args.method == "subvariants":
            gen_subvariants_main(pdb_fn=pdb_fn,
                                 seq=seq,
                                 seq_idxs=seq_idxs,
                                 chars=chars,
                                 target_num=args.target_num,
                                 max_num_subs=args.max_num_subs,
                                 min_num_subs=args.min_num_subs,
                                 seed=seed,
                                 out_dir=args.out_dir,
                                 db_fn=args.db_fn,
                                 db_mode=args.db_mode,
                                 db_pdb_fn=args.db_pdb_fn)

        elif args.method == "random":
            gen_random_main(pdb_fn, seq, seq_idxs, chars,
                            args.target_num, args.num_subs_list, args.num_replicates, seed, args.out_dir)

        elif args.method == "all":
            gen_all_main(pdb_fn=pdb_fn,
                         seq=seq,
                         seq_idxs=seq_idxs,
                         chars=chars,
                         num_subs_list=args.num_subs_list,
                         out_dir=args.out_dir,
                         db_fn=args.db_fn,
                         db_mode=args.db_mode,
                         db_pdb_fn=args.db_pdb_fn,
                         ignore_existing_out_file=args.ignore_existing_out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars="@")

    parser.add_argument("method",
                        help="what method to use to generate variants",
                        type=str,
                        choices=["all", "random", "subvariants"])

    # common args
    parser.add_argument("--pdb_fn",
                        help="the PDB file from which to generate variants. can specify multiple.",
                        type=str,
                        nargs="+")
    parser.add_argument("--chain_id",
                        help="the chain ID to use for the PDB file",
                        type=str,
                        default=None)
    parser.add_argument("--seq_idxs_range_start",
                        help="the start of the range where to mutate the pdb_fn sequence. 0-based indexing.",
                        type=int,
                        default=None)
    parser.add_argument("--seq_idxs_range_end",
                        help="the end of the range where to mutate the pdb_fn sequence, EXCLUSIVE. 0-based indexing",
                        type=int,
                        default=None)
    parser.add_argument("--target_num",
                        type=int,
                        help="target number of variants per pdb_fn")
    parser.add_argument("--seed",
                        type=int,
                        help="random seed, None for a random random seed",
                        default=None)
    parser.add_argument("--out_dir",
                        type=str,
                        help="output directory for variant lists",
                        default="variant_lists")
    parser.add_argument("--db_fn",
                        type=str,
                        help="database filename, if specified, will not generate variants already in database",
                        default=None)
    parser.add_argument("--db_mode",
                        type=str,
                        help="if 'filter', exclude variants that are in the given database. "
                             "if 'sample', only include variants that are in the given database",
                        default=None,
                        choices=["filter", "sample"])
    parser.add_argument("--db_pdb_fn",
                        type=str,
                        help="the PDB file to use for the database. if None, use the same PDB file as the one "
                             "being used to generate variants",
                        default=None)
    parser.add_argument("--ignore_existing_out_file",
                        action="store_true",
                        default=False,
                        help="ignore existing filename, create a new one with appended number")
    # random args
    parser.add_argument("--num_subs_list",
                        type=int,
                        help="for random and 'all' method, numbers of substitutions for variants",
                        nargs="+",
                        default=[1, 2])
    parser.add_argument("--num_replicates",
                        type=int,
                        help="for random method, the maximum number of replicates of each variant",
                        default=1)

    # subvariants args
    parser.add_argument("--max_num_subs",
                        type=int,
                        help="for subvariants method, the maximum number of substitutions for a variant",
                        default=5)
    parser.add_argument("--min_num_subs",
                        type=int,
                        help="for subvariants method, the minimum number of substitutions for a variant",
                        default=1)

    main(parser.parse_args())
