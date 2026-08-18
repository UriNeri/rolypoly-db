"""Build the external database bundle consumed by RolyPoly.
Not meant to be run directly in entirety.
"""

import datetime
import hashlib
import json
import logging
import os
import shutil
import sqlite3
import subprocess
import tarfile
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
import requests
from bbmapy import bbduk, bbmask, kcompress
from click import Choice
from rich_click import command, option

from rolypoly.utils.bio.alignments import (
    hmmdb_from_directory,
    mmseqs_profile_db_from_directory,
)
from rolypoly.utils.bio.polars_fastx import from_fastx_eager
from rolypoly.utils.bio.sequences import (
    filter_fasta_by_headers,
    remove_duplicates,
    write_fasta_file,
)

from rolypoly.utils.logging.loggit import setup_logging
from rolypoly.utils.various import (
    extract_tar,
    fetch_and_extract,
    read_fwf,
    run_command_comp,
    simple_fetch,
)

ICTV_RANKS = [
    "realm",
    "subrealm",
    "kingdom",
    "subkingdom",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "species",
]

BUILD_STEPS = (
    "genomad",
    "rdrp-scan",
    "rvmt-profiles",
    "rvmt-sequences",
    "neordrp",
    "ncbi-virus-nucleotide",
    "ncbi-virus-taxdb",
    "pfam",
    "vfam",
    "rvmt-motifs",
    "uniref50-virus",
    "contamination",
    "plastid",
    "trna",
    "rfam",
)

# Local debug arguments retained as a copy-paste reference.
data_dir = Path("/home/neri/Documents/GitHub/rps/rolypoly/data")
threads = 6
log_level = "DEBUG"
log_file = "rolypoly_build_data.log"
work_dir = Path("/run/media/neri/ssd2/ncbi_nr/ictv_filter")
selected_steps = set(BUILD_STEPS)
clustered_nr_db = Path(
    "/run/media/neri/ssd2/ncbi_nr/ictv_filter/ncbi_clustered_nr/nr_cluster_seq"
)
clustered_nr_metadata_db = Path(
    "/run/media/neri/ssd2/ncbi_nr/ictv_filter/ncbi_clustered_nr/nr_cluster_seq.sqlite3"
)
clustered_nr_taxonomy_db = Path(
    "/run/media/neri/ssd2/ncbi_nr/ictv_filter/ncbi_clustered_nr/taxonomy4blast.sqlite3"
)
clustered_nr_fetch_dir = Path(
    "/run/media/neri/ssd2/ncbi_nr/work/ncbi_clustered_nr"
)
clustered_nr_work_dir = Path("/run/media/neri/ssd2/ncbi_nr/ictv_filter")
clustered_nr_temp_dir = Path("/run/media/neri/ssd2/ncbi_nr/ictv_filter/temp")
bgzip_threads = 2
force = False
allow_missing_taxonomy = False

# Focused profile refresh:
# selected_steps = {"rfam", "vfam", "pfam"}
#
# Equivalent CLI invocation:
# pixi run build-data -- --data-dir ../rolypoly/data --threads 6 \
#     --step rfam --step vfam --step pfam --log-level DEBUG \
#     --log-file rolypoly_profile_refresh.log


@command()
@option("--data-dir", required=True, help="Path to the data directory")
@option("--threads", default=4, help="Number of threads to use")
@option(
    "--step",
    "steps",
    type=Choice(BUILD_STEPS, case_sensitive=False),
    multiple=True,
    help="Build only this data family; repeat as needed. By default, build all.",
)
@option(
    "--log-file",
    default="./prepare_external_data_logfile.txt",
    help="Path to the log file",
)
@option("--log-level", hidden=True, default="INFO", help="Log level")
@option(
    "--clustered-nr-db",
    envvar="NCBI_CLUSTERED_NR_DB",
    help="Path prefix of the downloaded NCBI nr_cluster_seq BLAST database.",
)
@option(
    "--clustered-nr-metadata-db",
    envvar="NCBI_CLUSTERED_NR_METADATA_DB",
    help=(
        "Path to the ClusteredNR membership SQLite database; defaults to "
        "nr_cluster_seq.sqlite3 beside --clustered-nr-db."
    ),
)
@option(
    "--clustered-nr-taxonomy-db",
    envvar="NCBI_CLUSTERED_NR_TAXONOMY_DB",
    help=(
        "Path to taxonomy4blast.sqlite3 from the NCBI ClusteredNR tarballs; "
        "defaults beside the membership database."
    ),
)
@option(
    "--clustered-nr-fetch-dir",
    envvar="NCBI_CLUSTERED_NR_FETCH_DIR",
    help=(
        "Directory containing pre-fetched nr_cluster_seq tarballs and MD5 files. "
        "By default, they are downloaded under --clustered-nr-work-dir."
    ),
)
@option(
    "--clustered-nr-work-dir",
    envvar="NCBI_CLUSTERED_NR_WORK_DIR",
    help="Working directory for clustered-NR extraction; defaults under --data-dir.",
)
@option(
    "--clustered-nr-temp-dir",
    envvar="NCBI_CLUSTERED_NR_TEMP_DIR",
    help="Temporary directory for MMseqs2 taxonomy construction.",
)
@option(
    "--bgzip-threads",
    default=1,
    type=int,
    envvar="ROLYPOLY_DB_BGZIP_THREADS",
    help="Threads used to bgzip the extracted representative FASTA.",
)
@option("--force", is_flag=True, help="Rebuild completed NCBI taxdb stages.")
@option(
    "--allow-missing-taxonomy",
    is_flag=True,
    help="Assign clustered representatives without an ICTV member mapping to root.",
)
def build_data(
    data_dir,
    threads,
    steps,
    log_file,
    log_level,
    clustered_nr_db,
    clustered_nr_metadata_db,
    clustered_nr_taxonomy_db,
    clustered_nr_fetch_dir,
    clustered_nr_work_dir,
    clustered_nr_temp_dir,
    bgzip_threads,
    force,
    allow_missing_taxonomy,
):
    """Build all or selected external-data families required for RolyPoly.

    With no ``--step``, this runs the complete umbrella workflow. Repeat
    ``--step`` to rebuild one or more selected data families.

    1. Build geNomad RNA viral HMMs
    2. Build protein HMMs RdRp-scan, RVMT, Neordrp_v2.1, tsa_2018 and Pfam 38
    3. Download and prepare rRNA databases SILVA_138.1_SSURef_NR99_tax_silva.fasta and SILVA_138.1_LSURef_NR99_tax_silva.fasta
    4. Download Rfam data.
    5. Download and prepare NCBI ribovirus nucleotide sequences and taxonomy.
    6. Download RVMT sequences and prepare MMseqs2 databases for searches.

    """

    global profile_dir
    global rrna_dir
    global hmmdb_dir
    global mmseqs_dbs
    global contam_dir
    global trna_dir

    logger = setup_logging(log_file, log_level)
    logger.info(f"Starting data preparation to : {data_dir}")

    contam_dir = os.path.join(data_dir, "contam")
    os.makedirs(contam_dir, exist_ok=True)

    rrna_dir = os.path.join(contam_dir, "rrna")
    os.makedirs(rrna_dir, exist_ok=True)

    trna_dir = os.path.join(contam_dir, "trna")
    os.makedirs(trna_dir, exist_ok=True)

    adapter_dir = os.path.join(contam_dir, "adapters")
    os.makedirs(adapter_dir, exist_ok=True)

    masking_dir = os.path.join(contam_dir, "masking")
    os.makedirs(masking_dir, exist_ok=True)

    # taxonomy_dir = os.path.join(data_dir, "taxdump")
    # os.makedirs(taxonomy_dir, exist_ok=True)

    reference_seqs = os.path.join(data_dir, "reference_seqs")
    os.makedirs(reference_seqs, exist_ok=True)

    mmseqs_ref_dir = os.path.join(reference_seqs, "mmseqs")
    os.makedirs(mmseqs_ref_dir, exist_ok=True)

    rvmt_dir = os.path.join(reference_seqs, "RVMT")
    os.makedirs(rvmt_dir, exist_ok=True)

    ncbi_virus_dir = os.path.join(reference_seqs, "ncbi_virus")
    os.makedirs(ncbi_virus_dir, exist_ok=True)

    profile_dir = os.path.join(data_dir, "profiles")
    hmmdb_dir = os.path.join(profile_dir, "hmmdbs")
    mmseqs_dbs = os.path.join(profile_dir, "mmseqs_dbs")

    os.makedirs(hmmdb_dir, exist_ok=True)
    os.makedirs(mmseqs_dbs, exist_ok=True)

    genomad_dir = os.path.join(profile_dir, "genomad")
    os.makedirs(genomad_dir, exist_ok=True)

    selected_steps = set(steps or BUILD_STEPS)

    if "genomad" in selected_steps:
        prepare_genomad_rna_viral_markers(data_dir, threads, logger)
        shutil.rmtree(
            genomad_dir
        )  # outputs have been moved into their final profile directories

    if "rdrp-scan" in selected_steps:
        prepare_rdrp_scan(data_dir, threads, logger)

    if "rvmt-profiles" in selected_steps:
        prepare_RVMT_profiles(data_dir, threads, logger)

    if "rvmt-sequences" in selected_steps:
        prepare_rvmt_mmseqs(data_dir, threads, logger)

    if "neordrp" in selected_steps:
        prepare_neordrp_profiles(data_dir, threads, logger)

    if "ncbi-virus-nucleotide" in selected_steps:
        prepare_ncbi_ribovirus(data_dir, threads, logger)

    if "ncbi-virus-taxdb" in selected_steps:
        prepare_ncbi_virus_taxdb(
            data_dir=data_dir,
            threads=threads,
            logger=logger,
            clustered_nr_db=clustered_nr_db,
            clustered_nr_metadata_db=clustered_nr_metadata_db,
            clustered_nr_taxonomy_db=clustered_nr_taxonomy_db,
            clustered_nr_fetch_dir=clustered_nr_fetch_dir,
            work_dir=clustered_nr_work_dir,
            temp_dir=clustered_nr_temp_dir,
            bgzip_threads=bgzip_threads,
            force=force,
            allow_missing_taxonomy=allow_missing_taxonomy,
        )

    if "pfam" in selected_steps:
        prepare_pfam_rdrps_rt(data_dir, threads, logger)

    if "vfam" in selected_steps:
        prepare_vfam(data_dir, logger)

    if "rvmt-motifs" in selected_steps:
        prepare_rvmt_motifs(data_dir, threads, logger)

    if "uniref50-virus" in selected_steps:
        prepare_uniref50_viral(data_dir, threads, logger)

    if "contamination" in selected_steps:
        prepare_contamination_seqs(data_dir, threads, logger)

    if "plastid" in selected_steps:
        prepare_plastid_data(data_dir, logger)

    if "trna" in selected_steps:
        prepare_trna_data(data_dir, logger)

    if "rfam" in selected_steps:
        download_and_extract_rfam(data_dir, logger)

    # subprocess.run(
    #     "cat NCBI_ribovirus/proteins/datasets_efetch_refseq_ribovirus_proteins_rmdup.faa RVMT/RVMT_allorfs_filtered_no_chimeras.faa | seqkit rmdup | seqkit seq -w0 > prots_for_masking.faa",
    #     shell=True,
    # )
    logger.info("Finished data preparation")


def prepare_rvmt_mmseqs(data_dir, threads, logger: logging.Logger):
    """Prepare RVMT database for seqs searches (mmseqs2 and diamond).

    Processes the RVMT database alignments
    and creates formatted databases for MMseqs2 searches.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors

    Note:
        Downloads RVMT contigs and metadata, filters out chimeric sequences,
        and creates MMseqs2 and compressed databases for fast searches.
    """

    logger.info("Preparing RVMT mmseqs database")

    # Create directories
    rvmt_dir = os.path.join(data_dir, "reference_seqs", "RVMT")
    mmdb_dir = os.path.join(rvmt_dir, "mmseqs")
    os.makedirs(rvmt_dir, exist_ok=True)
    os.makedirs(mmdb_dir, exist_ok=True)

    # Download RVMT contigs
    logger.info("Downloading RVMT contigs")
    contigs_fasta_path = fetch_and_extract(
        "https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/RiboV1.6_Contigs.fasta.gz",
        fetched_to=os.path.join(rvmt_dir, "RiboV1.6_Contigs.fasta.gz"),
        extract_to=rvmt_dir,
        expected_file="RiboV1.6_Contigs.fasta",
        logger=logger,
    )

    # Download and process RVMT info table to get chimeric sequences
    logger.info("Fetching RVMT metadata")
    chimera_ids = []

    info_df = pl.read_csv(
        "https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/RiboV1.6_Info.tsv",
        separator="\t",
        null_values=["NA", ""],
    )

    logger.info("Processing RVMT metadata to identify chimeric sequences")
    # Check for chimeric in `Note` column
    chimera_ids = (
        info_df.filter(
            pl.col("Note")
            .cast(pl.Utf8)
            .str.contains_any(["chim", "rRNA", "cell"], ascii_case_insensitive=True)
        )
        .select(pl.col("ND"))
        .to_series()
        .to_list()
    )
    # Filter for chimeric sequences

    logger.info(f"Found {len(chimera_ids)} chimeric sequences to exclude")

    # Filter out chimeric sequences using rolypoly's filter function
    cleaned_path = os.path.join(rvmt_dir, "RVMT_cleaned_contigs.fasta.gz")
    # logger.info("Filtering out chimeric sequences")
    # filter_fasta_by_headers(
    #     fasta_file=contigs_fasta_path,
    #     headers=chimera_ids,
    #     output_file=cleaned_path,
    #     invert=True,  # Keep sequences NOT in the chimera list
    # )

    ### run cmscan with rrnas to remove more potential chimeras
    logger.info("Running cmscan to identify rRNA sequences in RVMT contigs")
    run_command_comp(
        base_cmd="cmscan",
        positional_args_location="end",
        positional_args=[
            os.path.join(data_dir, "profiles", "cm", "rrna", "rrna.cm"),
            contigs_fasta_path,
        ],
        params={
            "cpu": 8,
            "tblout": "rvmt_rrna.tab",
            "cut_ga": True,
            "noali": True,
        },
        logger=logger,
    )
    rrna_df = read_fwf(
        "rvmt_rrna.tab"
    )  # ,widths=widths, columns=column_names,dtypes = "str")
    rrna_df = rrna_df.filter(pl.col("E-value").cast(pl.Float64) < 1e-3)
    chimera_rrna_ids = rrna_df.select(pl.col("query_name")).to_series().to_list()
    chimera_ids.extend(chimera_rrna_ids)
    logger.info("Filtering out known chimeric and potential rRNA containing sequences")
    filter_fasta_by_headers(
        fasta_file=contigs_fasta_path,
        headers=chimera_ids,
        output_file=cleaned_path,
        invert=True,  # Keep sequences NOT in the chimera list
    )

    # Create MMseqs2 database
    logger.info("Creating MMseqs2 database")
    run_command_comp(
        base_cmd="mmseqs createdb",
        positional_args_location="start",
        positional_args=[cleaned_path, os.path.join(mmdb_dir, "RVMT_cleaned")],
        params={"dbtype": "2"},
        logger=logger,
    )

    # # Create entropy-masked temporary file before compression
    # logger.info("Creating entropy-masked sequences")
    # entropy_masked_path = os.path.join(rvmt_dir, "RVMT_entropy_masked.fasta")
    # from bbmapy import bbmask

    # bbmask(
    #     in1=cleaned_path,
    #     out=entropy_masked_path,
    #     entropy=0.05,
    #     entropywindow=140,
    #     threads=threads,
    # )

    # now similarly, but getting the ORFs
    all_orf_info = pl.read_csv(
        "https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/Simplified_AllORFsInfo.tsv",
        separator="\t",
        null_values=["NA", ""],
    )
    not_chimeric_orfs = (
        all_orf_info.filter(~all_orf_info["seqid"].is_in(chimera_ids))
        .select(pl.col("ORFID"))
        .to_series()
        .to_list()
    )

    rvmt_orfs = fetch_and_extract(
        "https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/RiboV1.5_AllORFs.faa",
        fetched_to=os.path.join(rvmt_dir, "rvmt_orfs.faa"),
        extract_to=rvmt_dir,
        expected_file="rvmt_orfs",
        logger=logger,
    )

    cleaned_orfs_path = os.path.join(rvmt_dir, "RVMT_cleaned_orfs.faa.gz")
    logger.info("Filtering out ORFs from chimeric sequences")
    filter_fasta_by_headers(
        fasta_file=rvmt_orfs,
        headers=not_chimeric_orfs,
        output_file=cleaned_orfs_path,
        invert=False,  # Keep sequences in the non-chimeric ORF list
    )

    # Clean up temporary files, only keep the compressed
    try:
        os.remove(rvmt_dir + "/RiboV1.6_Contigs.fasta")
        os.remove(rvmt_dir + "/RiboV1.6_Contigs.fasta.gz")
        os.remove(rvmt_dir + "/rvmt_orfs.faa")
    except FileNotFoundError:
        logger.warning(
            "some temporary files for RVMT mmseqs preparation might not have been cleaned."
        )

    logger.info(f"RVMT databases created successfully in {rvmt_dir} and {mmdb_dir}")


def download_and_extract_rfam(data_dir, logger):
    """Download and process Rfam database files.

    Retrieves Rfam database files and processes them for use in RNA (structural)
    family identification and annotation.

    Args:
        data_dir (str): Base directory for data storage
        logger: Logger object for recording progress and errors

    Note:
        Downloads both the sequence database and covariance models,
        and processes them for use with Infernal.
    """

    rfam_url = "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
    rfam_cm_path = data_dir + "/Rfam.cm.gz"
    rfam_extract_path = data_dir + "/profiles/cm/"
    # subprocess.run("cmpress Rfam.cm", shell=True)

    logger.info("Downloading Rfam database")
    try:
        fetch_and_extract(rfam_url, extract_to=rfam_extract_path)
        logger.info("Rfam database downloaded and extracted successfully.")
    except requests.exceptions.RequestException as e:
        logger.error(f"Error downloading Rfam database: {e}")

    # press the cm file
    run_command_comp(
        base_cmd="cmpress",
        positional_args_location="start",
        positional_args=[os.path.join(rfam_extract_path, "Rfam.cm")],
        logger=logger,
    )


def prepare_vfam(data_dir, logger: logging.Logger):
    """Build filtered VFam HMM and MMseqs2 profile databases."""
    vfam_url = "https://fileshare.lisc.univie.ac.at/vog/latest"
    metadata_path = os.path.join(data_dir, "profiles", "vfam.annotations.tsv.gz")
    work_dir = os.path.join(hmmdb_dir, "vfam")
    archive_path = os.path.join(work_dir, "vfam.raw_algs.tar.gz")
    msa_dir = os.path.join(work_dir, "msa")
    os.makedirs(work_dir, exist_ok=True)

    logger.info("Downloading VFam annotations and alignments")
    fetch_and_extract(
        url=f"{vfam_url}/vfam.annotations.tsv.gz",
        fetched_to=metadata_path,
        extract=False,
        logger=logger,
    )
    fetch_and_extract(
        url=f"{vfam_url}/vfam.raw_algs.tar.gz",
        fetched_to=archive_path,
        extract=False,
        logger=logger,
    )
    extract_tar(Path(archive_path), Path(work_dir), logger)

    vfam_df = pl.read_csv(metadata_path, separator="\t")
    valid_accessions = set(
        vfam_df.filter(
            ~pl.col("ConsensusFunctionalDescription")
            .str.to_lowercase()
            .str.contains("hypothetical")
            & (pl.col("SpeciesCount") > 2)
        )["#GroupName"].to_list()
    )
    removed = 0
    for msa_path in Path(msa_dir).glob("*.msa"):
        if msa_path.stem not in valid_accessions:
            msa_path.unlink()
            removed += 1
    logger.info(
        f"Retained {len(valid_accessions):,} informative VFam accessions; "
        f"removed {removed:,} alignments"
    )

    hmmdb_from_directory(
        msa_dir=msa_dir,
        output=os.path.join(hmmdb_dir, "vfam.hmm"),
        msa_pattern="*.msa",
        info_table=metadata_path,
        name_col="#GroupName",
        accs_col="#GroupName",
        desc_col="ConsensusFunctionalDescription",
        gath_col=None,
    )
    mmseqs_profile_db_from_directory(
        msa_dir=msa_dir,
        output=os.path.join(mmseqs_dbs, "vfam", "vfam"),
        msa_pattern="*.msa",
        info_table=metadata_path,
        name_col="#GroupName",
        accs_col="#GroupName",
        desc_col="ConsensusFunctionalDescription",
    )

    shutil.rmtree(work_dir)
    logger.info("Created filtered VFam HMM and MMseqs2 databases")


def prepare_uniref50_viral(data_dir, threads, logger):
    """Download and prepare UniRef50 viral protein sequences."""
    logger.info("Downloading UniRef50 viral protein sequences")
    os.makedirs(os.path.join(data_dir, "reference_seqs/uniref"), exist_ok=True)
    uniref50_viral_fasta = os.path.join(
        data_dir, "reference_seqs/uniref/uniref50_viral.fasta"
    )
    uniref50_viral_fasta_gz = uniref50_viral_fasta + ".gz"
    fetch_and_extract(
        "https://rest.uniprot.org/uniref/stream?compressed=true&format=fasta&query=%28%28identity%3A0.5%29+AND+%28taxonomy_id%3A10239%29+AND+%28count%3A%5B1+TO+192133%5D%29%29",
        fetched_to=uniref50_viral_fasta_gz,
        extract_to=os.path.dirname(uniref50_viral_fasta),
        expected_file=os.path.basename(uniref50_viral_fasta),
        logger=logger,
        debug=True,
    )
    uniref50_viral_data = os.path.join(
        data_dir, "reference_seqs/uniref/uniref50_viral.tsv"
    )
    uniref50_viral_data_gz = uniref50_viral_data + ".gz"
    fetch_and_extract(
        "https://rest.uniprot.org/uniref/stream?compressed=true&fields=id%2Cname%2Ctypes%2Ccount%2Corganism%2Clength%2Cidentity%2Cmembers&format=tsv&query=%28%28identity%3A0.5%29+AND+%28taxonomy_id%3A10239%29+AND+%28count%3A%5B1+TO+192133%5D%29%29",
        fetched_to=uniref50_viral_data_gz,
        extract_to=os.path.dirname(uniref50_viral_fasta),
        extract=True,
        logger=logger,
        debug=True,
    )
    # why the hell is it still compressed??? 2 times gzipped???
    from rolypoly.utils.various import extract

    extract(uniref50_viral_data_gz, uniref50_viral_data)

    # remove sequences of uncharacterized/hypothetical proteins, or non-informative such as "polyprotein fragment"
    logger.info("Filtering UniRef50 viral protein sequences")
    uniref_df = pl.read_csv(uniref50_viral_data, separator="\t", null_values=["NA", ""])
    to_remove_terms = [
        "uncharacterized protein",
        "hypothetical protein",
        "putative protein",
        "predicted protein",
        "unnamed protein product",
        "polyprotein fragment",
        "genome polyprotein",
        "fragment",
    ]
    pattern = "|".join(to_remove_terms)
    filtered_uniref_df = uniref_df.filter(
        ~pl.col("Cluster Name").str.to_lowercase().str.contains(pattern)
    )
    filtered_ids = filtered_uniref_df.select(pl.col("Cluster ID")).to_series().to_list()
    logger.info(
        f"Removing {uniref_df.height - filtered_uniref_df.height} non-informative sequences from UniRef50 viral dataset"
    )
    # filter fasta
    # first load the fasta to a dataframe to get the original headers containing the Cluster IDs + defline and remove based on only the IDs
    unidf = from_fastx_eager(uniref50_viral_fasta).with_columns(
        pl.col("header").str.split(" ").list.first().alias("Cluster ID")
    )
    unidf_filtered = unidf.filter(pl.col("Cluster ID").is_in(filtered_ids))
    write_fasta_file(
        format="fasta",
        headers=unidf_filtered["header"].to_list(),
        seqs=unidf_filtered["sequence"].to_list(),
        output_file=uniref50_viral_fasta_gz,
    )
    uniref_df.write_csv(
        uniref50_viral_data_gz, separator="\t", compression="gzip"
    )

    # clean up temporary files, keeping the final compressed metadata table
    try:
        os.remove(uniref50_viral_fasta)
        os.remove(uniref50_viral_data)
    except FileNotFoundError:
        logger.warning(
            "some temporary files for UniRef50 viral preparation might not have been cleaned."
        )

    # https://rest.uniprot.org/uniref/stream?compressed=true&fields=id%2Cname%2Ctypes%2Ccount%2Corganism%2Clength%2Cidentity%2Cmembers&format=tsv&query=%28%28identity%3A0.5%29+AND+%28taxonomy_id%3A10239%29+AND+%28count%3A%5B1+TO+192133%5D%29%29

    # https://rest.uniprot.org/uniref/stream?compressed=true&format=fasta&query=%28%28identity%3A0.5%29+AND+%28taxonomy_id%3A10239%29+AND+%28count%3A%5B1+TO+192133%5D%29%29


def prepare_genomad_rna_viral_markers(data_dir, threads, logger: logging.Logger):
    """Download and prepare RNA viral HMMs from geNomad markers.

    Downloads the geNomad database, analyzes the marker metadata to identify
    RNA viral specific markers, and creates an HMM database from their alignments.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors
    """

    logger.info("Starting geNomad RNA viral HMM preparation")

    # Create directories
    genomad_dir = os.path.join(data_dir, "profiles/genomad")
    genomad_db_dir = os.path.join(genomad_dir, "genomad_db")
    genomad_markers_dir = os.path.join(genomad_db_dir, "markers")
    genomad_alignments_dir = os.path.join(genomad_markers_dir, "alignments")
    os.makedirs(genomad_dir, exist_ok=True)
    os.makedirs(genomad_db_dir, exist_ok=True)
    os.makedirs(genomad_markers_dir, exist_ok=True)
    os.makedirs(genomad_alignments_dir, exist_ok=True)
    genomad_info_table = os.path.join(
        data_dir, "profiles/genomad_rna_viral_markers_with_annotation.csv.gz"
    )

    # Download metadata and database
    genomad_data = "https://zenodo.org/api/records/14886553/files-archive"  # noqa
    db_url = (
        "https://zenodo.org/records/14886553/files/genomad_msa_v1.9.tar.gz?download=1"
    )
    metadata_url = "https://zenodo.org/records/14886553/files/genomad_metadata_v1.9.tsv.gz?download=1"
    # Download and read metadata
    logger.info("Downloading geNomad metadata")
    aria2c_command = (
        f"aria2c -c -d {genomad_dir} -o ./genomad_metadata_v1.9.tsv.gz {metadata_url}"
    )
    subprocess.run(aria2c_command, shell=True)

    metadata_df = pl.read_csv(
        f"{genomad_dir}/genomad_metadata_v1.9.tsv.gz",
        separator="\t",
        null_values=["NA"],
        infer_schema_length=10000,
    )
    # only virus specific markers
    # metadata_df = metadata_df.filter(pl.col("SPECIFICITY_CLASS") == "VV")
    # only RNA viral markers
    metadata_df = metadata_df.filter(
        pl.col("ANNOTATION_DESCRIPTION")
        .str.to_lowercase()
        .str.contains("rna-dependent rna polymerase")
        | pl.col("TAXONOMY").str.contains("Riboviria")
        | pl.col("SOURCE").str.contains("RVMT")
    )
    # Next, filling missing annotation from InterPro.
    # only ones without description
    to_fill = metadata_df.filter(pl.col("ANNOTATION_DESCRIPTION").is_null())
    # if multiple maybe split->explode->groupby->agg->majority vote, something like:
    # nah just using the first if multiple maybe split->explode->first:
    # for now, only using a single accession (but word cloud/majority vote would probably work for multiple accessions)
    # to_fill = to_fill.with_columns(pl.col("ANNOTATION_ACCESSIONS").str.split(";").list.first().alias("ANNOTATION_ACCESSIONS_first"))
    to_fill = to_fill.with_columns(
        pl.col("ANNOTATION_ACCESSIONS")
        .str.split(";")
        .alias("ANNOTATION_ACCESSIONS_struct")
    )
    to_fill = to_fill.explode("ANNOTATION_ACCESSIONS_struct")

    to_fill = to_fill.with_columns(
        pl.when(pl.col("ANNOTATION_ACCESSIONS").str.starts_with("PF"))
        .then(pl.lit("Pfam"))
        .when(pl.col("ANNOTATION_ACCESSIONS").str.starts_with("COG"))
        .then(pl.lit("COG"))
        .when(pl.col("ANNOTATION_ACCESSIONS").str.starts_with("K"))
        .then(pl.lit("KEGG"))
        .otherwise(pl.lit("unknown"))
        .alias("source_db")
    )

    # We (currently) only carte about viral specific markers, so filtering out the rest
    # Not sure Kegg is on interpro.
    to_fill = to_fill.filter(pl.col("source_db").str.contains("Pfam|COG"))

    def query_interpro(entry: str, source_db: str):
        """Fetch the InterPro description for a given entry."""
        # from bs4 import BeautifulSoup
        import requests

        if source_db == "unknown":
            return None
        url = f"https://www.ebi.ac.uk/interpro/api/entry/{source_db}/{entry}"
        # print(url)                     #  debugging

        response = requests.get(url)
        if response.status_code != 200:
            return None

        data = response.json()
        # print(data)                     #  debugging

        # desc = data.get("metadata", {}).get("description")
        desc = data.get("metadata", {}).get("name", {}).get("name", None)
        return desc

    filled_interpro = []
    from rich.progress import track

    tiny_fill = to_fill.select(["ANNOTATION_ACCESSIONS_struct", "source_db"]).unique()

    for row in track(tiny_fill.to_dicts()):
        if row["source_db"] == "unknown":
            filled_interpro.append(None)
        else:
            this_desc = query_interpro(
                row["ANNOTATION_ACCESSIONS_struct"], row["source_db"]
            )
            filled_interpro.append(this_desc)
            print(f"{row['ANNOTATION_ACCESSIONS_struct']}\t{this_desc}")  # debugging
            # filled_interpro.append(query_interpro(row["ANNOTATION_ACCESSIONS"], row["source_db"]))

    tiny_fill = tiny_fill.with_columns(pl.Series(filled_interpro).alias("interpro"))
    to_fill = to_fill.join(
        tiny_fill, on=["ANNOTATION_ACCESSIONS_struct", "source_db"], how="left"
    )
    to_fill = to_fill.with_columns(
        pl.coalesce(pl.col("ANNOTATION_DESCRIPTION"), pl.col("interpro")).alias(
            "ANNOTATION_DESCRIPTION"
        )
    )
    to_fill = to_fill.drop(
        "ANNOTATION_ACCESSIONS_struct", "interpro", "source_db"
    ).unique()
    to_fill = to_fill.filter(pl.col("ANNOTATION_DESCRIPTION").is_not_null())
    # hopefully now we have filled some of the missing descriptions, and we don't have any duplicate MARKERs

    metadata_df = metadata_df.filter(
        ~pl.col("MARKER").is_in(to_fill["MARKER"].implode())
    )
    metadata_df = metadata_df.vstack(to_fill)

    metadata_df.write_csv(genomad_info_table, compression="gzip")

    # Download MSAs
    logger.info("Downloading geNomad database")
    aria2c_command = f"aria2c -c -d {genomad_dir} -o ./genomad_msa_v1.9.tar.gz {db_url}"
    subprocess.run(aria2c_command, shell=True)

    # Extract RNA viral MSAs
    marker_ids = metadata_df["MARKER"].to_list()

    with tarfile.open(f"{genomad_dir}/genomad_msa_v1.9.tar.gz", "r") as tar:
        for member in tar.getmembers():
            if (
                member.name.removeprefix("genomad_msa_v1.9/").removesuffix(".faa")
                in marker_ids
            ):
                tar.extract(member, genomad_alignments_dir)
    # need to move all files in genomad/genomad_db/markers/alignments/genomad_msa_v1.9/* to genomad/genomad_db/markers/alignments/
    for file in os.listdir(genomad_alignments_dir + "/genomad_msa_v1.9"):
        shutil.move(
            genomad_alignments_dir + "/genomad_msa_v1.9/" + file,
            genomad_alignments_dir + "/" + file,
        )
    # remove the genomad_msa_v1.9 directory
    shutil.rmtree(genomad_alignments_dir + "/genomad_msa_v1.9")

    output_hmm = os.path.join(
        os.path.join(data_dir, "profiles/hmmdbs"),
        "genomad_rna_viral_markers.hmm",
    )
    hmmdb_from_directory(
        msa_dir=genomad_alignments_dir,
        output=output_hmm,
        msa_pattern="*.faa",
        info_table=genomad_info_table,
        name_col="MARKER",
        accs_col="ANNOTATION_ACCESSIONS",
        desc_col="ANNOTATION_DESCRIPTION",
        gath_col=None,  # no gathering theshold pre-defined for genomad
    )

    mmseqs_profile_db_from_directory(
        msa_dir=genomad_alignments_dir,
        output=os.path.join(
            data_dir, "profiles/mmseqs_dbs/genomad/", "rna_viral_markers"
        ),
        info_table=genomad_info_table,
        msa_pattern="*.faa",
        name_col="MARKER",
        accs_col="ANNOTATION_ACCESSIONS",
        desc_col="ANNOTATION_DESCRIPTION",
    )
    # clean up
    try:
        os.remove(f"{genomad_dir}/genomad_metadata_v1.9.tsv.gz")
        os.remove(f"{genomad_dir}/genomad_msa_v1.9.tar.gz")
        shutil.rmtree(genomad_db_dir)
    except Exception as e:
        logger.warning(f"Could not remove file: {e}")

    logger.info(f"Created genomad RNA viral HMM database at {output_hmm}")


def prepare_rdrp_scan(data_dir, threads, logger: logging.Logger):
    """Download and prepare RdRp profiles from RdRp-scan.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors
    """

    logger.info("Preparing RdRp-scan HMM and MMseqs databases")
    fetch_and_extract(
        "https://github.com/JustineCharon/RdRp-scan/archive/refs/heads/main.zip",
        fetched_to=hmmdb_dir + "/RdRp-scan.zip",
        extract_to=hmmdb_dir + "/RdRp-scan",
    )

    # Use utility function to build HMM database from MSAs
    rdrp_scan_msa_dir = os.path.join(
        hmmdb_dir, "RdRp-scan/RdRp-scan-main/Profile_db_and_alignments"
    )
    rdrp_scan_output = os.path.join(hmmdb_dir, "rdrp_scan.hmm")
    hmmdb_from_directory(
        msa_dir=rdrp_scan_msa_dir,
        output=rdrp_scan_output,
        msa_pattern="*.fasta.CLUSTALO",
    )
    # Also build an MMseqs profile DB from the RdRp-scan MSAs for fast searches
    mmseqs_profile_db_from_directory(
        msa_dir=rdrp_scan_msa_dir,
        output=os.path.join(mmseqs_dbs, "rdrp_scan/rdrp_scan"),
        msa_pattern="*.fasta.CLUSTALO",
        info_table=None,
    )
    # clean up
    shutil.rmtree(hmmdb_dir + "/RdRp-scan")
    os.remove(hmmdb_dir + "/RdRp-scan.zip")

    logger.debug("Finished preparing rdrp-scan databases")


def prepare_pfam_rdrps_rt(data_dir, threads, logger: logging.Logger):
    """Download and prepare RdRp profiles from PFAM.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors
    """

    logger.info("Preparing PFAM RdRps and RTs HMM database")

    # Pfam 38.2 RdRps and RTs
    # fetch Pfam-A.hmm.gz to the hmmdb directory and extract into that directory
    pfam_hmm_url = (
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.hmm.gz"
    )
    pfam_msa_url = (
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/Pfam-A.fasta.gz"
    )
    pfam_gz_path = os.path.join(hmmdb_dir, "Pfam-A.hmm.gz")
    os.makedirs(hmmdb_dir, exist_ok=True)
    fetch_and_extract(url=pfam_hmm_url, fetched_to=pfam_gz_path, extract_to=hmmdb_dir)

    RdRps_and_RTs = [
        "PF04197.17",
        "PF04196.17",
        "PF22212.1",
        "PF22152.1",
        "PF22260.1",
        "PF00680.25",
        "PF00978.26",
        "PF00998.28",
        "PF02123.21",
        "PF07925.16",
        "PF00078.32",
        "PF07727.19",
        "PF13456.11",
    ]

    # Use hmmfetch to extract the small set of Pfam HMMs we care about
    selected_pfam_output = os.path.join(hmmdb_dir, "pfam_rdrps_and_rts.hmm")
    from rolypoly.utils.bio.alignments import hmm_fetch

    hmm_fetch(
        accessions=RdRps_and_RTs,
        hmm_db=os.path.join(hmmdb_dir, "Pfam-A.hmm"),
        output=selected_pfam_output,
        strip_after_char=".",
        logger=logger,
    )

    # also prepare mmseqs profile db
    subprocess.run(
        "mmseqs databases Pfam-A.seed data/profiles/mmseqs_dbs/pfam_a/pfam_a_38_seed tmp",
        shell=True,
        check=True,
    )
    # TODO: Filter the pfam_A mmseqs db to remove "hypothetical protein" + similar entries?

    from rolypoly.utils.bio.alignments import fetchPfamMSA

    # fetch Pfam MSAs for these accessions
    pfam_msa_folder = os.path.join(hmmdb_dir, "pfam_msas")
    os.makedirs(pfam_msa_folder, exist_ok=True)
    for acc in [acc.split(".")[0] for acc in RdRps_and_RTs]:
        fetchPfamMSA(acc=acc, output_folder=pfam_msa_folder, logger=logger)
        logger.info(f"Fetched MSA for Pfam accession {acc}")
    # build mmseqs profile db from these msas
    mmseqs_profile_db_from_directory(
        msa_dir=pfam_msa_folder,
        output=os.path.join(mmseqs_dbs, "pfam_rdrps_and_rts/pfam_rdrps_and_rts"),
        msa_pattern="*.sth",
        info_table=None,
    )

    # clean up downloaded gz
    try:
        os.remove(pfam_gz_path)
        os.remove("tmp")

    except Exception:
        pass

    logger.debug("Finished preparing pfam and pfam rdrps+rts sub-database")


def prepare_neordrp_profiles(data_dir, threads, logger: logging.Logger):
    """Prepare NeoRdRp v2.1 RdRp profile HMM database
    NOTE: this is a USING A PRECOMPUTED HMM, not building from MSA! MIGHT NOT BE COMPATIBLE WITH FUTURE VERSIONS OF HMMER!
    """
    logger.info("Preparing NeoRdRp v2.1 HMM database")
    neordrp_url = (
        "https://zenodo.org/records/10851672/files/NeoRdRp.2.1.hmm.xz?download=1"
    )
    neordrp_path = os.path.join(hmmdb_dir, "NeoRdRp.2.1.hmm.xz")
    fetch_and_extract(url=neordrp_url, fetched_to=neordrp_path, extract_to=hmmdb_dir)
    shutil.move(
        os.path.join(hmmdb_dir, "NeoRdRp.2.1.hmm"),
        os.path.join(hmmdb_dir, "neordrp2.1.hmm"),
    )
    os.unlink(neordrp_path)


def prepare_RVMT_profiles(data_dir, threads, logger: logging.Logger):
    """Prepare RVMT RdRp profile + annotation db (NVPC)."""
    logger.info("Preparing RVMT HMM and MMseqs databases")
    rvmt_url = "https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/Alignments/zip.ali.220515.tgz"
    rvmt_path = os.path.join(hmmdb_dir, "zip.ali.220515.tgz")
    fetch_and_extract(
        url=rvmt_url,
        fetched_to=rvmt_path,
        extract_to=os.path.join(hmmdb_dir, "RVMT/"),
    )
    hmmdb_from_directory(
        msa_dir=os.path.join(hmmdb_dir, "RVMT/"),
        output=os.path.join(hmmdb_dir, "rvmt.hmm"),
        msa_pattern="ali*/*.FASTA",
        info_table=None,
        default_description="RdRp"
    )

    mmseqs_profile_db_from_directory(
        msa_dir=os.path.join(hmmdb_dir, "RVMT/"),
        output=os.path.join(mmseqs_dbs, "RVMT/RVMT"),
        msa_pattern="ali*/*.FASTA",
        info_table=None,
    )

    # now for nvpc
    fetch_and_extract(
        url="https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/zenodo/v4/RVMT_Zenodo_V4/Domains_Annotations/NVPC/msaFiles.tar.gz",
        fetched_to=os.path.join(hmmdb_dir, "NVPC_msaFiles.tar.gz"),
        extract_to=os.path.join(hmmdb_dir, "RVMT/NVPC/"),
    )

    # nvpc_info = pl.read_csv("https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/zenodo/v4/RVMT_Zenodo_V4/Domains_Annotations/NVPC/NVPC_info.tsv",
    #     separator="\t",
    #     null_values=["NA", ""],
    # ).rename({"New_Name":"Name"}).unique()
    # nvpc_descriptions = pl.read_excel("https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/zenodo/v4/RVMT_Zenodo_V4/Domains_Annotations/misc/NeoCM3_full.xlsx")
    # nvpc_descriptions = nvpc_descriptions.select(["profile_accession","New_Name","Comment"]).rename({"New_Name":"Name","Comment":"Description"}).unique()
    # # add a column with now manyu times each new_name appears
    # # nvpc_descriptions = nvpc_descriptions.filter(pl.col("Description").str.contains_any(["rdrp_fragment","caution","the profile should be split or ditched altogeth"],ascii_case_insensitive=True))
    # nvpc_info = nvpc_info.join(temp_df, on="Name", how="inner")
    # nvpc_descriptions = nvpc_descriptions.join(nvpc_info, on="Name", how="inner")
    # nvpc_descriptions = nvpc_descriptions.filter(pl.col("profile_accession").is_in(nvpc_info["profile_accession"].implode()))
    # #  sort by number of new_name (most duplicates first)
    # nvpc_descriptions.write_csv(os.path.join(hmmdb_dir, "RVMT/NVPC/NVPC_descriptions.csv"),include_header=True)

    nvpc_info_table = os.path.join(profile_dir, "NVPC_descriptions.csv.gz")
    info_table = pl.read_csv(nvpc_info_table)
    # remove any msa file that doesn't have a matching profile_accession in the info table
    import glob

    msa_files = glob.glob(os.path.join(hmmdb_dir, "RVMT/NVPC/msaFiles/*.afa"))
    valid_accessions = set(info_table["profile_accession"].to_list())
    for msa_file in msa_files:
        base_name = os.path.basename(msa_file)
        # get the NVPC.NNNn part (first two parts separated by .)
        acc = ".".join(os.path.splitext(base_name)[0].split(".")[:2])
        if acc not in valid_accessions:
            logger.info(f"Removed MSA file without matching accession: {msa_file}")
            os.remove(msa_file)

    hmmdb_from_directory(
        msa_dir=os.path.join(hmmdb_dir, "RVMT/NVPC/"),
        output=os.path.join(hmmdb_dir, "nvpc.hmm"),
        msa_pattern="msaFiles/*.afa",
        info_table=nvpc_info_table,
        accs_col="profile_accession",
        name_col="Name",
        desc_col="Description",
        missing_include=False,
        default_gath="5",
        debug=False,
    )

    mmseqs_profile_db_from_directory(
        msa_dir=os.path.join(hmmdb_dir, "RVMT/NVPC/msaFiles/"),
        output=os.path.join(mmseqs_dbs, "nvpc/nvpc"),
        msa_pattern="*.afa",
        info_table=nvpc_info_table,
        accs_col="profile_accession",
        name_col="Name",
        desc_col="Description",
        # missing_include=False,
    )

    # clean up
    shutil.rmtree(os.path.join(hmmdb_dir, "RVMT/"))
    os.remove(os.path.join(hmmdb_dir, "zip.ali.220515.tgz"))
    os.remove(os.path.join(hmmdb_dir, "NVPC_msaFiles.tar.gz"))
    logger.info("Finished preparing RVMT databases")


def prepare_ncbi_ribovirus(data_dir, threads, logger: logging.Logger):
    """Download and prepare NCBI ribovirus reference sequences (RefSeq only).

    Downloads complete RefSeq genomes for RNA viruses (Riboviria), processes them
    with entropy masking and compression for efficient searches.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors
    """

    logger.info("Preparing NCBI ribovirus reference sequences")

    ncbi_virus_dir = os.path.join(data_dir, "reference_seqs", "ncbi_virus")
    os.makedirs(ncbi_virus_dir, exist_ok=True)
    mmdb_dir = os.path.join(ncbi_virus_dir, "mmseqs")
    os.makedirs(mmdb_dir, exist_ok=True)

    # Define file paths
    raw_fasta_path = os.path.join(ncbi_virus_dir, "refseq_ribovirus_genomes.fasta")
    entropy_masked_path = os.path.join(
        ncbi_virus_dir, "refseq_ribovirus_genomes_entropy_masked.fasta"
    )  # noqa (F841)
    compressed_path = os.path.join(
        ncbi_virus_dir, "refseq_ribovirus_genomes_flat.fasta"
    )  # noqa (F841)

    # Riboviria taxid
    taxid = "2559587"

    # Use esearch and efetch to download complete RefSeq ribovirus genomes
    logger.info(f"Downloading RefSeq ribovirus genomes for taxid {taxid}")

    # if from_ena == True:
    #     # Use EBI/ENA REST API instead of NCBI E-utilities
    #     logger.info("Searching EBI/ENA for Riboviria sequences")

    #     # EBI/ENA API search for Riboviria complete genomes
    #     import requests

    #     # Search for Riboviria sequences in ENA
    #     search_url = "https://www.ebi.ac.uk/ena/portal/api/search"
    #     search_params = {
    #         "result": "sequence",
    #         "query": f'tax_tree({taxid}) AND mol_type="genomic RNA" AND base_count>1000',
    #         "fields": "accession,scientific_name,description,mol_type,tax_id,tax_lineage",
    #         "format": "json",
    #         "limit": "100"  # Get all results
    #     }

    #     logger.info("Querying EBI/ENA for sequence metadata")
    #     response = requests.get(search_url, params=search_params)
    #     response.raise_for_status()

    #     sequences_metadata = response.json()
    #     logger.info(f"Found {len(sequences_metadata)} Riboviria sequences")

    #     # Get FASTA sequences using EBI API
    #     logger.info("Downloading sequences from EBI/ENA")
    #     accessions = [seq["accession"] for seq in sequences_metadata[:1000]]  # Limit to avoid overwhelming

    #     with open(raw_fasta_path, "w") as fasta_out:
    #         for i, accession in enumerate(accessions):
    #             if i % 50 == 0:
    #                 logger.info(f"Downloaded {i}/{len(accessions)} sequences")

    #             # Get FASTA from EBI
    #             fasta_url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{accession}"
    #             fasta_response = requests.get(fasta_url)

    #             if fasta_response.status_code == 200:
    #                 fasta_out.write(fasta_response.text)
    #                 fasta_out.write("\n")
    #             else:
    #                 logger.warning(f"Failed to download {accession}: {fasta_response.status_code}")

    #     logger.info("Downloaded RefSeq ribovirus genomes from EBI/ENA")

    #     # Apply entropy masking first (consistent with RVMT approach)
    #     logger.info("Applying entropy masking")

    # if from_edirect == True:
    # esearch_query = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP] AND complete genome[title]"
    #     logger.info("Running esearch | efetch pipeline")
    #     pipeline_cmd = f"~/bin/edirect/esearch -db nuccore -query '{esearch_query}' | ~/bin/edirect/efetch -format fasta > {raw_fasta_path}"

    #     run_command_comp(
    #         base_cmd=pipeline_cmd,
    #         params={},
    #         output_file=raw_fasta_path,
    #         logger=logger,
    #         check_output=True
    #     )

    from_ncbi_ftp = True  # for now, above methods is 1. edirect dependent, 2. ENA API dependent which seems slow/limited (or I'm not filtering prorperly - very likely)
    if from_ncbi_ftp == True:
        # genomes
        fetch_and_extract(
            url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz",
            fetched_to=os.path.join(ncbi_virus_dir, "viral.1.1.genomic.fna.gz"),
            extract_to=ncbi_virus_dir,
            rename_extracted=raw_fasta_path,
        )
        # orfs
        fetch_and_extract(
            url="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz",
            fetched_to=os.path.join(ncbi_virus_dir, "viral.1.protein.faa.gz"),
            extract_to=ncbi_virus_dir,
            rename_extracted=raw_fasta_path.replace(".fasta", "_orfs.faa"),
        )

    logger.info("Downloaded NCBI ribovirus genomes")

    # Create MMseqs2 database
    logger.info("Creating MMseqs2 database")
    run_command_comp(
        base_cmd="mmseqs createdb",
        positional_args_location="start",
        positional_args=[
            raw_fasta_path,
            os.path.join(mmdb_dir, "ncbi_ribovirus_cleaned"),
        ],
        params={"dbtype": "2"},
        logger=logger,
    )

    # Clean up intermediate files
    try:
        os.remove(os.path.join(ncbi_virus_dir, "viral.1.1.genomic.fna.gz"))
        os.remove(os.path.join(ncbi_virus_dir, "viral.1.protein.faa.gz"))
    except FileNotFoundError:
        logger.warning("Some intermediate files not found for cleanup")

    logger.info(f"NCBI ribovirus preparation completed in {ncbi_virus_dir}")


def prepare_rvmt_motifs(data_dir, threads, logger):
    """Prepare RVMT motif sequences for profile-based searches.

    Extracts and processes the RVMT motif sequence library from the pre-downloaded
    tar.gz file, organizing motifs by type (A=mot.1, B=mot.2, C=mot.3) and taxon.
    Creates both HMM and MMseqs profile databases for fast searches.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors

    Note:
        This function assumes motif_sequence_library.tar.gz has been downloaded
        to data_dir/profiles/motif_sequence_library.tar.gz
    """

    logger.info("Preparing RVMT motif sequences")
    motif_dir = os.path.join(data_dir, "profiles/rvmt_motifs")
    motif_alignments_dir = os.path.join(motif_dir, "Sequence_Library")

    os.makedirs(motif_dir, exist_ok=True)
    os.makedirs(motif_alignments_dir, exist_ok=True)

    motif_archive = os.path.join(data_dir, "profiles/motif_sequence_library.tar.gz")
    # fetch motif archive
    logger.info(
        "fetching Motif archive from https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/rdrps/motif_sequence_library.tar.gz"
    )
    fetch_and_extract(
        url="https://portal.nersc.gov/dna/microbial/prokpubs/Riboviria/RiboV1.4/rdrps/motif_sequence_library.tar.gz",
        fetched_to=motif_archive,
        extract_to=motif_dir,
        logger=logger,
    )

    # The archive contains Sequence_Library/ with mot.1/, mot.2/, mot.3/ subdirectories
    sequence_library_dir = os.path.join(motif_dir, "Sequence_Library")

    if not os.path.exists(sequence_library_dir):
        logger.error("Expected Sequence_Library directory not found after extraction")
        return False

    # Process each motif type directory
    motif_metadata = {}

    for motif_type_dir in os.listdir(sequence_library_dir):
        if not motif_type_dir.startswith("mot."):
            continue

        motif_type_path = os.path.join(sequence_library_dir, motif_type_dir)
        if not os.path.isdir(motif_type_path):
            continue

        # Map motif directory names to letters (A, B, C, D)
        motif_type_map = {
            "mot.1": "A",
            "mot.2": "B",
            "mot.3": "C",
            "mot.4": "D",
        }
        motif_letter = motif_type_map.get(motif_type_dir, motif_type_dir)

        logger.info(f"Processing motif type {motif_letter} ({motif_type_dir})")

        # Copy alignment files to organized structure
        for afa_file in os.listdir(motif_type_path):
            if afa_file.endswith(".afa"):
                src = os.path.join(motif_type_path, afa_file)
                # Create descriptive filename: motifA_taxon_id.afa
                base_name = afa_file.replace(".afa", "")
                new_name = f"motif{motif_letter}_{base_name}.afa"
                dst = os.path.join(motif_alignments_dir, new_name)
                shutil.copy2(src, dst)

                # Store metadata
                motif_metadata[new_name.replace(".afa", "")] = {
                    "motif_type": motif_letter,
                    "original_name": afa_file,
                    "taxon": base_name,  # this is a placeholder for if eventually the taxon info is incorporated.
                    "file_path": dst,
                }

    # Save metadata
    metadata_file = os.path.join(motif_dir, "motif_metadata.json")
    with open(metadata_file, "w") as f:
        json.dump(motif_metadata, f, indent=2)

    logger.info(f"Processed {len(motif_metadata)} motif alignments")

    # Create HMM database
    output_hmm = os.path.join(data_dir, "profiles/hmmdbs", "rvmt_motifs.hmm")
    logger.info(f"Building HMM database: {output_hmm}")

    hmmdb_from_directory(
        msa_dir=motif_alignments_dir,
        output=output_hmm,
        msa_pattern="*.afa",
        info_table=None,  # We'll use our metadata file instead
        name_col=None,
        accs_col=None,
        desc_col=None,
    )

    # Create MMseqs profile database
    mmseqs_output = os.path.join(
        data_dir, "profiles/mmseqs_dbs/rvmt_motifs", "rvmt_motifs"
    )
    os.makedirs(os.path.dirname(mmseqs_output), exist_ok=True)
    logger.info(f"Building MMseqs profile database: {mmseqs_output}")

    mmseqs_profile_db_from_directory(
        msa_dir=motif_alignments_dir,
        output=mmseqs_output,
        info_table=None,
        msa_pattern="*.afa",
        name_col=None,
        accs_col=None,
        desc_col=None,
    )

    # place holder for dimanond db creation

    # Clean up extracted directory
    try:
        # move the metadata file out before removing
        shutil.move(metadata_file, os.path.join(data_dir, "profiles"))

        shutil.rmtree(motif_dir)
        os.remove(motif_archive)
        os.remove(os.path.join(motif_dir, "motif_sequence_library.tar.gz"))

    except Exception as e:
        logger.warning(f"Could not remove motif alignments directory: {e}")

    # logger.info(f"RVMT motif preparation completed. Metadata saved to: {metadata_file}")
    logger.info(f"HMM database: {output_hmm}")
    logger.info(f"MMseqs database: {mmseqs_output}")

    return True


def prepare_contamination_seqs(data_dir, threads, logger):
    """Prepare the masking and contamination sequence sets used in sequence filtering.
    The contamination sequences refers to adapters (from bbtools and Fire lab) and rRNA sequences (SILVA + NCBI).
    The masking refers to creating a compressed set of viral sequences (made from concatenating RVMT and NCBI ribovirus, applies entropy masking and compressesion) that can be used for masking potentail viral sequences in host data.

    Args:
        data_dir (str): Base directory for data storage
        threads (int): Number of CPU threads to use
        logger: Logger object for recording progress and errors
    returns:
        None
    """

    logger.info("Preparing masking sequences by combining RVMT and NCBI ribovirus")

    # Create directories (if not already existing)
    contam_dir = os.path.join(data_dir, "contam")
    rrna_dir = os.path.join(contam_dir, "rrna")
    adapter_dir = os.path.join(contam_dir, "adapters")
    masking_dir = os.path.join(contam_dir, "masking")
    os.makedirs(contam_dir, exist_ok=True)
    os.makedirs(rrna_dir, exist_ok=True)
    os.makedirs(adapter_dir, exist_ok=True)
    os.makedirs(masking_dir, exist_ok=True)

    # Masking sequences preparation
    rvmt_fasta_path = os.path.join(
        data_dir, "reference_seqs", "RVMT", "RVMT_cleaned_contigs.fasta.gz"
    )
    ncbi_ribovirus_fasta_path = os.path.join(
        data_dir,
        "reference_seqs",
        "ncbi_virus",
        "refseq_ribovirus_genomes.fasta",
    )

    # Deduplicate directly from multiple files (no concatenation needed)
    deduplicated_fasta = os.path.join(masking_dir, "combined_deduplicated.fasta")
    logger.info(
        f"Deduplicating sequences from {len([rvmt_fasta_path, ncbi_ribovirus_fasta_path])} files"
    )

    stats = remove_duplicates(
        input_file=[rvmt_fasta_path, ncbi_ribovirus_fasta_path],
        output_file=deduplicated_fasta,
        by="seq",
        revcomp_as_distinct=False,  # Treat reverse complement as duplicate
        return_stats=True,
        logger=logger,
    )
    #     (rolypoly_tk) ➜  rolypoly git:(main) ✗ time seqkit rmdup -i  -s code/rolypoly/data/reference_seqs/RVMT/RVMT_cleaned_contigs.fasta   code/rolypoly/data/reference_seqs/ncbi_virus/refseq_ribovirus_genomes.fasta > /dev/null
    # [INFO] 5399 duplicated records removed
    # seqkit rmdup -i -s   > /dev/null  14.70s user 0.53s system 75% cpu 20.095 total
    #     #(rolypoly_tk) ➜  rolypoly git:(main) ✗ time seqkit rmdup --quiet code/rolypoly/data/reference_seqs/RVMT/RVMT_cleaned_contigs.fasta   code/rolypoly/data/reference_seqs/ncbi_virus/refseq_ribovirus_genomes.fasta --quiet | seqkit stats
    # file  format  type  num_seqs        sum_len  min_len  avg_len    max_len
    # -     FASTA   DNA    397,135  1,582,230,847      136  3,984.1  2,473,870
    # seqkit rmdup --quiet   --quiet  0.85s user 0.58s system 10% cpu 13.454 total
    # seqkit stats  10.93s user 0.31s system 83% cpu 13.450 total
    # #In [10]: remove_duplicates(
    # ...:         input_file=[rvmt_fasta_path, ncbi_ribovirus_fasta_path],
    # ...:         output_file=deduplicated_fasta,
    # ...:         by="seq",
    # ...:         revcomp_as_distinct=False,  # Treat reverse complement as duplicate
    # ...:         return_stats=True,
    # ...:         logger=logger
    # ...:     )
    # INFO     2025-11-21 12:26:16 - Processing 2 input files                                                                               sequences.py:451
    # INFO     2025-11-21 12:26:25 - Processed 397135 records: 391736 unique, 5399 duplicates removed                                       sequences.py:597
    # Out[10]: {'total_records': 397135, 'unique_records': 391736, 'duplicates_removed': 5399}

    if stats:
        logger.info(
            f"Deduplication stats: {stats['unique_records']} unique sequences from {stats['total_records']} total, {stats['duplicates_removed']} duplicates removed"
        )

    # Apply entropy masking to the deduplicated sequences
    logger.info("Applying entropy masking to combined sequences")
    entropy_masked_path = os.path.join(
        masking_dir, "combined_entropy_masked.fasta.gz"
    )

    bbmask(
        in1=deduplicated_fasta,
        out=entropy_masked_path,
        entropy=0.1,
        entropywindow=30,
        threads=threads,
    )

    # reduce size with kcompress
    logger.info("Compressing sequences with kcompress")
    compressed_path = os.path.join(masking_dir, "combined_compressed.fasta.gz")

    kcompress(
        in1=entropy_masked_path,
        out=compressed_path,
        fuse=500,
        k=31,
        prealloc=True,
        threads=threads,
    )

    #  now complexity masing again just to be sure
    bbmask(
        in1=compressed_path,
        out=entropy_masked_path,
        entropy=0.2,
        entropywindow=25,
        threads=threads,
    )

    # now a similar process for the orfs
    rvmt_fasta_path = os.path.join(
        data_dir, "reference_seqs", "RVMT", "RVMT_cleaned_orfs.faa.gz"
    )
    ncbi_ribovirus_fasta_path = os.path.join(
        data_dir,
        "reference_seqs",
        "ncbi_virus",
        "refseq_ribovirus_genomes_orfs.faa",
    )

    # Deduplicate directly from multiple files (no concatenation needed)
    deduplicated_fasta = os.path.join(
        masking_dir, "combined_deduplicated_orfs.faa.gz"
    )

    logger.info(
        f"Deduplicating sequences from {len([rvmt_fasta_path, ncbi_ribovirus_fasta_path])} files"
    )

    stats = remove_duplicates(
        input_file=[rvmt_fasta_path, ncbi_ribovirus_fasta_path],
        output_file=deduplicated_fasta,
        by="seq",
        revcomp_as_distinct=False,  # Treat reverse complement as duplicate
        return_stats=True,
        logger=logger,
    )

    # clean up intermediate files
    try:
        os.remove(compressed_path)
    except Exception as e:
        logger.warning(f"Could not remove intermediate files: {e}")

    # Prepare adapter sequences
    logger.info("Fetching adapter sequences")
    fetch_and_extract(
        url="https://raw.githubusercontent.com/bbushnell/BBTools/refs/heads/master/resources/adapters.fa",
        fetched_to=os.path.join(adapter_dir, "bbmap_adapters.fa"),
        rename_extracted=os.path.join(adapter_dir, "bbmap_adapters.fa"),
    )
    fetch_and_extract(
        url="https://raw.githubusercontent.com/FireLabSoftware/CountRabbit/refs/heads/main/illuminatetritis1223wMultiN.fa",
        fetched_to=os.path.join(adapter_dir, "AFire_illuminatetritis1223.fa"),
        rename_extracted=os.path.join(adapter_dir, "AFire_illuminatetritis1223.fa"),
    )
    # remove the poly-monomer from Fire lab adapters
    filter_fasta_by_headers(
        fasta_file=os.path.join(adapter_dir, "AFire_illuminatetritis1223.fa"),
        output_file=os.path.join(adapter_dir, "AFire_illuminatetritis1223_filtered.fa"),
        headers=["A70", "T70"],
        invert=True,
    )
    shutil.move(
        os.path.join(adapter_dir, "AFire_illuminatetritis1223_filtered.fa"),
        os.path.join(adapter_dir, "AFire_illuminatetritis1223.fa"),
    )

    # ===== rRNA Database Preparation with Metadata =====
    # Download and prepare ribosomal RNA sequences (bacterial, archaeal, eukaryotic)
    # Creates a metadata table with taxonomy lineages and FTP download links for host genomes/transcriptomes
    # This eliminates dependencies on taxonkit, ncbi-datasets CLI, and taxdump files

    logger.info("Preparing rRNA database with metadata")
    rrna_dir = os.path.join(contam_dir, "rrna")
    os.makedirs(rrna_dir, exist_ok=True)

    silva_release = "138.2"

    # Download SILVA rRNA sequences (SSU and LSU)
    logger.info(f"Downloading SILVA {silva_release} rRNA sequences")
    silva_ssu_path = os.path.join(
        rrna_dir, f"SILVA_{silva_release}_SSURef_NR99_tax_silva.fasta"
    )
    silva_lsu_path = os.path.join(
        rrna_dir, f"SILVA_{silva_release}_LSURef_NR99_tax_silva.fasta"
    )

    fetch_and_extract(
        f"https://www.arb-silva.de/fileadmin/silva_databases/release_{silva_release.replace('.', '_')}/Exports/SILVA_{silva_release}_SSURef_NR99_tax_silva.fasta.gz",
        fetched_to=os.path.join(rrna_dir, "tmp_ssu.fasta.gz"),
        extract_to=rrna_dir,
        rename_extracted=silva_ssu_path,
        logger=logger,
    )
    fetch_and_extract(
        f"https://www.arb-silva.de/fileadmin/silva_databases/release_{silva_release.replace('.', '_')}/Exports/SILVA_{silva_release}_LSURef_NR99_tax_silva.fasta.gz",
        fetched_to=os.path.join(rrna_dir, "tmp_lsu.fasta.gz"),
        extract_to=rrna_dir,
        rename_extracted=silva_lsu_path,
        logger=logger,
    )

    # Download SILVA taxonomy mappings (maps accessions to NCBI taxids)
    logger.info("Fetching/making SILVA taxonomy mappings (to NCBI taxids)")

    silva_ssu_taxmap = pl.read_csv(
        "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/taxmap_embl-ebi_ena_ssu_ref_nr99_138.2.txt.gz",
        truncate_ragged_lines=True,
        separator="\t",
        infer_schema_length=123123,
    )
    silva_lsu_taxmap = pl.read_csv(
        "https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/ncbi/taxmap_embl-ebi_ena_lsu_ref_nr99_138.2.txt.gz",
        truncate_ragged_lines=True,
        separator="\t",
        infer_schema_length=123123,
    )
    silva_taxmap = pl.concat([silva_lsu_taxmap, silva_ssu_taxmap])

    # Parse SILVA headers and extract accessions
    logger.info("Parsing SILVA sequences and extracting metadata")

    silva_fasta_df = pl.concat(
        [
            from_fastx_eager(silva_ssu_path).with_columns(
                pl.lit("SSU").alias("rRNA_type")
            ),
            from_fastx_eager(silva_lsu_path).with_columns(
                pl.lit("LSU").alias("rRNA_type")
            ),
        ]
    )
    logger.info(f"total SILVA sequences {silva_fasta_df.height}")

    # Extract accession from header (format: >accession.version rest_of_header)
    silva_fasta_df = silva_fasta_df.with_columns(
        primaryAccession=pl.col("header").str.extract(
            r"^([A-Za-z0-9_]+)(?:\.\d+)*", 1
        ),  # DQ150555.1.2478 -> DQ150555
        accession=pl.col("header").str.extract(
            r"^([A-Za-z0-9_]+(?:\.\d+)?)", 1
        ),  # AY846379 or DQ150555.1
        taxonomy_raw=pl.col("header").str.replace(r"^\S+\s+", ""),
    )
    # silva_fasta_df = silva_fasta_df.with_columns(
    #     pl.col("sequence").str.len_chars().alias("seq_length")
    # )
    # silva_taxmap = silva_taxmap.with_columns(
    #     (pl.col("stop") - pl.col("start")).alias("seq_length")
    # )

    silva_df = silva_fasta_df.join(
        silva_taxmap.select(
            ["primaryAccession", "ncbi_taxonid", "submitted_path"]
        ).unique(),  # seq_length
        on=["primaryAccession"],
        how="inner",
    )
    silva_df.write_parquet(os.path.join(rrna_dir, "silva_rrna_sequences.parquet"))
    # silva_df.height
    # silva_df["ncbi_taxonid"].null_count()

    # Load SILVA taxonomy mappings
    logger.info(
        f"Merged taxonomy for {silva_df.filter(pl.col('ncbi_taxonid').is_not_null()).height} SILVA sequences"
    )

    unique_taxids = (
        silva_df.filter(pl.col("ncbi_taxonid").is_not_null())
        .select("ncbi_taxonid")
        .unique()["ncbi_taxonid"]
        .to_list()
    )
    logger.info(
        f"Total of {len(unique_taxids)} unique NCBI taxids found in SILVA sequences"
    )

    # Generate FTP download URLs for host genomes/transcriptomes
    fetch_and_extract(
        url="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt",
        fetched_to=os.path.join(rrna_dir, "assembly_summary_genbank.txt.gz"),
        extract=False,
    )
    logger.info("Loading NCBI GenBank assembly summary")
    # genbank_summary = pl.read_csv(os.path.join(rrna_dir, "assembly_summary_genbank.txt.gz",),
    # infer_schema_length=100020, separator="\t", skip_rows=1,
    # null_values=["na","NA","-"],ignore_errors=True,
    # has_header=True)
    # polars failed me, so using line by line iterator
    from gzip import open as gz_open

    with gz_open(os.path.join(rrna_dir, "assembly_summary_genbank.txt.gz"), "r") as f:
        header = None
        records = []
        i = 0
        for line in f:
            if i == 0:
                i += 1
                continue
            line = line.rstrip(b"\n")
            if i == 1:
                header = line.decode()[1:].strip().split("\t")
                i += 1
                continue
            fields = line.decode().strip().split("\t")
            record = dict(zip(header, fields))
            records.append(record)
    genbank_summary = pl.from_records(records).rename({"taxid": "ncbi_taxonid"})
    genbank_summary.collect_schema()
    # Schema([('assembly_accession', String),
    #         ('bioproject', String),
    #         ('biosample', String),
    #         ('wgs_master', String),
    #         ('refseq_category', String),
    #         ('ncbi_taxonid', String),
    #         ('species_taxid', String),
    #         ('organism_name', String),
    #         ('infraspecific_name', String),
    #         ('isolate', String),
    #         ('version_status', String),
    #         ('assembly_level', String),
    #         ('release_type', String),
    #         ('genome_rep', String),
    #         ('seq_rel_date', String),
    #         ('asm_name', String),
    #         ('asm_submitter', String),
    #         ('gbrs_paired_asm', String),
    #         ('paired_asm_comp', String),
    #         ('ftp_path', String),
    #         ('excluded_from_refseq', String),
    #         ('relation_to_type_material', String),
    #         ('asm_not_live_date', String),
    #         ('assembly_type', String),
    #         ('group', String),
    #         ('genome_size', String),
    #         ('genome_size_ungapped', String),
    #         ('gc_percent', String),
    #         ('replicon_count', String),
    #         ('scaffold_count', String),
    #         ('contig_count', String),
    #         ('annotation_provider', String),
    #         ('annotation_name', String),
    #         ('annotation_date', String),
    #         ('total_gene_count', String),
    #         ('protein_coding_gene_count', String),
    #         ('non_coding_gene_count', String),
    #         ('pubmed_id', String)])

    genbank_summary.write_csv(
        os.path.join(rrna_dir, "genbank_assembly_summary.tsv"), separator="\t"
    )
    genbank_summary = pl.read_csv(
        os.path.join(rrna_dir, "genbank_assembly_summary.tsv"),
        infer_schema_length=100020,
        separator="\t",
        null_values=["na", "NA", "-"],
        ignore_errors=True,
        has_header=True,
    )
    # In [91]: genbank_summary.collect_schema()
    # Out[91]:
    # Schema([('assembly_accession', String),
    #         ('bioproject', String),
    #         ('biosample', String),
    #         ('wgs_master', String),
    #         ('refseq_category', String),
    #         ('ncbi_taxonid', Int64),
    #         ('species_taxid', Int64),
    #         ('organism_name', String),
    #         ('infraspecific_name', String),
    #         ('isolate', String),
    #         ('version_status', String),
    #         ('assembly_level', String),
    #         ('release_type', String),
    #         ('genome_rep', String),
    #         ('seq_rel_date', String),
    #         ('asm_name', String),
    #         ('asm_submitter', String),
    #         ('gbrs_paired_asm', String),
    #         ('paired_asm_comp', String),
    #         ('ftp_path', String),
    #         ('excluded_from_refseq', String),
    #         ('relation_to_type_material', String),
    #         ('asm_not_live_date', String),
    #         ('assembly_type', String),
    #         ('group', String),
    #         ('genome_size', Int64),
    #         ('genome_size_ungapped', Int64),
    #         ('gc_percent', Float64),
    #         ('replicon_count', Int64),
    #         ('scaffold_count', Int64),
    #         ('contig_count', Int64),
    #         ('annotation_provider', String),
    #         ('annotation_name', String),
    #         ('annotation_date', String),
    #         ('total_gene_count', Int64),
    #         ('protein_coding_gene_count', Int64),
    #         ('non_coding_gene_count', Int64),
    #         ('pubmed_id', String)])

    genbank_summary.write_parquet(
        os.path.join(rrna_dir, "genbank_assembly_summary.parquet")
    )
    genbank_summary.write_csv(
        os.path.join(rrna_dir, "genbank_assembly_summary.tsv"), separator="\t"
    )

    # next, for every unique ncbi_taxonid, we select the one that has the most protein_coding_gene_count, then refseq_category, then tie breaking with non_coding_gene_count, tie breaking by latest assembly (by seq_rel_date).
    temp_genbank = genbank_summary.sort(
        by=[
            pl.col("protein_coding_gene_count").cast(pl.Int64).reverse(),
            pl.col("refseq_category").reverse(),
            pl.col("non_coding_gene_count").cast(pl.Int64).reverse(),
            pl.col("seq_rel_date").reverse(),
        ]
    ).unique(subset=["ncbi_taxonid"], keep="first")
    logger.info(
        f"Filtered GenBank summary to {temp_genbank.height} unique taxid entries for SILVA sequences"
    )
    temp_genbank = temp_genbank.filter(
        pl.col("ncbi_taxonid").is_in(unique_taxids)
    ).unique()
    temp_genbank.height
    # only 30482 out ok ~100k?
    fetch_and_extract(
        url="http://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz",
        fetched_to=os.path.join(rrna_dir, "gene2accession.gz"),
        extract=False,
    )
    gene2accession = pl.read_csv(
        os.path.join(rrna_dir, "gene2accession.gz"),
        separator="\t",
        # skip_rows=1,
        # infer_schema_length=100020,
        null_values=["na", "NA", "-"],
        ignore_errors=True,
        has_header=True,
        # n_rows=100
    )
    gene2accession.write_parquet(os.path.join(rrna_dir, "gene2accession.parquet"))
    # gene2accession = pl.read_parquet(os.path.join(rrna_dir, "gene2accession.parquet"))
    # gene2accession.collect_schema()
    # Schema([('#tax_id', Int64),
    #     ('GeneID', Int64),
    #     ('status', String),
    #     ('RNA_nucleotide_accession.version', String),
    #     ('RNA_nucleotide_gi', String),
    #     ('protein_accession.version', String),
    #     ('protein_gi', Int64),
    #     ('genomic_nucleotide_accession.version', String),
    #     ('genomic_nucleotide_gi', Int64),
    #     ('start_position_on_the_genomic_accession', Int64),
    #     ('end_position_on_the_genomic_accession', Int64),
    #     ('orientation', String),
    #     ('assembly', String),
    #     ('mature_peptide_accession.version', String),
    #     ('mature_peptide_gi', String),
    #     ('Symbol', String)])
    gene2accession = gene2accession.rename({"#tax_id": "ncbi_taxonid"})
    test_df = gene2accession.filter(pl.col("ncbi_taxonid").is_in(unique_taxids))
    test_df.height  # 148449745
    test_df2 = gene2accession.select(["ncbi_taxonid", "assembly"]).unique()
    test_df2.height  # 52548

    silva_df = silva_df.with_columns(
        ncbi_taxonid=pl.col("ncbi_taxonid").cast(pl.String)
    )

    silva_df1 = silva_df.join(
        genbank_summary.select(["ncbi_taxonid", "ftp_path"]),
        on=["ncbi_taxonid"],
        how="left",
    )
    silva_df1

    silva_df = silva_df.with_columns(
        genome_ftp_url=pl.when(pl.col("ncbi_taxonid").is_not_null())
        .then(
            pl.format(
                "https://ftp.ncbi.nlm.nih.gov/genomes/all/refseq/taxid_{}/",
                pl.col("ncbi_taxonid"),
            )
        )
        .otherwise(None),
        datasets_api_url=pl.when(pl.col("ncbi_taxonid").is_not_null())
        .then(
            pl.format(
                "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/download?include_annotation_type=GENOME_FASTA,RNA_FASTA",
                pl.col("ncbi_taxonid"),
            )
        )
        .otherwise(None),
    )

    # Save metadata table
    metadata_output = os.path.join(rrna_dir, "rrna_metadata.tsv")
    silva_df.write_csv(metadata_output, separator="\t")
    logger.info(
        f"Saved rRNA metadata table with {len(silva_df)} entries to {metadata_output}"
    )

    # Merge SILVA sequences and apply entropy masking
    logger.info("Merging and masking SILVA sequences")
    silva_merged = os.path.join(rrna_dir, "SILVA_merged.fasta")
    silva_masked = os.path.join(rrna_dir, "SILVA_merged_masked.fasta")

    # Concatenate SILVA files
    run_command_comp(
        base_cmd="cat",
        positional_args=[silva_ssu_path, silva_lsu_path],
        positional_args_location="end",
        params={},
        output_file=silva_merged,
        logger=logger,
    )

    # Apply entropy masking
    bbduk(
        in1=silva_merged,
        out=silva_masked,
        entropy=0.6,
        entropyk=4,
        entropywindow=24,
        maskentropy=True,
        ziplevel=9,
    )

    logger.info(f"Created masked SILVA rRNA database: {silva_masked}")

    # clean up
    try:
        os.remove(deduplicated_fasta)
        os.remove(compressed_path)
    except Exception as e:
        logger.warning(f"Could not remove intermediate files: {e}")

    logger.info(f"Masking sequences prepared in {masking_dir}")


def prepare_trna_data(data_dir, logger):
    trna_dir = os.path.join(data_dir, "tr", "trna")
    file_url = (
        "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/RF00005.fa.gz"
    )
    trna_seqs = os.path.join(trna_dir, "tRNA_sequences.fasta")
    gz_filename = "RF00005.fa.gz"
    deduplicated_fasta = os.path.join(trna_dir, "tRNA_sequences_deduplicated.fasta")

    fetch_and_extract(
        url=file_url,
        fetched_to=os.path.join(trna_dir, gz_filename),
        extract_to=trna_dir,
        expected_file=trna_seqs,
    )
    logger.info(f"Downloaded tRNA sequences to {trna_seqs}")
    # remove duplicates
    remove_duplicates(
        input_file=trna_seqs,
        output_file=deduplicated_fasta,
        return_stats=True,
        by="seq",
    )
    from rolypoly.utils.bio.polars_fastx import fasta_stats
    from rolypoly.utils.bio.sequences import write_fasta_file

    info_table = fasta_stats(deduplicated_fasta)
    info_table = info_table.filter(
        pl.col("length").is_between(60, 250), pl.col("gc_content") >= 0.01
    )

    write_fasta_file(
        seqs=info_table["sequence"].to_list(),
        headers=info_table["header"].to_list(),
        output_file=os.path.join(
            trna_dir, "tRNA_sequences_deduplicated_filtered.fasta"
        ),
    )
    logger.info(
        f"Wrote filtered tRNA sequences to {os.path.join(trna_dir, 'tRNA_sequences_deduplicated_filtered.fasta')}"
    )


def prepare_plastid_data(data_dir, logger):
    """Prepare plastid sequence data for contamination filtering.

    Downloads NCBI RefSeq plastid sequences, combines them, and removes duplicates.

    Args:
        data_dir (str): Base directory for data storage
        logger: Logger object for recording progress and errors
    returns:
        None
    """
    plastid_dir = os.path.join(data_dir, "reference_seqs", "plastid_refseq")
    os.makedirs(plastid_dir, exist_ok=True)

    logger.info("Downloading NCBI RefSeq plastid sequences")

    base_url = "https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid."
    suffix = ".genomic.fna.gz"
    files_to_get = ["1.1", "1.2", "2.1", "2.2", "3.1"]

    all_files = []
    downloaded_files = []

    for version in files_to_get:
        file_url = f"{base_url}{version}{suffix}"
        gz_filename = f"plastid.{version}.genomic.fna.gz"
        fasta_filename = f"plastid.{version}.genomic.fna"

        logger.info(f"Downloading plastid version {version}")

        # Download and extract the file
        try:
            extracted_path = fetch_and_extract(
                url=file_url,
                fetched_to=os.path.join(plastid_dir, gz_filename),
                extract_to=plastid_dir,
                expected_file=fasta_filename,
                logger=logger,
            )
            downloaded_files.append(extracted_path)
            all_files.append(extracted_path)
            all_files.append(os.path.join(plastid_dir, gz_filename))
            logger.info(f"Successfully downloaded and extracted {fasta_filename}")
        except Exception as e:
            logger.error(f"Failed to download plastid version {version}: {e}")
            continue

    if not downloaded_files:
        logger.error("No plastid files were successfully downloaded")
        return

    # Combine and deduplicate the sequences
    combined_fasta = os.path.join(plastid_dir, "combined_plastid_refseq.fasta")
    logger.info(f"Combining and deduplicating {len(downloaded_files)} plastid files")

    stats = remove_duplicates(
        input_file=downloaded_files,
        output_file=combined_fasta,
        by="seq",
        revcomp_as_distinct=False,  # Treat reverse complement as duplicate
        return_stats=True,
        logger=logger,
    )

    if stats:
        logger.info(
            f"Plastid deduplication stats: {stats['unique_records']} unique sequences from {stats['total_records']} total, {stats['duplicates_removed']} duplicates removed"
        )

    # Clean up individual files
    try:
        for file_path in all_files:
            if os.path.exists(file_path):
                os.remove(file_path)
                logger.debug(
                    f"Removed intermediate file: {os.path.basename(file_path)}"
                )
    except Exception as e:
        logger.warning(f"Could not remove intermediate plastid files: {e}")

    logger.info(f"Plastid sequences prepared in {plastid_dir}")


def parse_dmp_line(line: str) -> list[str]:
    """Parse one NCBI ``*.dmp`` line without discarding empty fields."""
    return [field.strip() for field in line.rstrip("\n").split("|")[:-1]]


def format_dmp_line(fields: list[str]) -> str:
    """Render fields using the conventional NCBI taxdump separators."""
    return "\t|\t".join(fields) + "\t|\n"


def renumber_taxdump(
    source_dir: Path,
    output_dir: Path,
    mapping_path: Path | None = None,
) -> dict[str, str]:
    """Copy a taxdump while replacing sparse taxids with dense sequential IDs.

    MMseqs2 allocates taxonomy structures by numeric taxid. Taxonkit-generated
    ICTV IDs can exceed two billion, making a small taxonomy consume gigabytes.
    The root remains ID 1, node order stays deterministic, and the source-to-
    dense mapping is written beside the output taxonomy.
    """
    source_dir = source_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    nodes_path = source_dir / "nodes.dmp"
    names_path = source_dir / "names.dmp"
    if not nodes_path.is_file() or not names_path.is_file():
        raise FileNotFoundError(
            f"Taxdump requires nodes.dmp and names.dmp under {source_dir}"
        )

    node_rows = [
        parse_dmp_line(line) for line in nodes_path.read_text().splitlines(True)
    ]
    if not node_rows:
        raise ValueError(f"Taxdump contains no nodes: {nodes_path}")
    root_taxid = next((row[0] for row in node_rows if row[0] == row[1]), None)
    if root_taxid is None:
        raise ValueError(f"Could not identify a self-parent root in {nodes_path}")

    ordered_taxids = [root_taxid]
    ordered_taxids.extend(row[0] for row in node_rows if row[0] != root_taxid)
    mapping = {
        source_taxid: str(index)
        for index, source_taxid in enumerate(ordered_taxids, start=1)
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "nodes.dmp").open("w", encoding="utf-8") as output:
        for row in node_rows:
            row[0] = mapping[row[0]]
            row[1] = mapping[row[1]]
            output.write(format_dmp_line(row))

    with (
        names_path.open(encoding="utf-8") as source,
        (output_dir / "names.dmp").open("w", encoding="utf-8") as output,
    ):
        for line in source:
            row = parse_dmp_line(line)
            row[0] = mapping[row[0]]
            output.write(format_dmp_line(row))

    merged_path = source_dir / "merged.dmp"
    with (output_dir / "merged.dmp").open("w", encoding="utf-8") as output:
        if merged_path.is_file():
            for line in merged_path.open(encoding="utf-8"):
                row = parse_dmp_line(line)
                if row[0] in mapping and row[1] in mapping:
                    output.write(format_dmp_line([mapping[row[0]], mapping[row[1]]]))

    delnodes_path = source_dir / "delnodes.dmp"
    with (output_dir / "delnodes.dmp").open("w", encoding="utf-8") as output:
        if delnodes_path.is_file():
            for line in delnodes_path.open(encoding="utf-8"):
                row = parse_dmp_line(line)
                if row[0] in mapping:
                    output.write(format_dmp_line([mapping[row[0]]]))

    for source_path in source_dir.iterdir():
        if source_path.is_file() and source_path.name not in {
            "nodes.dmp",
            "names.dmp",
            "merged.dmp",
            "delnodes.dmp",
        }:
            shutil.copy2(source_path, output_dir / source_path.name)

    mapping_path = mapping_path or output_dir / "taxid_map.tsv"
    mapping_path.parent.mkdir(parents=True, exist_ok=True)
    with mapping_path.open("w", encoding="utf-8") as output:
        output.write("source_taxid\tdense_taxid\n")
        for source_taxid in ordered_taxids:
            output.write(f"{source_taxid}\t{mapping[source_taxid]}\n")
    return mapping


def prepare_ncbi_virus_taxonomy_inputs(
    data_dir: str | Path,
    logger: logging.Logger,
    force: bool = False,
) -> tuple[Path, Path, Path]:
    """Build the ICTV taxdump and NCBI-virus-taxid to ICTV-taxid map.

    ClusteredNR membership metadata already carries NCBI taxids, so this step no
    longer downloads or scans the full accession2taxid table. It maps viral NCBI
    taxids to the nearest matching ICTV lineage name and keeps an all-viral taxid
    list for mapping representative taxids from the ClusteredNR SQLite tables.
    """
    output_root = Path(data_dir).expanduser().resolve() / "reference_seqs/ncbi_virus"
    input_root = output_root / "protein_taxdb"
    work_dir = output_root / "build/taxonomy_inputs"
    raw_ictv_taxdump = work_dir / "ictv_taxdump_raw"
    source_taxdump = input_root / "ictv_taxdump"
    ncbi_taxid_to_ictv = input_root / "ncbi_taxid2ictv.tsv"
    viral_taxids_txt = input_root / "ncbi_viral_taxids.txt"
    required = (
        ncbi_taxid_to_ictv,
        viral_taxids_txt,
        source_taxdump / "nodes.dmp",
        source_taxdump / "names.dmp",
    )
    if not force and all(path.is_file() and path.stat().st_size for path in required):
        logger.info("Reusing completed NCBI-virus taxonomy inputs")
        return ncbi_taxid_to_ictv, viral_taxids_txt, source_taxdump

    input_root.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    ictv_xlsx = work_dir / "ictv_vmr.xlsx"
    ictv_taxonomy_tsv = work_dir / "ictv_taxonomy.tsv"
    if force or not ictv_xlsx.is_file():
        simple_fetch("https://ictv.global/vmr/current", ictv_xlsx, logger=logger)

    ictv_df = pl.read_excel(ictv_xlsx, sheet_id=2)
    lineage_df = (
        ictv_df.select([rank.capitalize() for rank in ICTV_RANKS])
        .with_columns(pl.all().cast(pl.Utf8).str.strip_chars())
        .unique()
        .rename({rank.capitalize(): rank for rank in ICTV_RANKS})
    )
    lineage_df.write_csv(ictv_taxonomy_tsv, separator="\t")
    if force or not (raw_ictv_taxdump / "nodes.dmp").is_file():
        subprocess.run(
            [
                "taxonkit",
                "create-taxdump",
                "--out-dir",
                str(raw_ictv_taxdump),
                "--force",
                str(ictv_taxonomy_tsv),
            ],
            check=True,
        )
    renumber_taxdump(
        raw_ictv_taxdump,
        source_taxdump,
        source_taxdump / "taxid_map.tsv",
    )

    ncbi_taxdump_dir = work_dir / "ncbi_taxdump"
    if force or not (ncbi_taxdump_dir / "nodes.dmp").is_file():
        fetch_and_extract(
            url="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            fetched_to=str(work_dir / "taxdump.tar.gz"),
            extract_to=ncbi_taxdump_dir,
            logger=logger,
        )
    if force or not viral_taxids_txt.is_file():
        with viral_taxids_txt.open("w", encoding="utf-8") as output:
            subprocess.run(
                [
                    "taxonkit",
                    "list",
                    "--ids",
                    "10239",
                    "--data-dir",
                    str(ncbi_taxdump_dir),
                    "--indent",
                    "",
                ],
                stdout=output,
                text=True,
                check=True,
            )
    ncbi_taxids = pl.read_csv(
        viral_taxids_txt,
        has_header=False,
        new_columns=["taxid"],
        schema_overrides={"taxid": pl.String},
    ).filter(pl.col("taxid").str.contains(r"^\d+$"))
    logger.info(f"Found {ncbi_taxids.height:,} NCBI viral taxids")

    lineage_out_tsv = work_dir / "ncbi_lineage_ranked.tsv"
    if force or not lineage_out_tsv.is_file():
        with lineage_out_tsv.open("w", encoding="utf-8") as output:
            subprocess.run(
                [
                    "taxonkit",
                    "lineage",
                    "--data-dir",
                    str(ncbi_taxdump_dir),
                    "-R",
                    "-t",
                    str(viral_taxids_txt),
                ],
                stdout=output,
                text=True,
                check=True,
            )

    lineage_raw = pl.read_csv(
        lineage_out_tsv,
        separator="\t",
        has_header=False,
        new_columns=["taxid", "name_lineage", "taxid_lineage", "rank_lineage"],
        schema_overrides={
            "taxid": pl.String,
            "name_lineage": pl.String,
            "taxid_lineage": pl.String,
            "rank_lineage": pl.String,
        },
    ).drop_nulls(subset=["name_lineage", "rank_lineage"])
    ranked_long = (
        lineage_raw.with_columns(
            pl.col("name_lineage").str.split(";").alias("name"),
            pl.col("rank_lineage").str.split(";").alias("rank"),
        )
        .explode(["name", "rank"])
        .filter(pl.col("rank").is_in(ICTV_RANKS))
    )
    wide_lineage = ranked_long.pivot(
        index="taxid",
        on="rank",
        values="name",
        aggregate_function="first",
    )
    for rank in ICTV_RANKS:
        if rank not in wide_lineage.columns:
            wide_lineage = wide_lineage.with_columns(
                pl.lit(None, dtype=pl.String).alias(rank)
            )

    ictv_name_rows = []
    with (source_taxdump / "names.dmp").open(encoding="utf-8") as names:
        for line in names:
            fields = parse_dmp_line(line)
            if len(fields) >= 4 and fields[3] == "scientific name":
                ictv_name_rows.append({"ictv_taxid": fields[0], "name": fields[1]})
    ictv_names = pl.DataFrame(
        ictv_name_rows,
        schema={"ictv_taxid": pl.String, "name": pl.String},
    ).unique(subset=["name"], keep="first")

    candidate_columns = []
    for rank in reversed(ICTV_RANKS):
        column = f"{rank}_ictv_taxid"
        wide_lineage = wide_lineage.join(
            ictv_names.rename({"name": rank, "ictv_taxid": column}),
            on=rank,
            how="left",
        )
        candidate_columns.append(column)
    ncbi_to_ictv = (
        wide_lineage.with_columns(pl.coalesce(candidate_columns).alias("ictv_taxid"))
        .drop_nulls(subset=["ictv_taxid"])
        .select("taxid", "ictv_taxid")
        .unique()
        .sort("taxid")
    )
    logger.info(
        f"Mapped {ncbi_to_ictv.height:,}/{ncbi_taxids.height:,} "
        "NCBI viral taxids to ICTV"
    )
    ncbi_to_ictv.write_csv(ncbi_taxid_to_ictv, separator="\t")
    return ncbi_taxid_to_ictv, viral_taxids_txt, source_taxdump


def outputs_are_current(paths: Iterable[Path], inputs: Iterable[Path] = ()) -> bool:
    """Return whether outputs are nonempty and no older than their inputs."""
    output_paths = tuple(paths)
    input_paths = tuple(inputs)
    if not output_paths or not all(
        path.is_file() and path.stat().st_size for path in output_paths
    ):
        return False
    if not input_paths:
        return True
    if not all(path.is_file() for path in input_paths):
        return False
    newest_input = max(path.stat().st_mtime_ns for path in input_paths)
    return min(path.stat().st_mtime_ns for path in output_paths) >= newest_input


def prepare_ncbi_virus_taxdb(
    data_dir: str | Path,
    threads: int,
    logger: logging.Logger,
    clustered_nr_db: str | Path | None = None,
    clustered_nr_metadata_db: str | Path | None = None,
    clustered_nr_taxonomy_db: str | Path | None = None,
    clustered_nr_fetch_dir: str | Path | None = None,
    work_dir: str | Path | None = None,
    temp_dir: str | Path | None = None,
    bgzip_threads: int = 1,
    force: bool = False,
    allow_missing_taxonomy: bool = False,
) -> dict[str, Path]:
    """Build the NCBI-virus protein taxonomy bundle from ClusteredNR.

    The step can use an existing decompressed NCBI directory or fetch/extract
    the official tarballs. It queries ``nr_cluster_seq.sqlite3`` through Polars
    with the viral taxid list from the taxonomy-input stage, extracts official
    representative accessions with ``blastdbcmd``, and builds both MMseqs2 and
    DIAMOND from that representative FASTA.
    """
    import taxopy

    data_root = Path(data_dir).expanduser().resolve()
    output_root = data_root / "reference_seqs/ncbi_virus"
    work_dir = (
        Path(work_dir).expanduser().resolve()
        if work_dir
        else output_root / "build/clustered_nr"
    )
    temp_dir = (
        Path(temp_dir).expanduser().resolve() if temp_dir else work_dir / "backend_tmp"
    )
    work_dir.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    blast_db = Path(clustered_nr_db).expanduser().resolve() if clustered_nr_db else None
    if blast_db and blast_db.is_dir():
        blast_db = blast_db / "nr_cluster_seq"
    metadata_db = (
        Path(clustered_nr_metadata_db).expanduser().resolve()
        if clustered_nr_metadata_db
        else next(
            (
                candidate
                for candidate in (
                    Path(str(blast_db) + ".sqlite3"),
                    blast_db.parent / "cluster_data.sqlite3",
                )
                if candidate.is_file()
            ),
            Path(str(blast_db) + ".sqlite3"),
        )
        if blast_db is not None
        else None
    )
    taxonomy_db = (
        Path(clustered_nr_taxonomy_db).expanduser().resolve()
        if clustered_nr_taxonomy_db
        else metadata_db.parent / "taxonomy4blast.sqlite3"
        if metadata_db is not None
        else None
    )
    if blast_db is None or metadata_db is None or taxonomy_db is None:
        fetch_dir = (
            Path(clustered_nr_fetch_dir).expanduser().resolve()
            if clustered_nr_fetch_dir
            else work_dir / "downloads"
        )
        extract_dir = work_dir / "ncbi_clustered_nr"
        fetch_dir.mkdir(parents=True, exist_ok=True)
        extract_dir.mkdir(parents=True, exist_ok=True)

        archive_names = [f"nr_cluster_seq.{index:03d}.tar.gz" for index in range(160)]
        archive_paths = [fetch_dir / name for name in archive_names]
        checksum_paths = [fetch_dir / f"{name}.md5" for name in archive_names]
        required_fetches = archive_paths + checksum_paths
        if clustered_nr_fetch_dir:
            missing = [path.name for path in required_fetches if not path.is_file()]
            if missing:
                preview = ", ".join(missing[:5])
                suffix = f" and {len(missing) - 5} more" if len(missing) > 5 else ""
                raise FileNotFoundError(
                    f"Pre-fetched ClusteredNR directory {fetch_dir} is incomplete: "
                    f"missing {preview}{suffix}"
                )
            logger.info(f"Using pre-fetched ClusteredNR archives from {fetch_dir}")
        elif not all(
            path.is_file() and path.stat().st_size for path in required_fetches
        ):
            if shutil.which("aria2c") is None:
                raise RuntimeError("aria2c is required to download ClusteredNR")
            base_url = "https://ftp.ncbi.nlm.nih.gov/blast/db/experimental"
            url_file = fetch_dir / "clustered_nr_urls.txt"
            urls = [
                f"{base_url}/{filename}"
                for archive_name in archive_names
                for filename in (archive_name, f"{archive_name}.md5")
            ]
            urls.extend(
                [
                    f"{base_url}/nr_cluster_seq-prot-metadata.json",
                    f"{base_url}/README.md",
                ]
            )
            url_file.write_text("\n".join(urls) + "\n", encoding="utf-8")
            logger.info(f"Downloading ClusteredNR archives to {fetch_dir}")
            subprocess.run(
                [
                    "aria2c",
                    f"--dir={fetch_dir}",
                    f"--input-file={url_file}",
                    "--continue=true",
                    "--max-concurrent-downloads=3",
                    "--split=4",
                    "--max-connection-per-server=4",
                    "--min-split-size=64M",
                    "--retry-wait=30",
                    "--max-tries=0",
                    "--auto-file-renaming=false",
                    "--allow-overwrite=true",
                ],
                check=True,
            )

        missing = [path.name for path in required_fetches if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                f"ClusteredNR fetch did not produce {len(missing)} required files"
            )

        verified_dir = work_dir / "verified_archives"
        verified_dir.mkdir(parents=True, exist_ok=True)
        logger.info("Verifying ClusteredNR archive checksums")
        for archive_path, checksum_path in zip(archive_paths, checksum_paths):
            verified = verified_dir / f"{archive_path.name}.complete"
            if outputs_are_current((verified,), (archive_path, checksum_path)):
                continue
            checksum_fields = checksum_path.read_text(encoding="utf-8").split()
            if not checksum_fields:
                raise ValueError(f"Empty checksum file: {checksum_path}")
            expected_checksum = checksum_fields[0].lower()
            digest = hashlib.md5()
            with archive_path.open("rb") as archive_file:
                for chunk in iter(lambda: archive_file.read(8 * 1024 * 1024), b""):
                    digest.update(chunk)
            if digest.hexdigest() != expected_checksum:
                verified.unlink(missing_ok=True)
                raise ValueError(f"Checksum mismatch for {archive_path}")
            verified.write_text(expected_checksum + "\n", encoding="utf-8")

        extracted_dir = extract_dir / ".extracted_archives"
        extracted_dir.mkdir(parents=True, exist_ok=True)
        pending_extractions = []
        for archive_path in archive_paths:
            extracted = extracted_dir / f"{archive_path.name}.complete"
            if outputs_are_current((extracted,), (archive_path,)):
                continue
            pending_extractions.append((archive_path, extracted))

        if pending_extractions:
            extraction_threads = min(threads, len(pending_extractions))
            logger.info(
                f"Extracting {len(pending_extractions)} ClusteredNR archives "
                f"with {extraction_threads} threads"
            )
            with ThreadPoolExecutor(max_workers=extraction_threads) as executor:
                futures = {
                    executor.submit(extract_tar, archive_path, extract_dir, logger): (
                        archive_path,
                        extracted,
                    )
                    for archive_path, extracted in pending_extractions
                }
                for future in as_completed(futures):
                    archive_path, extracted = futures[future]
                    future.result()
                    extracted.write_text("complete\n", encoding="utf-8")
                    logger.info(f"Extracted {archive_path.name}")

        if blast_db is None:
            blast_db = extract_dir / "nr_cluster_seq"
        if metadata_db is None:
            metadata_db = extract_dir / "nr_cluster_seq.sqlite3"
        if taxonomy_db is None:
            taxonomy_db = extract_dir / "taxonomy4blast.sqlite3"

    if blast_db is not None and not list(blast_db.parent.glob(f"{blast_db.name}*.pin")):
        raise FileNotFoundError(f"ClusteredNR BLAST database not found: {blast_db}")
    if metadata_db is not None and not metadata_db.is_file():
        raise FileNotFoundError(
            f"ClusteredNR metadata database not found: {metadata_db}"
        )
    if taxonomy_db is not None and not taxonomy_db.is_file():
        raise FileNotFoundError(
            f"ClusteredNR taxonomy database not found: {taxonomy_db}"
        )

    ncbi_taxid_to_ictv, viral_taxids_txt, source_taxdump = (
        prepare_ncbi_virus_taxonomy_inputs(data_root, logger, force=force)
    )

    final_fasta = work_dir / "ictv_clustered_nr_representatives.fasta.bgz"
    source_metadata = work_dir / "clustered_nr_blastdb_metadata.json"
    representative_viral_taxids = (
        work_dir / "clustered_nr_representative_viral_member_taxids.tsv"
    )
    representative_mapping = work_dir / "clustered_nr_viral_representatives_ictv.tsv"
    representatives = work_dir / "clustered_nr_representatives.txt"
    extracted_representatives = work_dir / "clustered_nr_extracted_representatives.tsv"
    missing_representatives = work_dir / "clustered_nr_missing_representatives.tsv"

    representative_viral_taxids_columns = set()
    if representative_viral_taxids.is_file():
        try:
            representative_viral_taxids_columns = set(
                pl.read_csv(
                    representative_viral_taxids,
                    separator="\t",
                    n_rows=0,
                ).columns
            )
        except Exception:
            representative_viral_taxids_columns = set()
    representative_mapping_columns = set()
    if representative_mapping.is_file():
        try:
            representative_mapping_columns = set(
                pl.read_csv(
                    representative_mapping,
                    separator="\t",
                    n_rows=0,
                ).columns
            )
        except Exception:
            representative_mapping_columns = set()

    if (
        force
        or not {
            "representative_accession",
            "member_accession",
            "ncbi_taxid",
        }.issubset(representative_viral_taxids_columns)
        or not outputs_are_current(
            (representative_viral_taxids, representatives),
            (metadata_db, viral_taxids_txt),
        )
    ):
        if metadata_db is None or not metadata_db.is_file():
            raise FileNotFoundError(
                "ClusteredNR membership SQLite database is required; pass "
                "--clustered-nr-metadata-db or place nr_cluster_seq.sqlite3 beside "
                "--clustered-nr-db"
            )
        logger.info(
            "Selecting ClusteredNR representatives from clusters with viral member taxids"
        )
        viral_taxids = (
            pl.read_csv(
                viral_taxids_txt,
                has_header=False,
                new_columns=["taxid"],
                schema_overrides={"taxid": pl.String},
            )
            .filter(pl.col("taxid").str.contains(r"^\d+$"))
            .select(pl.col("taxid").cast(pl.Int64))
            .unique()
            .sort("taxid")
        )
        representative_sql = """
SELECT DISTINCT
    R.accession AS representative_accession,
    C.member_accession AS member_accession,
    CAST(C.member_taxid AS TEXT) AS ncbi_taxid
FROM ClusterInfo C INDEXED BY ClusterInfoIdx_MembTaxid
JOIN Representative R ON C.representative_id = R.id
WHERE C.member_taxid IN (SELECT taxid FROM ViralTaxid)
"""
        chunk_dir = work_dir / "clustered_nr_sqlite_chunks"
        if chunk_dir.exists():
            shutil.rmtree(chunk_dir)
        chunk_dir.mkdir(parents=True, exist_ok=True)
        with sqlite3.connect(f"file:{metadata_db}?mode=ro", uri=True) as connection:
            connection.execute("PRAGMA temp_store = MEMORY")
            connection.execute("CREATE TEMP TABLE ViralTaxid(taxid INTEGER PRIMARY KEY)")
            connection.executemany(
                "INSERT OR IGNORE INTO ViralTaxid(taxid) VALUES (?)",
                viral_taxids.iter_rows(),
            )
            main_tables = set(
                pl.read_database(
                    "SELECT name FROM sqlite_master WHERE type IN ('table', 'view')",
                    connection,
                )["name"].to_list()
            )
            missing_schema = {
                "ClusterInfo",
                "Representative",
            }.difference(main_tables)
            if missing_schema:
                raise RuntimeError(
                    "Unexpected ClusteredNR SQLite schema. Missing "
                    f"{sorted(missing_schema)}"
                )
            chunk_count = 0
            selected_rows = 0
            for batch in pl.read_database(
                representative_sql,
                connection,
                iter_batches=True,
                batch_size=1_000_000,
            ):
                selected = (
                    batch.select(
                        pl.col("representative_accession").cast(pl.String),
                        pl.col("member_accession").cast(pl.String),
                        pl.col("ncbi_taxid").cast(pl.String),
                    )
                    .drop_nulls()
                    .unique()
                )
                if selected.is_empty():
                    continue
                selected_rows += selected.height
                selected.write_parquet(
                    chunk_dir / f"representatives_{chunk_count:06d}.parquet"
                )
                chunk_count += 1

        if chunk_count == 0:
            shutil.rmtree(chunk_dir)
            raise RuntimeError(
                "ClusteredNR metadata contained no representatives with viral members"
            )

        partial_representative_viral_taxids = representative_viral_taxids.with_suffix(
            ".tsv.partial"
        )
        partial_representatives = representatives.with_suffix(".txt.partial")
        pl.scan_parquet(str(chunk_dir / "*.parquet")).unique().sort(
            "representative_accession", "member_accession", "ncbi_taxid"
        ).sink_csv(partial_representative_viral_taxids, separator="\t")
        (
            pl.scan_csv(partial_representative_viral_taxids, separator="\t")
            .select("representative_accession")
            .unique()
            .sort("representative_accession")
            .sink_csv(partial_representatives, include_header=False)
        )
        partial_representative_viral_taxids.replace(representative_viral_taxids)
        if (
            representatives.is_file()
            and representatives.read_bytes() == partial_representatives.read_bytes()
        ):
            partial_representatives.unlink()
        else:
            partial_representatives.replace(representatives)
        shutil.rmtree(chunk_dir)
        representative_count = (
            pl.scan_csv(representatives, has_header=False, new_columns=["accession"])
            .select(pl.len())
            .collect(engine="streaming")
            .item()
        )
        logger.info(
            f"Selected {representative_count:,} ClusteredNR representatives "
            f"from {selected_rows:,} viral member taxid rows"
        )
    else:
        logger.info("Reusing completed ClusteredNR representative selection")

    if (
        force
        or not {
            "representative_accession",
            "member_accession",
            "member_is_refseq",
            "ncbi_taxid",
            "ictv_taxid",
        }.issubset(representative_mapping_columns)
        or not outputs_are_current(
            (representative_mapping,),
            (representative_viral_taxids, ncbi_taxid_to_ictv),
        )
    ):
        partial_representative_mapping = representative_mapping.with_suffix(
            ".tsv.partial"
        )
        (
            pl.scan_csv(
                representative_viral_taxids,
                separator="\t",
                schema_overrides={
                    "representative_accession": pl.String,
                    "member_accession": pl.String,
                    "ncbi_taxid": pl.String,
                },
            )
            .with_columns(
                pl.col("member_accession")
                .str.contains(r"^[A-Z]{2}_[A-Z0-9]+(?:\.\d+)?$")
                .fill_null(False)
                .alias("member_is_refseq")
            )
            .join(
                pl.scan_csv(
                    ncbi_taxid_to_ictv,
                    separator="\t",
                    schema_overrides={"taxid": pl.String, "ictv_taxid": pl.String},
                ).rename({"taxid": "ncbi_taxid"}),
                on="ncbi_taxid",
                how="left",
            )
            .select(
                "representative_accession",
                "member_accession",
                "member_is_refseq",
                "ncbi_taxid",
                "ictv_taxid",
            )
            .sort(
                ["representative_accession", "member_is_refseq", "member_accession"],
                descending=[False, True, False],
            )
            .sink_csv(partial_representative_mapping, separator="\t")
        )
        partial_representative_mapping.replace(representative_mapping)
    else:
        logger.info("Reusing completed representative-to-ICTV mapping")

    if force or not outputs_are_current(
        (
            final_fasta,
            source_metadata,
            extracted_representatives,
            missing_representatives,
        ),
        (representatives,),
    ):
        if blast_db is None:
            raise ValueError(
                "--clustered-nr-db (or NCBI_CLUSTERED_NR_DB) is required "
                "because the representative FASTA cannot be reused"
            )
        if shutil.which("blastdbcmd") is None:
            raise RuntimeError("blastdbcmd is required for ClusteredNR extraction")
        if shutil.which("bgzip") is None:
            raise RuntimeError("bgzip is required for ClusteredNR extraction")
        logger.info(f"Extracting ClusteredNR representatives from {blast_db}")
        subprocess.run(
            ["blastdbcmd", "-db", str(blast_db), "-dbtype", "prot", "-info"],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        partial_metadata = source_metadata.with_suffix(".json.partial")
        with partial_metadata.open("w", encoding="utf-8") as output:
            subprocess.run(
                ["blastdbcmd", "-db", str(blast_db), "-dbtype", "prot", "-metadata"],
                check=True,
                stdout=output,
                text=True,
            )
        partial_metadata.replace(source_metadata)

        partial_fasta = final_fasta.with_name(final_fasta.name + ".partial")
        partial_extracted = extracted_representatives.with_suffix(".tsv.partial")
        partial_missing = missing_representatives.with_suffix(".tsv.partial")
        log_path = work_dir / "clustered_nr_extract.log"
        chunk_dir = work_dir / "clustered_nr_extract_chunks"
        shutil.rmtree(chunk_dir, ignore_errors=True)
        chunk_dir.mkdir(parents=True, exist_ok=True)

        with representatives.open(encoding="utf-8") as accessions:
            representative_total = sum(
                1 for accession in accessions if accession.strip()
            )
        if representative_total == 0:
            raise RuntimeError(
                f"No ClusteredNR representatives were written to {representatives}"
            )

        extraction_jobs = min(max(1, int(threads)), representative_total)
        chunk_paths = [
            chunk_dir / f"representatives_{index:03d}.txt"
            for index in range(extraction_jobs)
        ]
        chunk_files = [path.open("w", encoding="utf-8") for path in chunk_paths]
        try:
            chunk_index = 0
            with representatives.open(encoding="utf-8") as accessions:
                for accession in accessions:
                    accession = accession.strip()
                    if not accession:
                        continue
                    chunk_files[chunk_index].write(f"{accession}\n")
                    chunk_index = (chunk_index + 1) % extraction_jobs
        finally:
            for chunk_file in chunk_files:
                chunk_file.close()
        chunk_paths = [path for path in chunk_paths if path.stat().st_size > 0]
        logger.info(
            f"Extracting {representative_total:,} ClusteredNR representatives "
            f"with {len(chunk_paths)} parallel blastdbcmd jobs"
        )

        processes = []
        try:
            for chunk_path in chunk_paths:
                fasta_path = chunk_path.with_suffix(".fasta")
                chunk_log = chunk_path.with_suffix(".log")
                stdout = fasta_path.open("w", encoding="utf-8")
                stderr = chunk_log.open("w", encoding="utf-8")
                try:
                    process = subprocess.Popen(
                        [
                            "blastdbcmd",
                            "-db",
                            str(blast_db),
                            "-dbtype",
                            "prot",
                            "-entry_batch",
                            str(chunk_path),
                            "-outfmt",
                            "%f",
                        ],
                        stdout=stdout,
                        stderr=stderr,
                        text=True,
                    )
                except Exception:
                    stdout.close()
                    stderr.close()
                    raise
                processes.append(
                    (process, stdout, stderr, chunk_path, fasta_path, chunk_log)
                )
        except Exception:
            for (
                process,
                stdout,
                stderr,
                _chunk_path,
                _fasta_path,
                _chunk_log,
            ) in processes:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
                stdout.close()
                stderr.close()
            raise

        fasta_parts = []
        failed_chunks = []
        extracted_count = 0
        missing_count = 0
        skipped_prefix = "Error: [blastdbcmd] Skipped "
        allowed_missing_error = (
            "Error: [blastdbcmd] Entry or entries not found in BLAST database"
        )
        partial_extracted.unlink(missing_ok=True)
        partial_missing.unlink(missing_ok=True)
        with (
            partial_extracted.open("w", encoding="utf-8") as extracted,
            partial_missing.open("w", encoding="utf-8") as missing_output,
        ):
            extracted.write("representative_accession\tfasta_accession\n")
            missing_output.write("representative_accession\treason\n")
            for (
                process,
                stdout,
                stderr,
                chunk_path,
                fasta_path,
                chunk_log,
            ) in processes:
                code = process.wait()
                stdout.close()
                stderr.close()

                skipped_accessions = []
                unexpected_log_lines = []
                with chunk_log.open(encoding="utf-8") as log_input:
                    for log_line in log_input:
                        log_line = log_line.strip()
                        if not log_line:
                            continue
                        if log_line.startswith(skipped_prefix):
                            skipped_accessions.append(
                                log_line.removeprefix(skipped_prefix)
                            )
                        elif log_line != allowed_missing_error:
                            unexpected_log_lines.append(log_line)

                if unexpected_log_lines:
                    failed_chunks.append(
                        (
                            fasta_path,
                            chunk_log,
                            "; ".join(unexpected_log_lines[:3]),
                        )
                    )
                    continue
                if code and not skipped_accessions:
                    failed_chunks.append((fasta_path, chunk_log, f"exit {code}"))
                    continue

                skipped_in_chunk = set(skipped_accessions)
                with chunk_path.with_suffix(".missing.tsv").open(
                    "w", encoding="utf-8"
                ) as chunk_missing:
                    for skipped_accession in skipped_accessions:
                        chunk_missing.write(
                            f"{skipped_accession}\tnot_found_in_blastdb\n"
                        )
                        missing_output.write(
                            f"{skipped_accession}\tnot_found_in_blastdb\n"
                        )
                missing_count += len(skipped_accessions)

                expected_accessions = []
                with fasta_path.with_suffix(".extracted.tsv").open(
                    "w", encoding="utf-8"
                ) as chunk_extracted:
                    with chunk_path.open(encoding="utf-8") as accessions:
                        for accession in accessions:
                            accession = accession.strip()
                            if accession and accession not in skipped_in_chunk:
                                expected_accessions.append(accession)

                    fasta_accessions = []
                    if fasta_path.stat().st_size > 0:
                        with fasta_path.open(encoding="utf-8") as fasta_input:
                            for fasta_line in fasta_input:
                                if fasta_line.startswith(">"):
                                    fasta_accessions.append(
                                        fasta_line[1:].split(maxsplit=1)[0]
                                    )

                    if len(expected_accessions) != len(fasta_accessions):
                        failed_chunks.append(
                            (
                                fasta_path,
                                chunk_log,
                                "extracted FASTA count did not match requested "
                                f"non-skipped accessions "
                                f"({len(fasta_accessions):,} != "
                                f"{len(expected_accessions):,})",
                            )
                        )
                        continue

                    for representative, fasta_accession in zip(
                        expected_accessions, fasta_accessions, strict=True
                    ):
                        chunk_extracted.write(
                            f"{representative}\t{fasta_accession}\n"
                        )
                        extracted.write(f"{representative}\t{fasta_accession}\n")
                    extracted_count += len(fasta_accessions)
                    if fasta_accessions:
                        fasta_parts.append(fasta_path)
        if failed_chunks:
            partial_fasta.unlink(missing_ok=True)
            partial_extracted.unlink(missing_ok=True)
            partial_missing.unlink(missing_ok=True)
            failed_summary = ", ".join(
                f"{fasta_path.name} ({reason}; log {chunk_log.name})"
                for fasta_path, chunk_log, reason in failed_chunks[:5]
            )
            raise RuntimeError(
                f"ClusteredNR FASTA extraction failed for {len(failed_chunks)} "
                f"chunk(s): {failed_summary}; see {chunk_dir}"
            )
        if not fasta_parts:
            partial_extracted.unlink(missing_ok=True)
            partial_missing.unlink(missing_ok=True)
            raise RuntimeError(
                f"No ClusteredNR representative sequences were extracted; see {chunk_dir}"
            )
        if missing_count:
            logger.warning(
                f"Skipped {missing_count:,} ClusteredNR representatives absent from "
                f"the BLAST database; see {missing_representatives}"
            )

        with (
            partial_fasta.open("wb") as compressed,
            log_path.open("w", encoding="utf-8") as log,
        ):
            log.write(
                f"Compressing {len(fasta_parts)} extracted ClusteredNR FASTA chunks\n"
            )
            compressor = subprocess.Popen(
                ["bgzip", "-@", str(max(1, int(bgzip_threads))), "-c"],
                stdin=subprocess.PIPE,
                stdout=compressed,
                stderr=log,
            )
            if compressor.stdin is None:
                partial_fasta.unlink(missing_ok=True)
                partial_extracted.unlink(missing_ok=True)
                partial_missing.unlink(missing_ok=True)
                raise RuntimeError("Could not open the bgzip input stream")
            try:
                for fasta_path in fasta_parts:
                    with fasta_path.open("rb") as fasta_input:
                        shutil.copyfileobj(
                            fasta_input, compressor.stdin, length=8 * 1024 * 1024
                        )
                compressor.stdin.close()
            except Exception:
                if compressor.poll() is None:
                    compressor.kill()
                    compressor.wait()
                partial_fasta.unlink(missing_ok=True)
                partial_extracted.unlink(missing_ok=True)
                partial_missing.unlink(missing_ok=True)
                raise
            compressor_code = compressor.wait()
        if compressor_code:
            partial_fasta.unlink(missing_ok=True)
            partial_extracted.unlink(missing_ok=True)
            partial_missing.unlink(missing_ok=True)
            raise RuntimeError(
                f"ClusteredNR FASTA compression failed; see {log_path}"
            )
        try:
            subprocess.run(["bgzip", "-t", str(partial_fasta)], check=True)
        except subprocess.CalledProcessError:
            partial_fasta.unlink(missing_ok=True)
            partial_extracted.unlink(missing_ok=True)
            partial_missing.unlink(missing_ok=True)
            raise
        partial_fasta.replace(final_fasta)
        partial_extracted.replace(extracted_representatives)
        partial_missing.replace(missing_representatives)
        logger.info(
            f"Extracted {extracted_count:,} ClusteredNR representative sequences"
        )
        shutil.rmtree(chunk_dir)
    else:
        logger.info("Reusing extracted ClusteredNR representative FASTA")

    mapping_dir = output_root / "build/taxonomy_mappings"
    mapping_dir.mkdir(parents=True, exist_ok=True)
    taxonomy_dir = output_root / "taxonomy"
    assignments_path = mapping_dir / "primary_sequence_taxonomy.tsv"
    mmseqs_mapping = mapping_dir / "accession2taxid_mmseqs.tsv"
    diamond_mapping = mapping_dir / "accession2taxid_diamond.tsv"
    audit_path = mapping_dir / "clustered_nr_taxonomy_provenance.tsv"
    assignment_outputs = (
        assignments_path,
        mmseqs_mapping,
        diamond_mapping,
        audit_path,
        taxonomy_dir / "nodes.dmp",
        taxonomy_dir / "names.dmp",
    )
    taxonomy_inputs = (
        representative_mapping,
        extracted_representatives,
        source_taxdump / "nodes.dmp",
        source_taxdump / "names.dmp",
    )
    audit_columns = set()
    if audit_path.is_file():
        try:
            audit_columns = set(
                pl.read_csv(audit_path, separator="\t", n_rows=0).columns
            )
        except Exception:
            audit_columns = set()
    if (
        force
        or not {
            "assignment_scope",
            "refseq_viral_member_taxids",
            "refseq_mapped_ictv_taxids",
        }.issubset(audit_columns)
        or not outputs_are_current(assignment_outputs, taxonomy_inputs)
    ):
        taxid_remap = renumber_taxdump(
            source_taxdump,
            taxonomy_dir,
            taxonomy_dir / "taxid_map.tsv",
        )
        taxdb = taxopy.TaxDb(
            nodes_dmp=str(taxonomy_dir / "nodes.dmp"),
            names_dmp=str(taxonomy_dir / "names.dmp"),
            keep_files=True,
        )
        grouped = (
            pl.read_csv(
                representative_mapping,
                separator="\t",
                schema_overrides={
                    "representative_accession": pl.String,
                    "member_accession": pl.String,
                    "member_is_refseq": pl.Boolean,
                    "ncbi_taxid": pl.String,
                    "ictv_taxid": pl.String,
                },
            )
            .join(
                pl.read_csv(
                    extracted_representatives,
                    separator="\t",
                    schema_overrides={
                        "representative_accession": pl.String,
                        "fasta_accession": pl.String,
                    },
                ),
                on="representative_accession",
                how="inner",
            )
            .with_columns(
                pl.when(pl.col("ictv_taxid") == "")
                .then(None)
                .otherwise(pl.col("ictv_taxid"))
                .alias("ictv_taxid"),
                pl.col("member_is_refseq")
                .cast(pl.Boolean)
                .fill_null(False)
                .alias("member_is_refseq"),
            )
            .with_columns(
                pl.when(pl.col("member_is_refseq"))
                .then(pl.col("ncbi_taxid"))
                .otherwise(None)
                .alias("refseq_ncbi_taxid"),
                pl.when(pl.col("member_is_refseq"))
                .then(pl.col("ictv_taxid"))
                .otherwise(None)
                .alias("refseq_ictv_taxid"),
            )
            .group_by("fasta_accession")
            .agg(
                pl.col("representative_accession")
                .unique()
                .alias("representative_accessions"),
                pl.col("ncbi_taxid").unique().alias("ncbi_taxids"),
                pl.col("refseq_ncbi_taxid")
                .drop_nulls()
                .unique()
                .alias("refseq_ncbi_taxids"),
                pl.col("ictv_taxid").drop_nulls().unique().alias("ictv_taxids"),
                pl.col("refseq_ictv_taxid")
                .drop_nulls()
                .unique()
                .alias("refseq_ictv_taxids"),
                pl.len().alias("viral_member_count"),
                pl.col("member_is_refseq").sum().alias("refseq_viral_member_count"),
                pl.col("ictv_taxid").is_not_null().sum().alias("mapped_member_count"),
                (
                    pl.col("member_is_refseq")
                    & pl.col("ictv_taxid").is_not_null()
                )
                .sum()
                .alias("refseq_mapped_member_count"),
                pl.col("ictv_taxid").null_count().alias("unmapped_member_count"),
            )
            .sort("fasta_accession")
        )
        partials = {
            path: path.with_suffix(path.suffix + ".partial")
            for path in (
                assignments_path,
                mmseqs_mapping,
                diamond_mapping,
                audit_path,
            )
        }
        taxon_cache = {}

        def taxon_for(taxid: str):
            dense_taxid = str(taxid)
            if dense_taxid not in taxon_cache:
                taxon_cache[dense_taxid] = taxopy.Taxon(dense_taxid, taxdb)
            return taxon_cache[dense_taxid]

        missing = 0
        refseq_assigned = 0
        genus_rank_index = ICTV_RANKS.index("genus")
        with (
            partials[assignments_path].open("w", encoding="utf-8") as assignments,
            partials[mmseqs_mapping].open("w", encoding="utf-8") as mmseqs,
            partials[diamond_mapping].open("w", encoding="utf-8") as diamond,
            partials[audit_path].open("w", encoding="utf-8") as audit,
        ):
            assignments.write(
                "primary_name\tcanonical_accession\ttaxid\ttaxon_name\t"
                "represented_taxids\trepresented_taxid_count\tassignment_method\n"
            )
            diamond.write("accession\taccession.version\ttaxid\tgi\n")
            audit.write(
                "representative_accession\tfasta_accession\t"
                "ncbi_viral_member_taxids\t"
                "refseq_viral_member_taxids\t"
                "mapped_ictv_taxids\trefseq_mapped_ictv_taxids\t"
                "selected_ictv_taxids\tassigned_ictv_taxid\t"
                "viral_member_count\trefseq_viral_member_count\t"
                "mapped_member_count\trefseq_mapped_member_count\t"
                "unmapped_member_count\tassignment_scope\tassignment_method\n"
            )
            for (
                fasta_accession,
                representative_accessions,
                ncbi_taxids,
                refseq_ncbi_taxids,
                ictv_taxids,
                refseq_ictv_taxids,
                viral_member_count,
                refseq_viral_member_count,
                mapped_member_count,
                refseq_mapped_member_count,
                unmapped,
            ) in grouped.iter_rows():
                representative_accessions = sorted(
                    str(accession)
                    for accession in representative_accessions
                    if accession is not None and str(accession)
                )
                ncbi_taxids = sorted(
                    str(taxid)
                    for taxid in ncbi_taxids
                    if taxid is not None and str(taxid)
                )
                refseq_ncbi_taxids = sorted(
                    str(taxid)
                    for taxid in refseq_ncbi_taxids
                    if taxid is not None and str(taxid)
                )
                ictv_taxids = sorted(
                    str(taxid)
                    for taxid in ictv_taxids
                    if taxid is not None and str(taxid)
                )
                refseq_ictv_taxids = sorted(
                    str(taxid)
                    for taxid in refseq_ictv_taxids
                    if taxid is not None and str(taxid)
                )
                if refseq_ncbi_taxids:
                    source_taxids = refseq_ictv_taxids
                    assignment_scope = "refseq_viral_members"
                else:
                    source_taxids = ictv_taxids
                    assignment_scope = "all_viral_members"
                represented_taxids = sorted(
                    {
                        taxid_remap.get(str(taxid), str(taxid))
                        for taxid in source_taxids
                        if taxid is not None and str(taxid)
                    },
                    key=int,
                )
                if not represented_taxids:
                    missing += 1
                    if not allow_missing_taxonomy:
                        continue
                    assigned = taxon_for("1")
                    method = "unresolved_root"
                    assignment_scope = "unresolved"
                elif len(represented_taxids) == 1:
                    assigned = taxon_for(represented_taxids[0])
                    method = (
                        "refseq_direct"
                        if assignment_scope == "refseq_viral_members"
                        else "direct"
                    )
                    if assignment_scope == "refseq_viral_members":
                        refseq_assigned += 1
                else:
                    assigned = taxopy.find_lca(
                        [taxon_for(taxid) for taxid in represented_taxids],
                        taxdb,
                    )
                    method = (
                        "refseq_lca"
                        if assignment_scope == "refseq_viral_members"
                        else "lca"
                    )
                    if assignment_scope == "refseq_viral_members":
                        refseq_assigned += 1
                rank_index = (
                    ICTV_RANKS.index(assigned.rank)
                    if assigned.rank in ICTV_RANKS
                    else -1
                )
                if rank_index > genus_rank_index:
                    for lineage_taxid in assigned.taxid_lineage:
                        if taxdb.taxid2rank.get(lineage_taxid) == "genus":
                            assigned = taxopy.Taxon(lineage_taxid, taxdb)
                            method += "_genus_cap"
                            break
                assignments.write(
                    f"{fasta_accession}\t{fasta_accession}\t{assigned.taxid}\t"
                    f"{assigned.name}\t{';'.join(represented_taxids)}\t"
                    f"{len(represented_taxids)}\t{method}\n"
                )
                mmseqs.write(f"{fasta_accession}\t{assigned.taxid}\n")
                diamond.write(
                    f"{fasta_accession}\t{fasta_accession}\t{assigned.taxid}\t\n"
                )
                audit.write(
                    f"{';'.join(representative_accessions)}\t{fasta_accession}\t"
                    f"{';'.join(ncbi_taxids)}\t"
                    f"{';'.join(refseq_ncbi_taxids)}\t"
                    f"{';'.join(ictv_taxids)}\t"
                    f"{';'.join(refseq_ictv_taxids)}\t"
                    f"{';'.join(represented_taxids)}\t{assigned.taxid}\t"
                    f"{viral_member_count}\t{refseq_viral_member_count}\t"
                    f"{mapped_member_count}\t{refseq_mapped_member_count}\t"
                    f"{unmapped}\t{assignment_scope}\t{method}\n"
                )
        if missing and not allow_missing_taxonomy:
            for partial in partials.values():
                partial.unlink(missing_ok=True)
            raise RuntimeError(
                f"{missing:,} ClusteredNR representatives lacked an ICTV member "
                "mapping; rerun with --allow-missing-taxonomy to assign them to root"
            )
        for final_path, partial_path in partials.items():
            partial_path.replace(final_path)
        assigned_count = (
            grouped.height if allow_missing_taxonomy else grouped.height - missing
        )
        logger.info(
            f"Assigned taxonomy for {assigned_count:,} representatives "
            f"({refseq_assigned:,} using RefSeq viral member mappings; "
            f"{missing:,} unresolved)"
        )
    else:
        logger.info("Reusing completed primary-sequence taxonomy assignments")

    mmseqs_dir = output_root / "mmseqs"
    diamond_dir = output_root / "diamond"
    mmseqs_dir.mkdir(parents=True, exist_ok=True)
    diamond_dir.mkdir(parents=True, exist_ok=True)
    mmseqs_db = mmseqs_dir / "ncbi_virus"
    diamond_db = diamond_dir / "ncbi_virus"
    marker = output_root / "build/ncbi_virus_taxdb.complete.json"
    mmseqs_dbtype = Path(str(mmseqs_db) + ".dbtype")
    diamond_output = Path(str(diamond_db) + ".dmnd")
    backend_fasta = final_fasta
    if final_fasta.suffix == ".bgz":
        gzip_alias = temp_dir / f"{final_fasta.stem}.gz"
        if gzip_alias.exists() or gzip_alias.is_symlink():
            if not gzip_alias.is_symlink() or gzip_alias.resolve() != final_fasta:
                raise FileExistsError(
                    f"Cannot create gzip compatibility alias: {gzip_alias}"
                )
        else:
            gzip_alias.symlink_to(final_fasta)
        backend_fasta = gzip_alias

    if force or not outputs_are_current((mmseqs_dbtype,), (final_fasta,)):
        if not run_command_comp(
            base_cmd="mmseqs createdb",
            positional_args=[str(backend_fasta), str(mmseqs_db)],
            positional_args_location="start",
            params={"dbtype": 1},
            output_file=str(mmseqs_db),
            logger=logger,
        ):
            raise RuntimeError("MMseqs2 createdb failed")

    mmseqs_taxonomy = Path(str(mmseqs_db) + "_taxonomy")
    mmseqs_taxonomy_inputs = (
        mmseqs_mapping,
        taxonomy_dir / "nodes.dmp",
        taxonomy_dir / "names.dmp",
    )
    if force or not outputs_are_current((mmseqs_taxonomy,), mmseqs_taxonomy_inputs):
        if not run_command_comp(
            base_cmd="mmseqs createtaxdb",
            positional_args=[str(mmseqs_db), str(temp_dir / "mmseqs-taxonomy")],
            positional_args_location="start",
            params={
                "ncbi-tax-dump": str(taxonomy_dir),
                "tax-mapping-file": str(mmseqs_mapping),
                "threads": threads,
            },
            check_output=False,
            logger=logger,
        ):
            raise RuntimeError("MMseqs2 createtaxdb failed")

    diamond_inputs = (
        final_fasta,
        diamond_mapping,
        taxonomy_dir / "nodes.dmp",
        taxonomy_dir / "names.dmp",
    )
    if force or not outputs_are_current((diamond_output,), diamond_inputs):
        if not run_command_comp(
            base_cmd="diamond makedb",
            params={
                "in": str(backend_fasta),
                "db": str(diamond_db),
                "taxonmap": str(diamond_mapping),
                "taxonnodes": str(taxonomy_dir / "nodes.dmp"),
                "taxonnames": str(taxonomy_dir / "names.dmp"),
                "threads": threads,
                "no-parse-seqids": True,
            },
            output_file=str(diamond_output),
            logger=logger,
        ):
            raise RuntimeError("DIAMOND makedb failed")

    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(
        json.dumps(
            {
                "completed_at": datetime.datetime.now(
                    datetime.timezone.utc
                ).isoformat(),
                "fasta": str(final_fasta),
                "clustered_nr_db": str(blast_db) if blast_db else None,
                "clustered_nr_metadata_db": str(metadata_db) if metadata_db else None,
                "clustered_nr_taxonomy_db": str(taxonomy_db) if taxonomy_db else None,
                "clustered_nr_blastdb_metadata": str(source_metadata),
                "representative_viral_taxids": str(representative_viral_taxids),
                "representative_taxonomy": str(representative_mapping),
                "representatives": str(representatives),
                "extracted_representatives": str(extracted_representatives),
                "missing_representatives": str(missing_representatives),
                "ncbi_taxid_to_ictv": str(ncbi_taxid_to_ictv),
                "source_taxdump": str(source_taxdump),
                "taxonomy": str(taxonomy_dir),
                "mmseqs": str(mmseqs_db),
                "diamond": str(diamond_output),
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    logger.info("Finished the complete NCBI-virus protein taxonomy build step")
    return {
        "fasta": final_fasta,
        "representatives": representatives,
        "representative_viral_taxids": representative_viral_taxids,
        "representative_taxonomy": representative_mapping,
        "taxonomy": taxonomy_dir,
        "mmseqs": mmseqs_db,
        "diamond": diamond_output,
        "completion_manifest": marker,
    }


if __name__ == "__main__":
    build_data()
