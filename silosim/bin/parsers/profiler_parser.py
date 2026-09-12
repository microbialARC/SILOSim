"""Parser for profiler command"""
import argparse
import os

# The parser for the overall profiler command
def add_profiler_parser(subparsers):
    profiler_parser = subparsers.add_parser(
        "profiler",
        formatter_class=argparse.RawTextHelpFormatter,
        help="Infers the probability of substitutions and mobile genetic elements."
    )

    profiler_parser.add_argument(
        "--input_genome",
        required=True,
        help="Path to genome assembly in FASTA format"
    )

    profiler_parser.add_argument(
        "--output",
        required=False,
        help="Output directory of the profiling results. If not provided, defaults to silosim_profiler_output_<YYYY_MM_DD_HHMMSS> under the current working directory."
    )
    
    profiler_parser.add_argument(
        "--species",
        type=str,
        required=False,
        default=None,
        choices=["sau", "sepi", "cdiff", "kp", "bp"],
        help="""Bacteria species.
Available options: [sau, sepi, cdiff, kp, bp]
sau: Staphylococcus aureus
sepi: Staphylococcus epidermidis
cdiff: Clostridium difficile
kp: Klebsiella pneumoniae
bp: Bordetella pertussis"""
    )

    profiler_parser.add_argument(
        "--local_query_dir",
        type=str,
        required=False,
        help="""Path to the local query directory containing genome assemblies in FASTA format.
If provided, the WhatsGNU-based approach to fetch top genomes will be skipped.
The pipeline will use these local genomes for profiling."""
    )

    profiler_parser.add_argument(
        "--top_genomes",
        type=int,
        default=1000,
        help="Number of initial top genomes for profiling before applying the ANI exclusion threshold (default: 1000)."
    )

    profiler_parser.add_argument(
        "--min_ctg_len",
        type=int,
        required=False,
        default=1000,
        help="""Minimum contig length in bp. Contigs shorter than this are excluded 
from the concatenated sequence used for profiling (default: 1000)."""
    )

    profiler_parser.add_argument(
        "--cov_cutoff",
        type=float,
        default=0.7,
        help="""Coverage cutoff for MGE inference.
Genomic regions with coverage above this cutoff will be inferred as MGEs.
Default is 0.7."""
    )

    profiler_parser.add_argument(
        "--ani_threshold",
        type=float,
        default=95,
        help="""Average Nucleotide Identity (ANI) exclusion threshold (default: 95).
Genomes with ANI below this value will be removed from the profiling analysis."""
    )

    profiler_parser.add_argument(
        "--bakta_db_type",
        required=False,
        default="full",
        help="""Bakta database.
Available options: [full, light]
Default is full"""
    )

    profiler_parser.add_argument(
        "--bakta_db_path",
        required=False,
        type=str,
        help="""The path of the directory where the existing Bakta database locates.
If provided, the Bakta database will not be downloaded.
If not provided, defaults to <OUTPUT>/bakta/db""",
    )

    profiler_parser.add_argument(
        "--whatsgnu_db_path",
        type=str,
        required=False,
        help="""The path to the existing WhatsGNU database.
If provided, the WhatsGNU database will not be downloaded.
If not provided, defaults to <OUTPUT>/whatsgnu/db.""",
    )
    profiler_parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help = "Thread number. Default is 1."
    )

    profiler_parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help="Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS"
    )

    profiler_parser.add_argument(
        "--conda_prefix",
        default=None,
        help="Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>"
    )

    # Override system compatibility checks (OS and minimum RAM) and run pipeline regardless
    profiler_parser.add_argument(
        "--force",
        action="store_true",
        help="""Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
This may cause instability or failures."""
    )
