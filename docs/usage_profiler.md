# SILOsim Profiler
Infers substitution probabilities and mobile genetic element (MGE) dynamics leveraging publicly available genomes.
```
silosim profiler -h
usage:  profiler [-h] --input_genome INPUT_GENOME [--output OUTPUT] --species {sau,sepi,cdiff,kp}
                 [--top_genomes TOP_GENOMES] [--cov_cutoff COV_CUTOFF] [--ani_threshold ANI_THRESHOLD]
                 [--bakta_db_type BAKTA_DB_TYPE] [--bakta_db_path BAKTA_DB_PATH]
                 [--whatsgnu_db_path WHATSGNU_DB_PATH] [-t THREADS] [--prefix PREFIX] [--conda_prefix CONDA_PREFIX]
                 [--force]

  -h, --help            show this help message and exit
  --input_genome INPUT_GENOME
                        Path to genome assembly in FASTA format
  --output OUTPUT       Output directory of the profiling results. If not provided, defaults to silosim_profiler_output_<YYYY_MM_DD_HHMMSS> under the current working directory.
  --species {sau,sepi,cdiff,kp}
                        Bacteria species.
                        Available options: [sau, sepi, cdiff, kp]
                        sau: Staphylococcus aureus
                        sepi: Staphylococcus epidermidis
                        cdiff: Clostridium difficile
                        kp: Klebsiella pneumoniae
  --local_query_dir LOCAL_QUERY_DIR
                        Path to the local query directory containing genome assemblies in FASTA format.If provided, the pipeline will use these local genomes for profiling instead of downloading them from NCBI.The WhatsGNU-based approach to fetch top genomes will be skipped if this option is provided.
  --top_genomes TOP_GENOMES
                        Number of initial top genomes for profiling before applying the ANI exclusion threshold (default: 1000).
  --min_ctg_len MIN_CTG_LEN
                        Minimum contig length in bp. Contigs shorter than this are excluded 
                        from the concatenated sequence used for profiling (default: 1000).
  --cov_cutoff COV_CUTOFF
                        Coverage cutoff for MGE inference.
                        Genomic regions with coverage above this cutoff will be inferred as MGEs.
                        Default is 0.7.
  --ani_threshold ANI_THRESHOLD
                        Average Nucleotide Identity (ANI) exclusion threshold (default: 95).
                        Genomes with ANI below this value will be removed from the profiling analysis.
  --bakta_db_type BAKTA_DB_TYPE
                        Bakta database.
                        Available options: [full, light]
                        Default is full
  --bakta_db_path BAKTA_DB_PATH
                        The path of the directory where the existing Bakta database locates.
                        If provided, the Bakta database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/bakta/db
  --whatsgnu_db_path WHATSGNU_DB_PATH
                        The path to the existing WhatsGNU database.
                        If provided, the WhatsGNU database will not be downloaded.
                        If not provided, defaults to <OUTPUT>/whatsgnu/db.
  -t THREADS, --threads THREADS
                        Thread number. Default is 1.
  --prefix PREFIX       Prefix for config file, output files, and analysis naming. If not provided, defaults to timestamp: YYYY_MM_DD_HHMMSS
  --conda_prefix CONDA_PREFIX
                        Directory for conda environments needed for this analysis. If not provided, defaults to OUTPUT/conda_envs_<YYYY_MM_DD_HHMMSS>
  --force               Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline.
                        This may cause instability or failures.
```

## Required Input
1. **Input Genome(--input_genome):**
    
    Path to the genome assembly in FASTA format for profiling.

2. **Species(--species):**

    Indicate the species being analyzed. Currently supported species include:
      - `sau` (*Staphylococcus aureus*)
      - `sepi` (*Staphylococcus epidermidis*)
      - `cdiff` (*Clostridium difficile*)
      - `kp` (*Klebsiella pneumoniae*)

## Optional Input
1. **Output Directory(--output):**
    - Path to the output directory. If not provided, defaults to `silosim_profiler_output_<YYYY_MM_DD_HHMMSS>` under the current working directory.
2. **Local Query Directory(--local_query_dir):**
    - Path to the local query directory containing genome assemblies in FASTA format.
    - If provided, the pipeline will use these local genomes for profiling instead of downloading them from NCBI.
    - The WhatsGNU-based approach to fetch top genomes will be skipped if this option is provided.
3. **Top Genomes(--top_genomes):**
    - Number of initial top genomes for profiling before applying the ANI exclusion threshold (default: 1000).
4. **ANI Exclusion Threshold(--ani_threshold):**
    - Average Nucleotide Identity (ANI) exclusion threshold (default: 95).
    - Genomes with ANI below this value will be removed from the profiling analysis.
5. **Bakta Database Type(--bakta_db_type) and Path(--bakta_db_path):**
    - Specify the type of Bakta database (`full` or `light`, default is `full`).
    - Path to an existing Bakta database. If not provided, the database will be downloaded to `<OUTPUT>/bakta/db`.
6. **WhatsGNU Database Path(--whatsgnu_db_path):**

    Path to an existing WhatsGNU database. If not provided, the database will be downloaded to `<OUTPUT>/whatsgnu/db`.
7. **Minimum Contig Length(--min_ctg_len):**

    Minimum contig length in bp. Contigs shorter than this are excluded from the concatenated sequence used for profiling (default: 1000).
8. **Coverage Cutoff(--cov_cutoff):**

    Coverage cutoff for MGE inference. Genomic regions with coverage above this cutoff will be inferred as MGEs.
    Default is 0.7, which means genomic regions with coverage above 70% will be inferred as MGEs.
9. **Threads(--threads / -t):**

    Number of threads to use (default is 1). It is highly recommended to increase the thread count using this option to improve performance, regardless of dataset size.
10. **Prefix(--prefix):**

    Prefix for config file, output files, and analysis naming. If not provided, defaults to a timestamp in the format `YYYY_MM_DD_HHMMSS`.
11. **Conda Prefix(--conda_prefix):**

    Directory for conda environments needed for this analysis. If not provided, defaults to `<OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>`.
12. **Force Execution(--force):**

    Bypass system compatibility checks (operating system and available RAM) and force execution of the pipeline. This may cause instability or failures.  

## Output
1. **Bakta Annotation:**
  - `bakta_annotation/{genome_name}/`
2. **WhatsGNU:**
  - `whatsgnu/{genome_name}/`
3. **FastANI:**
 - Path to all query genomes used in FastANI:
 
    `fastani/{genome_name}_topgenomes_all.txt`
 - FastANI results between the input genome and all query genomes by top_genomes option:
 
    `fastani/{genome_name}_fastani.csv `
 - Path to filtered query genomes after applying ANI exclusion threshold:

    `fastani/{genome_name}_topgenomes_filtered.txt`

4. **Pair-wise comparison using Mummer4:**
 - Scripts for running Mummer4 pair-wise comparison: `mummer4/scripts/`
 - Raw Mummer4 Output: `mummer4/output/`

5. **Public genomes downloaded using Datasets:**
 - `datasets_topgenomes/*.fna`

6. **SILOSim Profiler Results:**
 - Bin summary of concatenated input genome:

    `profiler/{genome_name}_bin_summary.csv`
    - Columns included in the bin summary csv file:
      - bin_index: Index of the bin
      - start: The starting position of the bin in the concatenated genome
      - end: The ending position of the bin in the concatenated genome
      - length: The length of the bin in base pairs
      - gene: The annotation of the bin by Bakta
      - locus_tag: The locus tag of the bin by Bakta
      - original_contig: The original contig name where the bin is located
      - coverage_max: The maximum coverage of the bin across all query genomes
      - coverage_min: The minimum coverage of the bin across all query genomes
      - coverage_mean: The mean coverage of the bin across all query genomes
      - coverage_sd: The standard deviation of the coverage of the bin across all query genomes
      - coverage_median: The median coverage of the bin across all query genomes
      - coverage_q1: The first quartile (25th percentile) of the coverage of the bin across all query genomes
      - coverage_q3: The third quartile (75th percentile) of the coverage of the bin across all query genomes
      - entropy_max: The maximum substitution probability of the bin across all query genomes
      - entropy_min: The minimum substitution probability of the bin across all query genomes
      - entropy_mean: The mean substitution probability of the bin across all query genomes
      - entropy_sd: The standard deviation of the substitution probability of the bin across all query genomes
      - entropy_median: The median substitution probability of the bin across all query genomes
      - entropy_q1: The first quartile (25th percentile) of the substitution probability
      - entropy_q3: The third quartile (75th percentile) of the substitution probability
      - total_snp_count: The total number of SNPs in the bin across all query genomes
      - mge_id: Whether the bin is identified as a mobile genetic element (MGE) or not. MGE ID is provided if identified as MGE.

 - Summary of MGEs identified in the input genome:

    `profiler/{genome_name}_mges.csv`
    - Columns included in the MGE summary csv file:
      - mge_index: Unique identifier for each MGE
      - start: The starting position of the MGE in the concatenated genome
      - end: The ending position of the MGE in the concatenated genome
      - length: The length of the MGE in base pairs
      - bin: The bin indexes that the MGE spans
      - bin_gene: The genes annotated in the bins that the MGE spans

 - Plots indicating the MGEs and substitution probabilities across the genome:

    `profiler/{genome_name}_profiler_plot.pdf`

    `profiler/{genome_name}_profiler_plot.RDS` (R object)
    - Figure description:
        - Shared X-axis: Position in the concatenated genome (in base pairs)
        - Upper panel: Coverage across the genome based on the query genomes after applying the ANI exclusion threshold.
            - Y-axis: Coverage Percentage (Percentage of query genomes covering the position in the concatenated genome)
            - Blue blocks: Identified MGEs in the concatenated genome
        - Lower panel: Substitution probabilities across the genome based on the query genomes after applying the ANI exclusion threshold.
            - Y-axis: Substitution Probability (Entropy)
            - Colored line: Substitution probabilities at each position across the genome
                - Red: Intergenic regions
                - Green: Coding sequences (CDS)
                
- Substitution probability table across each site within the genome:

    `profiler/{genome_name}_entropy.csv`

    `profiler/{genome_name}_entropy.RDS` (R object)