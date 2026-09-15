# This script is invoked by the `profiler` rule in profiler.smk through
# Snakemake's `script:` directive, so every parameter arrives on the
# `snakemake` object. Nothing is read from the command line.
#
# To run it standalone instead, replace the block below with your own parsing
# of commandArgs(trailingOnly = TRUE) and supply the same eight values. A
# reasonable flag set would be:
#
#   Rscript profiler.R \
#     --input_genome <ref.fna>   reference assembly (FASTA; .fasta/.fna)  -> input_genome_path
#     --gff <ref.gff3>           annotation (GFF/GFF3)                    -> gff_path
#     --snp <dir>                directory of MUMmer *.snps and *.coords  -> snp_dir
#     --output <dir>             directory for all outputs                -> output_dir
#     --min_ctg_len <int>        drop contigs shorter than this, in bp    -> min_ctg_len
#     --cov_cutoff <float>       MGE threshold, fraction of max coverage  -> cov_cutoff
#     --min_mge_size <int>       shortest run of low coverage to call     -> min_mge_size
#     --cpus <int>               cores to use                             -> ncores
#
# profiler.R -- reference-genome entropy and mobile-element profiling
# Part of the SILOSim pipeline. Vectorised + parallel implementation.
#
# Outputs (all written to output_dir)
#   <name>_concat.fasta                   concatenated reference
#   <name>_new_pos.RDS                    contig -> concatenated coordinate map
#   <name>_position_coverage.RDS          per-site query-genome coverage
#   <name>_snps_sum.RDS                   all SNP calls
#   <name>_chr_bins.RDS                   CDS / intergenic bins
#   <name>_entropy.RDS / .csv             per-site entropy (CSV feeds the simulator)
#   <name>_mges.RDS / .csv                inferred MGE intervals
#   <name>_mges_seq.fasta                 MGE sequences
#   mge_entropy/<name>_MGE_N_entropy.csv  per-MGE entropy
#   <name>_bin_summary.csv                per-bin coverage / entropy / SNP summary
#   <name>_profiler_plot.pdf              coverage + entropy figure
#   <name>_plot_df.RDS                    the two data frames behind the figure


# Libraries ----
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(parallel))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(entropy))
suppressMessages(library(cowplot))
suppressMessages(library(ggnewscale))
suppressMessages(library(ggExtra))
suppressMessages(library(IRanges))

# Input from arguments ----
# Command-line Rscript tool. Accepts these flags when run standalone:
#   --input_genome <path>   : reference assembly (FASTA; .fasta/.fna). Required.
#   --min_ctg_len <int>       : minimum contig length in bp. Contigs shorter than this are excluded from the concatenated sequence used for profiling (optional; defaults to 1000).
#   --gff <path>            : annotation file (GFF/GFF3). Required.
#   --snp <dir>             : directory with SNP/coord outputs (expects *.snps and *.coords from MUMmer). Required.
#   --cpus <int>            : number of CPU cores to use (optional; defaults to detectCores()).
#   --output <dir>          : directory where all outputs (RDS/CSV/PDF) will be written. Required.
#
# Snakemake integration:
#   When run inside a Snakemake rule, Snakemake should invoke Rscript with the same flags
#   (see the pipeline's profiler.smk rule). The script reads commandArgs(trailingOnly=TRUE).
#
# Outputs (written to --output):
#   *_concat.fasta, *_new_pos.RDS, *_position_coverage.RDS, *_snps_sum.RDS, *_chr_bins.RDS,
#   *_entropy.RDS / .csv, *_mges.RDS / .csv / _mges_seq.fasta, *_bin_summary.csv, *_profiler_plot.pdf

# When run inside the Thresher pipeline, Snakemake calls this script and provides the necessary command-line arguments.
# Therefore the argument-parsing code (args and related blocks) is commented out.
# To run the script standalone, edit the file and uncomment those sections or supply the arguments manually.
#args <- commandArgs(trailingOnly = TRUE)
#for (i in seq_along(args)) {
#  if (args[i] == "--input_genome") {
#    input_genome_path <- args[i + 1]
#  } else if (args[i] == "--min_ctg_len") {
#    min_ctg_len <- as.integer(args[i + 1])
#  } else if (args[i] == "--gff") {
#    gff_path <- args[i + 1]
#  }else if (args[i] == "--snp") {
#    snp_dir <- args[i + 1]
#  }else if (args[i] == "--cpus") {
#    ncores <- as.integer(args[i + 1])
#  }else if (args[i] == "--output") {
#    output_dir <- args[i + 1]
#  }
#}

# Directly get the input from snakemake when run inside the pipeline
# Path to the input genome assembly file
input_genome_path <- snakemake@input[["fna_path"]]
# Minimum contig length for concatenation
min_ctg_len <- as.integer(snakemake@params[["min_ctg_len"]])
# Minimum length (bp) of a genomic region to be considered an MGE
min_mge_size <- as.integer(snakemake@params[["min_mge_len"]])
# Path to the gff file of the input genome assembly
gff_path <- snakemake@input[["gff3_path"]]
# Path to directory where the output files will be saved
output_dir <- snakemake@params[["output_dir"]]
# Coverage cutoff for MGE inference
cov_cutoff <- as.numeric(snakemake@params[["cov_cutoff"]])
# Make directory if not exists
if(!dir.exists(output_dir)){
  dir.create(output_dir,recursive = TRUE)
}
# Path to the directory of corresponding snp results 
snp_dir <- snakemake@params[["snp_dir"]]
# how many cores will be used
ncores <- snakemake@threads

# Check if required arguments were provided

## Validate before doing any work ----
stop_if_missing <- function(argument_value, parameter_name) {
  if (is.null(argument_value) || length(argument_value) == 0 ||
      is.na(argument_value) || !nzchar(as.character(argument_value))) {
    stop("Missing required parameter: ", parameter_name, call. = FALSE)
  }
}

stop_if_missing(input_genome_path, "input_genome")
stop_if_missing(gff_path, "gff")
stop_if_missing(snp_dir, "snp")
stop_if_missing(output_dir, "output")
stop_if_missing(min_ctg_len, "min_ctg_len")
stop_if_missing(cov_cutoff, "cov_cutoff")

if (!file.exists(input_genome_path)) stop("Input genome not found: ", input_genome_path)
if (!file.exists(gff_path)) stop("GFF file not found: ", gff_path)
if (!dir.exists(snp_dir)) stop("SNP directory not found: ", snp_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)


message("")
message("Input genome: ", input_genome_path)
message("Minimum contig length: ", min_ctg_len)
message("Coverage cutoff: ", cov_cutoff)
message("Minimum MGE length: ", min_mge_size)
message("GFF file: ", gff_path)
message("SNP directory: ", snp_dir)
message("Output directory: ", output_dir)
message("Threads: ", ncores)
message("")


# Shared helpers ----

## Run a list of independent tasks across cores ----

run_tasks_in_parallel <- function(task_list, worker_function, n_cores) {
  
  if (length(task_list) == 0) return(list())
  if (n_cores <= 1L) return(lapply(task_list, worker_function))
  
  task_results <- parallel::mclapply(
    task_list,
    worker_function,
    mc.cores = min(n_cores, length(task_list)),
    mc.preschedule = FALSE
  )
  
  worker_failed <- vapply(task_results, inherits, logical(1), "try-error")
  if (any(worker_failed)) {
    stop("Parallel worker failed on ", sum(worker_failed), " task(s). First error: ",
         as.character(task_results[[which(worker_failed)[1]]]), call. = FALSE)
  }
  
  task_results
}

## Split 1..n_items into at most n_chunks contiguous blocks ----
# Used when the per-item work is tiny (single bins, single entropy values):
# forking once per item would cost more than the work itself.
split_into_chunks <- function(n_items, n_chunks) {
  if (n_items == 0) return(list())
  n_chunks <- max(1L, min(as.integer(n_chunks), n_items))
  split(seq_len(n_items), cut(seq_len(n_items), breaks = n_chunks, labels = FALSE))
}

## Turn a (position, value) table into a plain vector indexed by position ----
# Every table in this script covers positions 1..genome_length exactly once, so
# a dense vector is both smaller and far faster to query than repeated
# `table$value[table$position >= a & table$position <= b]` scans.
scatter_to_position_vector <- function(positions, values, genome_length, fill = 0) {
  position_vector <- rep(fill, genome_length)
  position_vector[positions] <- values
  position_vector
}

## Pull one attribute out of a GFF3 column 9 field ----
# Vectorised over the whole column. The leading ";" lets a fixed-width
# lookbehind anchor the key, so `product=Fake Name=x` cannot masquerade as a
# real `Name=` attribute.
extract_gff_attribute <- function(attribute_field, attribute_key, default = NA_character_) {
  
  padded_field <- paste0(";", attribute_field)
  attribute_pattern <- paste0("(?<=;)", attribute_key, "=[^;]*")
  
  match_position <- regexpr(attribute_pattern, padded_field, perl = TRUE)
  attribute_values <- rep(default, length(attribute_field))
  
  key_was_found <- match_position > 0
  attribute_values[key_was_found] <- sub(
    paste0("^", attribute_key, "="), "",
    regmatches(padded_field, match_position)
  )
  
  attribute_values
}

## Read the CDS records out of a GFF/GFF3 file ----
# Bakta the assembly itself after a `##FASTA` marker. Feeding
# those lines to a TSV reader produces junk rows, so cut the file there first.
read_gff_cds <- function(gff_path) {
  
  gff_lines <- readLines(gff_path, warn = FALSE)
  
  fasta_marker_line <- which(startsWith(gff_lines, "##FASTA"))
  if (length(fasta_marker_line) > 0) {
    gff_lines <- gff_lines[seq_len(fasta_marker_line[1] - 1L)]
  }
  gff_lines <- gff_lines[nzchar(gff_lines) & !startsWith(gff_lines, "#")]
  
  if (length(gff_lines) == 0) stop("No annotation records found in ", gff_path)
  
  # quote = "" because GFF3 attribute values legitimately contain quote marks.
  gff_table <- data.table::fread(text = gff_lines, sep = "\t", header = FALSE,
                                 quote = "", showProgress = FALSE)
  gff_table <- gff_table[gff_table$V3 == "CDS", ]
  
  if (nrow(gff_table) == 0) {
    warning("No CDS features in ", gff_path, "; every bin will be intergenic.")
    return(data.frame(contig = character(0), start = integer(0), end = integer(0),
                      gene_name = character(0), locus_tag = character(0),
                      stringsAsFactors = FALSE))
  }
  
  locus_tags <- extract_gff_attribute(gff_table$V9, "locus_tag", default = NA_character_)
  gene_names <- extract_gff_attribute(gff_table$V9, "Name", default = NA_character_)
  
  # Many Bakta CDS carry no Name=. Fall back to the locus tag rather than
  # dropping the bin, which is what the original did when the grep came back
  # empty.
  gene_names[is.na(gene_names)] <- locus_tags[is.na(gene_names)]
  gene_names[is.na(gene_names)] <- "unnamed_cds"
  locus_tags[is.na(locus_tags)] <- "unnamed_cds"
  
  data.frame(
    contig = as.character(gff_table$V1),
    start = as.integer(gff_table$V4),
    end = as.integer(gff_table$V5),
    gene_name = gene_names,
    locus_tag = locus_tags,
    stringsAsFactors = FALSE
  )
}

# Step 1: concatenate the reference genome ----
#
# Draft assemblies come as many contigs. Everything downstream works in a
# single "concatenated" coordinate system, so we glue the surviving contigs
# together and keep a map back to the original contig coordinates.

concat_genome <- function(input_genome_name,
                          input_genome_path,
                          min_ctg_len,
                          output_dir) {
  
  input_genome_fasta <- Biostrings::readDNAStringSet(input_genome_path)
  
  if (max(width(input_genome_fasta)) < min_ctg_len) {
    stop("The longest contig (", max(width(input_genome_fasta)),
         " bp) is shorter than min_ctg_len (", min_ctg_len, " bp).")
  }
  
  contig_passes_filter <- width(input_genome_fasta) >= min_ctg_len
  
  message("Contigs kept: ", sum(contig_passes_filter),
          " (", sum(width(input_genome_fasta)[contig_passes_filter]), " bp)")
  message("Contigs filtered: ", sum(!contig_passes_filter),
          " (", sum(width(input_genome_fasta)[!contig_passes_filter]), " bp)")
  
  input_genome_fasta <- input_genome_fasta[contig_passes_filter]
  
  # nucmer keeps only the first whitespace-delimited word of a FASTA header, so
  # contig names must be trimmed the same way here or the joins below miss.
  names(input_genome_fasta) <- sub("\\s.*$", "", names(input_genome_fasta))
  
  contig_lengths <- as.integer(width(input_genome_fasta))
  contig_ends <- cumsum(contig_lengths)
  
  ctg_pos_df <- data.frame(
    ctg = names(input_genome_fasta),
    length = contig_lengths,
    new_start = c(1L, head(contig_ends, -1) + 1L),
    new_end = contig_ends,
    stringsAsFactors = FALSE
  )
  
  if (length(input_genome_fasta) > 1) {
    input_genome_fasta_concat <- Biostrings::DNAStringSet(unlist(input_genome_fasta))
    names(input_genome_fasta_concat) <- input_genome_name
  } else {
    input_genome_fasta_concat <- input_genome_fasta
  }
  
  Biostrings::writeXStringSet(
    input_genome_fasta_concat,
    file.path(output_dir, paste0(input_genome_name, "_concat.fasta"))
  )
  
  new_pos <- list(
    concat_size = as.integer(width(input_genome_fasta_concat)),
    ctg_pos_df = ctg_pos_df
  )
  
  saveRDS(new_pos, file.path(output_dir, paste0(input_genome_name, "_new_pos.RDS")))
  
  new_pos
}

# Step 2: per-site query-genome coverage ----
#
# "Coverage" here means: how many DISTINCT query genomes have an alignment
# spanning this reference position. Not how many alignment blocks -- one query
# can align to the same region in several fragments and must still count once.
#
# The original asked that question one position at a time, scanning the whole
# alignment table for each. This version answers it for the entire contig in
# one pass:
#   1. merge each query's own overlapping intervals (IRanges::reduce on a list
#      grouped by query) so every query contributes at most +1 per site;
#   2. stack all merged intervals and let IRanges::coverage() count them.
# Result is identical, cost drops from O(genome x alignments) to O(n log n).

get_position_coverage <- function(input_genome_name,
                                  new_pos,
                                  snp_dir,
                                  output_dir,
                                  ncores) {
  
  contig_positions <- new_pos$ctg_pos_df
  
  coords_files <- list.files(snp_dir, pattern = "\\.coords$", full.names = TRUE)
  if (length(coords_files) == 0) {
    stop("No *.coords files found in ", snp_dir)
  }
  
  ## Read every .coords file; one query genome per file ----
  # Columns used (MUMmer show-coords): V1 = reference start, V2 = reference end,
  # V3 = reference contig tag. The file index is the query genome identity.
  read_one_coords_file <- function(file_idx) {
    
    data.table::setDTthreads(1L)
    
    coords_path <- coords_files[file_idx]
    if (file.size(coords_path) == 0) return(NULL)
    
    coords_table <- data.table::fread(coords_path, header = FALSE,
                                      select = 1:3, showProgress = FALSE)
    if (nrow(coords_table) == 0) return(NULL)
    
    data.table::data.table(
      ref_start = pmin(as.integer(coords_table$V1), as.integer(coords_table$V2)),
      ref_end = pmax(as.integer(coords_table$V1), as.integer(coords_table$V2)),
      ref_contig = as.character(coords_table$V3),
      query_id = file_idx
    )
  }
  
  all_alignments <- data.table::rbindlist(
    run_tasks_in_parallel(as.list(seq_along(coords_files)), read_one_coords_file, ncores)
  )
  
  message("Loaded ", format(nrow(all_alignments), big.mark = ","),
          " alignment blocks from ", length(coords_files), " query genome(s).")
  
  ## Count distinct covering query genomes, one contig at a time ----
  coverage_for_one_contig <- function(contig_row) {
    
    contig_name <- contig_positions$ctg[contig_row]
    contig_length <- contig_positions$length[contig_row]
    contig_offset <- contig_positions$new_start[contig_row] - 1L
    
    contig_alignments <- all_alignments[all_alignments$ref_contig == contig_name, ]
    
    if (nrow(contig_alignments) == 0) {
      
      coverage_depth <- integer(contig_length)
      
    } else {
      
      # Clip alignments to the contig
      # drop anything left entirely outside.
      
      clipped_start <- pmax(contig_alignments$ref_start, 1L)
      clipped_end <- pmin(contig_alignments$ref_end, contig_length)
      row_is_usable <- clipped_start <= clipped_end
      
      if (!any(row_is_usable)) {
        coverage_depth <- integer(contig_length)
      } else {
        
        # Merge each query's blocks with its own before counting, so a query
        # aligning in three fragments still contributes +1 per covered base.
        # Sorted by (query, start), a block opens a new merged interval only
        # when it starts past every block seen so far for that query.
        query_blocks <- data.table::data.table(
          query_id = contig_alignments$query_id[row_is_usable],
          block_start = clipped_start[row_is_usable],
          block_end = clipped_end[row_is_usable]
        )
        
        data.table::setorder(query_blocks, query_id, block_start)
        
        query_blocks[, highest_end_so_far := data.table::shift(cummax(block_end), fill = 0L),
                     by = query_id]
        
        query_blocks[, interval_id := cumsum(block_start > highest_end_so_far),
                     by = query_id]
        
        merged_intervals <- query_blocks[, .(interval_start = block_start[1],
                                             interval_end = max(block_end)),
                                         by = .(query_id, interval_id)]
        
        coverage_depth <- as.integer(IRanges::coverage(
          IRanges::IRanges(start = merged_intervals$interval_start,
                           end = merged_intervals$interval_end),
          width = contig_length
        ))
      }
    }
    
    data.table::data.table(
      position = contig_offset + seq_len(contig_length),
      coverage = coverage_depth
    )
  }
  
  position_coverage <- data.table::rbindlist(
    run_tasks_in_parallel(as.list(seq_len(nrow(contig_positions))),
                          coverage_for_one_contig, ncores)
  )
  
  saveRDS(position_coverage,
          file.path(output_dir, paste0(input_genome_name, "_position_coverage.RDS")))
  
  position_coverage
}


# Step 3: collect every SNP call ----
#
# One row per (SNP position, reference base, query base) per query genome.
# Duplicates are deliberate: the number of query genomes carrying a given
# alternate base is exactly what the entropy calculation counts.
#
# The original built a one-row data.frame per SNP and rbind-ed them. Here the
# contig-offset lookup and the indel filter are applied to the whole file at
# once.

get_snps_sum <- function(input_genome_name,
                         new_pos,
                         snp_dir,
                         output_dir,
                         ncores) {
  
  contig_positions <- new_pos$ctg_pos_df
  
  # Contig name -> how much to add to turn a contig-local position into a
  # concatenated-genome position. Contigs dropped by min_ctg_len are absent,
  # so their SNPs fall out as NA and are filtered.
  contig_offset_by_name <- setNames(contig_positions$new_start - 1L,
                                    contig_positions$ctg)
  
  snps_files <- list.files(snp_dir, pattern = "\\.snps$", full.names = TRUE)
  if (length(snps_files) == 0) {
    stop("No *.snps files found in ", snp_dir)
  }
  
  ## Read one show-snps file ----
  # Columns used: V1 = position in the reference, V2 = reference base,
  # V3 = query base, V5 = reference contig tag.
  read_one_snps_file <- function(snps_path) {
    
    data.table::setDTthreads(1L)
    
    if (!file.exists(snps_path) || file.size(snps_path) == 0) return(NULL)
    
    snps_table <- data.table::fread(snps_path, header = FALSE, showProgress = FALSE)
    if (nrow(snps_table) == 0 || ncol(snps_table) < 5) return(NULL)
    
    contig_offsets <- contig_offset_by_name[as.character(snps_table$V5)]
    
    # Insertions and deletions are marked "." on one side and are out of scope
    # for a per-site substitution entropy.
    row_is_usable <- !is.na(contig_offsets) &
      snps_table$V2 != "." &
      snps_table$V3 != "."
    
    if (!any(row_is_usable)) return(NULL)
    
    data.table::data.table(
      snp_pos = as.integer(contig_offsets[row_is_usable]) + as.integer(snps_table$V1[row_is_usable]),
      ref_site = as.character(snps_table$V2[row_is_usable]),
      new_site = as.character(snps_table$V3[row_is_usable])
    )
  }
  
  snps_sum <- data.table::rbindlist(
    run_tasks_in_parallel(as.list(snps_files), read_one_snps_file, ncores)
  )
  
  if (nrow(snps_sum) == 0) {
    warning("No usable SNP calls found in ", snp_dir, "; entropy will be zero everywhere.")
    snps_sum <- data.table::data.table(snp_pos = integer(0),
                                       ref_site = character(0),
                                       new_site = character(0))
  }
  
  message("Collected ", format(nrow(snps_sum), big.mark = ","), " SNP calls at ",
          format(data.table::uniqueN(snps_sum$snp_pos), big.mark = ","), " position(s).")
  
  saveRDS(snps_sum, file.path(output_dir, paste0(input_genome_name, "_snps_sum.RDS")))
  
  snps_sum
}

# Step 4: bin the genome ----
#
# Two kinds of bin, covering the concatenated genome without gaps:
#   - one bin per annotated CDS, spanning that CDS;
#   - one bin per stretch of a contig that no CDS covers ("non-cds").
#
# The original walked every base pair and asked "which CDS am I inside?", which
# is O(genome x genes). The same partition is just the CDS intervals plus the
# complement of their union within each contig, which IRanges computes directly.

get_chr_bins <- function(input_genome_name,
                         new_pos,
                         gff_path,
                         output_dir) {
  
  contig_positions <- new_pos$ctg_pos_df
  
  cds_table <- read_gff_cds(gff_path)
  
  ## Move CDS coordinates into concatenated space ----
  contig_offset_by_name <- setNames(contig_positions$new_start - 1L,
                                    contig_positions$ctg)
  cds_table$contig_offset <- contig_offset_by_name[cds_table$contig]
  
  cds_on_dropped_contig <- is.na(cds_table$contig_offset)
  if (any(cds_on_dropped_contig)) {
    message("Ignoring ", sum(cds_on_dropped_contig),
            " CDS on contigs removed by the min_ctg_len filter.")
    cds_table <- cds_table[!cds_on_dropped_contig, , drop = FALSE]
  }
  
  cds_table$new_start <- cds_table$contig_offset + cds_table$start
  cds_table$new_end <- cds_table$contig_offset + cds_table$end
  
  ## Coding bins: one per CDS ----
  coding_bins <- data.frame(
    start = cds_table$new_start,
    end = cds_table$new_end,
    gene = cds_table$gene_name,
    locus_tag = cds_table$locus_tag,
    original_contig = cds_table$contig,
    # Preserves GFF row order when two CDS share a start, matching the order in
    # which the original loop would have created them.
    source_order = seq_len(nrow(cds_table)),
    stringsAsFactors = FALSE
  )
  
  ## Intergenic bins: whatever the CDS leave uncovered, contig by contig ----
  cds_ranges <- if (nrow(cds_table) > 0) {
    IRanges::reduce(IRanges::IRanges(start = cds_table$new_start,
                                     end = cds_table$new_end))
  } else {
    IRanges::IRanges()
  }
  
  intergenic_bins <- data.table::rbindlist(lapply(seq_len(nrow(contig_positions)), function(contig_row) {
    
    contig_range <- IRanges::IRanges(start = contig_positions$new_start[contig_row],
                                     end = contig_positions$new_end[contig_row])
    uncovered_ranges <- IRanges::setdiff(contig_range, cds_ranges)
    
    if (length(uncovered_ranges) == 0) return(NULL)
    
    data.frame(
      start = IRanges::start(uncovered_ranges),
      end = IRanges::end(uncovered_ranges),
      gene = "non-cds",
      locus_tag = "non-cds",
      original_contig = contig_positions$ctg[contig_row],
      source_order = 0L,
      stringsAsFactors = FALSE
    )
  }))
  
  ## Assemble, order by position, number the bins ----
  chr_bins <- rbind(coding_bins, as.data.frame(intergenic_bins))
  chr_bins <- chr_bins[order(chr_bins$start, chr_bins$source_order), , drop = FALSE]
  
  chr_bins$bin_index <- seq_len(nrow(chr_bins))
  chr_bins$length <- chr_bins$end - chr_bins$start + 1L
  
  chr_bins <- chr_bins[, c("bin_index", "start", "end", "length",
                           "gene", "locus_tag", "original_contig")]
  rownames(chr_bins) <- NULL
  
  message("Built ", nrow(chr_bins), " bins (",
          sum(chr_bins$gene != "non-cds"), " coding, ",
          sum(chr_bins$gene == "non-cds"), " intergenic).")
  
  saveRDS(chr_bins, file.path(output_dir, paste0(input_genome_name, "_chr_bins.RDS")))
  
  chr_bins
}

# Step 5: per-site Shannon entropy ----
#
# At each position the allele counts across the query panel are:
#   reference allele : 1 (the reference itself) + coverage - total SNP calls
#   each alternate   : number of query genomes calling that base
#
# Entropy is the James-Stein shrinkage estimator from the `entropy` package
# (Hausser & Strimmer 2009, JMLR 10:1469). Positions with no SNP have entropy 0
# by definition, which is the overwhelming majority of the genome.
#
# Two changes carry the speed-up:
#   - only SNP-bearing positions are touched at all, via one data.table
#     group-by, instead of scanning the whole SNP table once per base pair;
#   - entropy depends only on the multiset of counts, so identical count
#     vectors are computed once and reused. On a real panel this collapses
#     hundreds of thousands of calls into a few thousand.

## Shrinkage entropy for one vector of allele counts ----
shrinkage_entropy <- function(observed_counts) {
  tryCatch(
    as.numeric(entropy::entropy(observed_counts, method = "shrink", verbose = FALSE)[1]),
    error = function(e) NA_real_
  )
}

get_entropy <- function(input_genome_name,
                        new_pos,
                        snps_sum,
                        position_coverage,
                        output_dir,
                        ncores) {
  
  data.table::setDTthreads(ncores)
  
  genome_length <- new_pos$concat_size
  coverage_vector <- scatter_to_position_vector(position_coverage$position,
                                                as.integer(position_coverage$coverage),
                                                genome_length, fill = 0L)
  
  entropy_vector <- numeric(genome_length)
  
  snp_table <- data.table::as.data.table(snps_sum)
  
  if (nrow(snp_table) > 0) {
    snp_out_of_range <- snp_table$snp_pos < 1L | snp_table$snp_pos > genome_length
    if (any(snp_out_of_range)) {
      warning("Dropping ", sum(snp_out_of_range),
              " SNP call(s) outside the concatenated genome.")
      snp_table <- snp_table[!snp_out_of_range, ]
    }
  }
  
  if (nrow(snp_table) > 0) {
    
    ## How many query genomes call each alternate base at each position ----
    allele_counts <- snp_table[, .(alt_count = .N), by = .(snp_pos, ref_site, new_site)]
    data.table::setorder(allele_counts, snp_pos, -alt_count)
    
    ## One row per SNP position ----
    position_summary <- allele_counts[, .(
      n_reference_bases = data.table::uniqueN(ref_site),
      total_alt_count = sum(alt_count),
      alt_count_key = paste(alt_count, collapse = ",")
    ), by = snp_pos]
    
    position_summary[, ref_count := 1L + coverage_vector[snp_pos] - total_alt_count]
    position_summary[, count_key := paste(ref_count, alt_count_key, sep = ",")]
    
    ## Ambiguous positions ----
    # More than one reference base at a single position means the alignments
    # disagree about the reference itself. The original silently produced NA
    # here (its tryCatch fired); we do the same but say so.
    position_is_ambiguous <- position_summary$n_reference_bases > 1L
    if (any(position_is_ambiguous)) {
      warning(sum(position_is_ambiguous),
              " position(s) report more than one reference base; entropy set to NA.")
    }
    
    ## Compute each distinct count vector once ----
    unique_count_keys <- unique(position_summary$count_key[!position_is_ambiguous])
    unique_count_vectors <- lapply(strsplit(unique_count_keys, ",", fixed = TRUE), as.numeric)
    
    message("Computing entropy at ",
            format(nrow(position_summary), big.mark = ","), " SNP position(s) using ",
            format(length(unique_count_keys), big.mark = ","), " distinct count vector(s).")
    
    entropy_chunks <- run_tasks_in_parallel(
      split_into_chunks(length(unique_count_keys), ncores * 4L),
      function(key_indices) {
        vapply(unique_count_vectors[key_indices], shrinkage_entropy, numeric(1))
      },
      ncores
    )
    entropy_by_key <- unlist(entropy_chunks, use.names = FALSE)
    
    ## Scatter the results back onto the genome ----
    entropy_vector[position_summary$snp_pos[!position_is_ambiguous]] <-
      entropy_by_key[match(position_summary$count_key[!position_is_ambiguous],
                           unique_count_keys)]
    
    entropy_vector[position_summary$snp_pos[position_is_ambiguous]] <- NA_real_
  }
  
  entropy_df <- data.table::data.table(
    position = seq_len(genome_length),
    entropy = entropy_vector
  )
  
  saveRDS(entropy_df, file.path(output_dir, paste0(input_genome_name, "_entropy.RDS")))
  
  # CSV feeds the evolution simulator.
  write.csv(entropy_df,
            file.path(output_dir, paste0(input_genome_name, "_entropy.csv")),
            quote = FALSE, row.names = FALSE)
  
  entropy_df
}

# Step 6: infer mobile genetic elements ----
#
# An MGE is called where coverage stays below `cov_cutoff` of the genome-wide
# maximum for at least `min_mge_size` consecutive bases: query genomes are
# systematically failing to align there, which is what accessory content looks
# like. Runs are found per contig so a call can never straddle the artificial
# junction between two concatenated contigs.

get_mges <- function(input_genome_name,
                     entropy_df,
                     new_pos,
                     position_coverage,
                     chr_bins,
                     min_mge_size = 100,
                     cov_cutoff,
                     output_dir,
                     ncores) {
  
  data.table::setDTthreads(ncores)
  
  contig_positions <- new_pos$ctg_pos_df
  genome_length <- new_pos$concat_size
  
  coverage_vector <- scatter_to_position_vector(position_coverage$position,
                                                as.integer(position_coverage$coverage),
                                                genome_length, fill = 0L)
  entropy_vector <- scatter_to_position_vector(entropy_df$position,
                                               entropy_df$entropy,
                                               genome_length, fill = 0)
  
  max_coverage <- max(coverage_vector)
  if (max_coverage == 0) {
    stop("No reference position is covered by any query genome; cannot infer MGEs.")
  }
  
  position_is_low_coverage <- (coverage_vector / max_coverage) < cov_cutoff
  
  ## Find runs of low coverage, contig by contig ----
  mge_intervals <- data.table::rbindlist(lapply(seq_len(nrow(contig_positions)), function(contig_row) {
    
    contig_start <- contig_positions$new_start[contig_row]
    contig_end <- contig_positions$new_end[contig_row]
    
    coverage_runs <- rle(position_is_low_coverage[contig_start:contig_end])
    
    run_end_local <- cumsum(coverage_runs$lengths)
    run_start_local <- c(1L, head(run_end_local, -1) + 1L)
    
    run_is_mge <- coverage_runs$values & coverage_runs$lengths >= min_mge_size
    if (!any(run_is_mge)) return(NULL)
    
    data.table::data.table(
      mge_contig = contig_positions$ctg[contig_row],
      mge_start = run_start_local[run_is_mge] + contig_start - 1L,
      mge_end = run_end_local[run_is_mge] + contig_start - 1L
    )
  }))
  
  ## No MGEs is a legitimate result, not an error ----
  if (nrow(mge_intervals) == 0) {
    
    message("No MGE met the size and coverage criteria.")
    
    mges_df <- data.frame(mge_index = character(0),
                          start = integer(0),
                          end = integer(0),
                          length = integer(0),
                          bin = character(0), bin_gene = character(0),
                          stringsAsFactors = FALSE)
    
    saveRDS(mges_df, file.path(output_dir, paste0(input_genome_name, "_mges.RDS")))
    write.csv(mges_df, file.path(output_dir, paste0(input_genome_name, "_mges.csv")),
              quote = FALSE, row.names = FALSE)
    writeLines(character(0),
               file.path(output_dir, paste0(input_genome_name, "_mges_seq.fasta")))
    
    return(mges_df)
  }
  
  ## Annotate each MGE with the bins and genes it spans ----
  mge_ranges <- IRanges::IRanges(start = mge_intervals$mge_start,
                                 end = mge_intervals$mge_end)
  bin_ranges <- IRanges::IRanges(start = chr_bins$start, end = chr_bins$end)
  
  mge_bin_overlaps <- IRanges::findOverlaps(mge_ranges, bin_ranges)
  overlapping_bins_by_mge <- split(S4Vectors::subjectHits(mge_bin_overlaps),
                                   S4Vectors::queryHits(mge_bin_overlaps))
  
  describe_bin_span <- function(bin_rows) {
    if (length(bin_rows) == 0) return(NA_character_)
    bin_indices <- chr_bins$bin_index[bin_rows]
    if (length(bin_indices) == 1) as.character(bin_indices)
    else paste0(min(bin_indices), "-", max(bin_indices))
  }
  
  describe_bin_genes <- function(bin_rows) {
    if (length(bin_rows) == 0) return("")
    paste(setdiff(chr_bins$gene[bin_rows], "non-cds"), collapse = "; ")
  }
  
  # split() only names the MGEs that hit a bin. Scattering into a pre-sized
  # list leaves the rest NULL, which describe_* read as "no overlap".
  bin_rows_per_mge <- vector("list", nrow(mge_intervals))
  bin_rows_per_mge[as.integer(names(overlapping_bins_by_mge))] <- overlapping_bins_by_mge
  
  mges_df <- data.frame(
    mge_index = paste0(input_genome_name, "_MGE_", seq_len(nrow(mge_intervals))),
    start = mge_intervals$mge_start,
    end = mge_intervals$mge_end,
    length = mge_intervals$mge_end - mge_intervals$mge_start + 1L,
    bin = vapply(bin_rows_per_mge, describe_bin_span, character(1)),
    bin_gene = vapply(bin_rows_per_mge, describe_bin_genes, character(1)),
    stringsAsFactors = FALSE
  )
  
  message("Inferred ", nrow(mges_df), " MGE(s) totalling ",
          format(sum(mges_df$length), big.mark = ","), " bp.")
  
  ## Write the table once, not once per MGE ----
  saveRDS(mges_df, file.path(output_dir, paste0(input_genome_name, "_mges.RDS")))
  write.csv(mges_df, file.path(output_dir, paste0(input_genome_name, "_mges.csv")),
            quote = FALSE, row.names = FALSE)
  
  ## Per-MGE entropy tables ----
  mge_entropy_dir <- file.path(output_dir, "mge_entropy")
  dir.create(mge_entropy_dir, recursive = TRUE, showWarnings = FALSE)
  
  write_one_mge_entropy <- function(mge_row) {
    write.csv(
      data.frame(
        mge_position = seq_len(mges_df$length[mge_row]),
        entropy = entropy_vector[mges_df$start[mge_row]:mges_df$end[mge_row]]
      ),
      file.path(mge_entropy_dir,
                paste0(input_genome_name, "_MGE_", mge_row, "_entropy.csv")),
      quote = FALSE, row.names = FALSE
    )
    invisible(NULL)
  }
  
  run_tasks_in_parallel(as.list(seq_len(nrow(mges_df))), write_one_mge_entropy, ncores)
  
  ## MGE sequences ----
  # The genome is converted to a character string exactly once. The original
  # did this inside the loop, i.e. once per MGE.
  concat_fasta <- Biostrings::readDNAStringSet(
    file.path(output_dir, paste0(input_genome_name, "_concat.fasta"))
  )
  concat_sequence_text <- as.character(concat_fasta[[1]])
  
  mge_sequences <- substring(concat_sequence_text, mges_df$start, mges_df$end)
  
  # rbind + as.vector interleaves header, sequence, header, sequence, ...
  fasta_lines <- as.vector(rbind(paste0(">", mges_df$mge_index), mge_sequences))
  writeLines(fasta_lines,
             file.path(output_dir, paste0(input_genome_name, "_mges_seq.fasta")))
  
  mges_df
}

# Step 7: per-bin summary table ----
#
# Coverage and entropy statistics for every bin, plus its SNP count and any MGE
# it touches.
#
# The original filtered the full genome-length vectors once per bin, which is
# O(bins x genome). Because both tables cover positions 1..genome_length in
# order, a bin's values are simply vector[start:end], and the SNP count comes
# from a prefix sum in constant time.

get_bin_sum <- function(input_genome_name,
                        chr_bins,
                        position_coverage,
                        snps_sum,
                        entropy_df,
                        mges_df,
                        output_dir,
                        ncores) {
  
  genome_length <- max(position_coverage$position)
  
  coverage_vector <- scatter_to_position_vector(position_coverage$position,
                                                as.numeric(position_coverage$coverage),
                                                genome_length, fill = 0)
  entropy_vector <- scatter_to_position_vector(entropy_df$position,
                                               as.numeric(entropy_df$entropy),
                                               genome_length, fill = 0)
  
  # Prefix sum: SNP calls in [a, b] == cumulative[b + 1] - cumulative[a].
  snp_calls_per_position <- tabulate(snps_sum$snp_pos, nbins = genome_length)
  cumulative_snp_calls <- c(0L, cumsum(snp_calls_per_position))
  
  ## Which MGEs does each bin touch? ----
  # Reported as MGE row numbers, matching the original. A bin overlapping
  # several MGEs gets them semicolon-joined -- the original returned a vector
  # here, which silently duplicated the bin's row in the output table.
  mge_label_per_bin <- rep("non-mge", nrow(chr_bins))
  
  if (nrow(mges_df) > 0) {
    bin_mge_overlaps <- IRanges::findOverlaps(
      IRanges::IRanges(start = chr_bins$start, end = chr_bins$end),
      IRanges::IRanges(start = mges_df$start, end = mges_df$end)
    )
    if (length(bin_mge_overlaps) > 0) {
      mge_rows_by_bin <- split(S4Vectors::subjectHits(bin_mge_overlaps),
                               S4Vectors::queryHits(bin_mge_overlaps))
      mge_label_per_bin[as.integer(names(mge_rows_by_bin))] <-
        vapply(mge_rows_by_bin, paste, character(1), collapse = ";")
    }
  }
  
  ## Statistics for one bin ----
  summarise_one_bin <- function(bin_row) {
    
    bin_start <- chr_bins$start[bin_row]
    bin_end <- chr_bins$end[bin_row]
    
    bin_coverage <- coverage_vector[bin_start:bin_end]
    bin_entropy <- entropy_vector[bin_start:bin_end]
    
    bin_has_entropy <- length(bin_entropy) > 0
    
    data.frame(
      bin_index = chr_bins$bin_index[bin_row],
      coverage_max = max(bin_coverage),
      coverage_min = min(bin_coverage),
      coverage_mean = mean(bin_coverage),
      coverage_sd = sd(bin_coverage),
      coverage_median = median(bin_coverage),
      coverage_q1 = as.numeric(quantile(bin_coverage, 0.25)),
      coverage_q3 = as.numeric(quantile(bin_coverage, 0.75)),
      entropy_max = if (bin_has_entropy) max(bin_entropy) else NA_real_,
      entropy_min = if (bin_has_entropy) min(bin_entropy) else NA_real_,
      entropy_mean = if (bin_has_entropy) mean(bin_entropy) else NA_real_,
      entropy_sd = if (bin_has_entropy) sd(bin_entropy) else NA_real_,
      entropy_median = if (bin_has_entropy) median(bin_entropy) else NA_real_,
      entropy_q1 = if (bin_has_entropy) as.numeric(quantile(bin_entropy, 0.25)) else NA_real_,
      entropy_q3 = if (bin_has_entropy) as.numeric(quantile(bin_entropy, 0.75)) else NA_real_,
      total_snp_count = cumulative_snp_calls[bin_end + 1L] - cumulative_snp_calls[bin_start],
      mge_id = mge_label_per_bin[bin_row],
      stringsAsFactors = FALSE
    )
  }
  
  # Per-bin work is small, so hand each worker a block of bins rather than one.
  bin_stats <- data.table::rbindlist(
    run_tasks_in_parallel(
      split_into_chunks(nrow(chr_bins), ncores * 4L),
      function(bin_rows) data.table::rbindlist(lapply(bin_rows, summarise_one_bin)),
      ncores
    )
  )
  
  bin_summary <- merge(chr_bins, bin_stats, by = "bin_index", all = TRUE)
  
  write.csv(bin_summary,
            file.path(output_dir, paste0(input_genome_name, "_bin_summary.csv")),
            quote = FALSE, row.names = FALSE)
  
  invisible(bin_summary)
}

# Step 8: the figure ----
#
# Top panel: % coverage along the concatenated genome, MGEs shaded blue.
# Bottom panel: per-site entropy, coloured by coding vs intergenic, with
# alternating grey blocks marking contig boundaries. Both panels carry a
# marginal histogram of their y values.
#
# The only change is how each position is labelled Coding or Intergenic. The
# original ran a which() over all bins for every base pair; here one coverage
# mask over the coding bins answers it for the whole genome at once. A position
# is Coding when at least one coding bin covers it, identical to the original
# rule "Intergenic unless some overlapping bin is not non-cds".

get_visual <- function(input_genome_name,
                       new_pos,
                       chr_bins,
                       position_coverage,
                       entropy_df,
                       mges_df,
                       output_dir) {
  
  ctg_pos_df <- new_pos$ctg_pos_df
  genome_length <- new_pos$concat_size
  
  ## Coverage panel ----
  coverage_plot_df <- position_coverage %>%
    mutate(cov_pct = 100 * coverage / max(position_coverage$coverage))
  
  coverage_plot <- ggplot() +
    geom_rect(data = mges_df,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "#41b6e6", color = "transparent", linewidth = 1, alpha = 0.65) +
    # Invisible points: ggMarginal needs an x/y mapping it can pick up.
    geom_point(data = coverage_plot_df, aes(x = position, y = cov_pct), alpha = 0) +
    geom_line(data = coverage_plot_df, aes(x = position, y = cov_pct),
              linewidth = 0.1, color = "black", alpha = 0.85) +
    scale_y_continuous(name = "% Coverage") +
    scale_x_continuous(expand = c(0, 0)) +
    theme(
      axis.title.y.left = element_text(colour = "black", size = 20, face = "bold"),
      axis.text.y.left = element_text(colour = "black", size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "transparent"),
      panel.background = element_rect(fill = "transparent"),
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
  
  coverage_plot_margin <- ggMarginal(
    p = coverage_plot, type = "histogram", margins = "y",
    groupFill = FALSE, groupColour = FALSE, size = 5,
    yparams = list(linewidth = 0.15, binwidth = 1, position = "identity")
  )
  
  ## Label every position Coding or Intergenic in one pass ----
  coding_bin_rows <- chr_bins$gene != "non-cds"
  
  if (any(coding_bin_rows)) {
    coding_ranges <- IRanges::IRanges(start = chr_bins$start[coding_bin_rows],
                                      end = chr_bins$end[coding_bin_rows])
    position_is_coding <- as.vector(IRanges::coverage(coding_ranges, width = genome_length)) > 0
  } else {
    position_is_coding <- rep(FALSE, genome_length)
  }
  
  entropy_plot_df <- as.data.frame(entropy_df)
  entropy_plot_df$category <- ifelse(position_is_coding[entropy_plot_df$position],
                                     "Coding", "Intergenic")
  
  # The shrinkage estimator can return values a hair below zero through
  # floating-point rounding; clamp so the bars start at the axis.
  entropy_plot_df$entropy[!is.na(entropy_plot_df$entropy) & entropy_plot_df$entropy < 0] <- 0
  
  max_entropy <- max(entropy_plot_df$entropy, na.rm = TRUE)
  if (!is.finite(max_entropy) || max_entropy <= 0) max_entropy <- 1
  
  ## Entropy panel ----
  contig_shading_df <- ctg_pos_df %>%
    mutate(
      x_start = lag(new_end, default = 0),
      x_end = new_end,
      fill_color = ifelse(row_number() %% 2 == 1, "#808285", "white")
    )
  
  entropy_plot <- ggplot() +
    geom_rect(data = contig_shading_df,
              aes(xmin = x_start, xmax = x_end, ymin = -Inf, ymax = 0, fill = fill_color),
              color = "transparent", alpha = 0.75) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    geom_col(data = entropy_plot_df,
             aes(x = position, y = entropy, fill = category),
             position = "identity", width = 10, alpha = 0.75) +
    geom_point(data = entropy_plot_df,
               aes(x = position, y = entropy, color = category),
               size = 0.1, alpha = 0, show.legend = FALSE) +
    geom_segment(data = ctg_pos_df,
                 aes(x = new_end, xend = new_end, y = 0, yend = -Inf),
                 linetype = "solid", color = "transparent",
                 linewidth = 0.5, alpha = 0.75) +
    scale_color_manual(values = c("Intergenic" = "#E2A4C6", "Coding" = "#91a01e")) +
    scale_fill_manual(values = c("Intergenic" = "#E2A4C6", "Coding" = "#91a01e")) +
    labs(x = "Concatenated Contig Position (bp)", y = "Entropy", fill = "Position Type") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous() +
    theme(
      axis.text.x = element_text(size = 15),
      axis.line.x = element_line(),
      axis.line.y = element_line(),
      axis.text.y = element_text(angle = 0, size = 15),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 20, face = "bold"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 7.5),
      legend.title = element_text(size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
      legend.position = "inside",
      legend.position.inside = c(0.95, 0.9),
      legend.background = element_rect(colour = NA, fill = NA)
    )
  
  entropy_plot_margin <- ggMarginal(
    p = entropy_plot, type = "histogram", margins = "y",
    groupFill = TRUE, groupColour = TRUE, size = 5,
    yparams = list(linewidth = 0.15, binwidth = max_entropy / 200, position = "identity")
  )
  
  combine_plot <- plot_grid(coverage_plot_margin, entropy_plot_margin,
                            ncol = 1, align = "v", axis = "tblr",
                            rel_heights = c(1, 3))
  
  pdf(file = file.path(output_dir, paste0(input_genome_name, "_profiler_plot.pdf")),
      width = 15, height = 7.5)
  print(combine_plot)
  dev.off()
  
  # Kept so the figure can be restyled later without re-running the pipeline.
  plot_df <- list(coverage_plot_df = coverage_plot_df,
                  entropy_plot_df = entropy_plot_df)
  
  saveRDS(plot_df, file.path(output_dir, paste0(input_genome_name, "_plot_df.RDS")))
  
  # ggMarginal leaves a stray default device file behind.
  if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
  
  invisible(plot_df)
}

# Main driver ----
#
# Every step checks for its own cached output first. Deleting one RDS file and
# re-running recomputes that step and everything after it that depends on it,
# which is how to iterate on a single step without repeating the whole run.

## Derive the short genome name used in every output filename ----
derive_genome_name <- function(input_genome_path) {
  if (grepl("GCA_", input_genome_path)) {
    # GenBank assemblies: keep just the accession, e.g. GCA_000013425.1
    paste0("GCA_", strsplit(gsub("\\.fna", "", basename(input_genome_path)),
                            split = "_")[[1]][2])
  } else {
    gsub("\\.fna$|\\.fasta$", "", basename(input_genome_path))
  }
}

## Add the file.path to the files

output_path <- function(suffix) {
  file.path(output_dir, paste0(input_genome_name, suffix))
}

profiler <- function(input_genome_path,
                     output_dir,
                     gff_path,
                     snp_dir,
                     min_ctg_len,
                     cov_cutoff,
                     min_mge_size = 100,
                     ncores) {
  
  data.table::setDTthreads(ncores)
  setwd(output_dir)
  
  input_genome_name <- derive_genome_name(input_genome_path)
  
  
  ## Step 1: concatenate ----
  if (file.exists(output_path("_new_pos.RDS"))) {
    message("new_pos found. Loading RDS.")
    new_pos <- readRDS(output_path("_new_pos.RDS"))
  } else {
    message("new_pos not found. Generating.")
    new_pos <- concat_genome(input_genome_name = input_genome_name,
                             input_genome_path = input_genome_path,
                             min_ctg_len = min_ctg_len,
                             output_dir = output_dir)
    message("Finished generating new_pos.")
  }
  
  ## Step 2: coverage ----
  if (file.exists(output_path("_position_coverage.RDS"))) {
    message("position_coverage found. Loading RDS.")
    position_coverage <- readRDS(output_path("_position_coverage.RDS"))
  } else {
    message("position_coverage not found. Generating.")
    position_coverage <- get_position_coverage(input_genome_name = input_genome_name,
                                               new_pos = new_pos,
                                               snp_dir = snp_dir,
                                               output_dir = output_dir,
                                               ncores = ncores)
    message("Finished generating position_coverage.")
  }
  
  ## Step 3: SNPs ----
  if (file.exists(output_path("_snps_sum.RDS"))) {
    message("snps_sum found. Loading RDS.")
    snps_sum <- readRDS(output_path("_snps_sum.RDS"))
  } else {
    message("snps_sum not found. Generating.")
    snps_sum <- get_snps_sum(input_genome_name = input_genome_name,
                             new_pos = new_pos,
                             snp_dir = snp_dir,
                             output_dir = output_dir,
                             ncores = ncores)
    message("Finished generating snps_sum.")
  }
  
  ## Step 4: bins ----
  if (file.exists(output_path("_chr_bins.RDS"))) {
    message("chr_bins found. Loading RDS.")
    chr_bins <- readRDS(output_path("_chr_bins.RDS"))
  } else {
    message("chr_bins not found. Generating.")
    chr_bins <- get_chr_bins(input_genome_name = input_genome_name,
                             new_pos = new_pos,
                             gff_path = gff_path,
                             output_dir = output_dir)
    message("Finished generating chr_bins.")
  }
  
  ## Step 5: entropy ----
  if (file.exists(output_path("_entropy.RDS"))) {
    message("entropy_df found. Loading RDS.")
    entropy_df <- readRDS(output_path("_entropy.RDS"))
  } else {
    message("entropy_df not found. Generating.")
    entropy_df <- get_entropy(input_genome_name = input_genome_name,
                              new_pos = new_pos,
                              snps_sum = snps_sum,
                              position_coverage = position_coverage,
                              output_dir = output_dir,
                              ncores = ncores)
    message("Finished generating entropy_df.")
  }
  
  ## Step 6: MGEs ----
  if (file.exists(output_path("_mges.RDS"))) {
    message("mges_df found. Loading RDS.")
    mges_df <- readRDS(output_path("_mges.RDS"))
  } else {
    message("mges_df not found. Generating.")
    mges_df <- get_mges(input_genome_name = input_genome_name,
                        entropy_df = entropy_df,
                        new_pos = new_pos,
                        position_coverage = position_coverage,
                        chr_bins = chr_bins,
                        min_mge_size = min_mge_size,
                        cov_cutoff = cov_cutoff,
                        output_dir = output_dir,
                        ncores = ncores)
    message("Finished generating mges_df.")
  }
  
  ## Step 7: bin summary ----
  if (!file.exists(output_path("_bin_summary.csv"))) {
    message("bin_summary.csv not found. Generating.")
    get_bin_sum(input_genome_name = input_genome_name,
                chr_bins = chr_bins,
                position_coverage = position_coverage,
                snps_sum = snps_sum,
                entropy_df = entropy_df,
                mges_df = mges_df,
                output_dir = output_dir,
                ncores = ncores)
    message("Finished generating bin_summary.csv.")
  }
  
  ## Step 8: figure ----
  if (!file.exists(output_path("_profiler_plot.pdf"))) {
    message("Visualization not found. Generating.")
    get_visual(input_genome_name = input_genome_name,
               new_pos = new_pos,
               chr_bins = chr_bins,
               position_coverage = position_coverage,
               entropy_df = entropy_df,
               mges_df = mges_df,
               output_dir = output_dir)
    message("Finished visualization.")
  }
  
  invisible(TRUE)
}

# Execute the function ----

profiler(input_genome_path = input_genome_path,
         output_dir = output_dir,
         gff_path = gff_path,
         snp_dir = snp_dir,
         min_ctg_len = min_ctg_len,
         cov_cutoff = cov_cutoff,
         min_mge_size = min_mge_size,
         ncores = ncores)