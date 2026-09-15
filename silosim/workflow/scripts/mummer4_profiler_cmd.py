# This script generates MUMmer4 commands for profiling genomes against their top genomes filtered by FastANI.
# There are 7 steps for each genome
# The purpose of each step is explained below:
# 1. Change directory to the output directory so that all output files are saved there
# 2. nucmer to generate the delta file
# 3. show-coords to get the coverage information to get the denominator for calculating entropy
# 4. awk to extract relevant columns from the coords file

    # show-coords output columns:
    # [S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [TAGS]
    # [S1] Start of the alignment region in the reference sequence.
    # [E1] End of the alignment region in the reference sequence.
    # [S2] Start of the alignment region in the query sequence.
    # [E2] End of the alignment region in the query sequence.
    # [LEN 1] Length of the alignment region in the reference sequence,
    # measured in nucleotides.
    # [LEN 2] Length of the alignment region in the query sequence, measured
    # in nucleotides.
    # [TAGS] The reference FastA ID and the query FastA ID.

    # 1st column: [S1] Start of the alignment region in the reference sequence.
    # 2nd column: [E1] End of the alignment region in the reference sequence.
    # 7th column: [TAG] The reference FastA ID of the alignment.
# 5. show-snps to get the coordinates of snps
# 6. awk to extract relevant columns from the snps file
    # show-snps output columns:
    # [P1]    [SUB]   [SUB]   [P2]    [BUFF]  [DIST]  [FRM]   [TAGS]
    # [P1]  SNP position in the reference.
    # [SUB] Character in the reference.
    # [SUB] Character in the query.
    # [P2] SNP position in the query.
    # [BUFF] Distance from this SNP to the nearest mismatch (end of 
    # alignment, indel, SNP, etc) in the same alignment.
    # [DIST] Distance from this SNP to the nearest sequence end.
    # [FRM] Reading frame for the reference sequence and the
    # reading frame for the query sequence respectively. Simply
    # 'forward'  1, or 'reverse' -1 for nucmer data.
    # [TAGS] The reference FastA ID and the query FastA ID.
    
    # 1st column: [P1] SNP position in the reference.
    # 2nd column: [SUB] Character in the reference.
    # 3rd column: [SUB] Character in the query.
    # 4th column: [P2] SNP position in the query.
    # 5th column: [BUFF] Distance from this SNP to the nearest mismatch (end of
    # alignment, indel, SNP, etc) in the same alignment.
    # We use the Good SNPs with sufficient match buffer (20)
    # https://github.com/mummer4/mummer/blob/828c3df5e6c0cdb8dd272969ac54cfebb83f67ac/scripts/dnadiff.pl#L421
    # "if ( $A[4] >= $SNPBuff ) {"
    # 9th column: The reference FastA ID
# 7. Remove the delta file to save space
import os

# Import from snakemake
# Input
reference_genome_path = snakemake.input.input_fna
reference_genome_name = snakemake.params.reference_genome_name
topgenomes_filtered_list = snakemake.input.topgenomes_filtered
pipeline_output_dir = snakemake.params.output_dir

def get_mummer4_commands(reference_genome_name, reference_genome_path, topgenomes_filtered_list, pipeline_output_dir):
    # Path to the scripts dir
    scripts_dir = os.path.join(pipeline_output_dir, 'mummer4', "scripts")
    # Path to mummer4 output dir
    output_dir = os.path.join(pipeline_output_dir, 'mummer4', "output")
    # Make directories if not exist
    os.makedirs(scripts_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    # Use the iteration to create MUMmer4 commands for each filtered top genome against the reference genome
    with open(topgenomes_filtered_list, 'r') as topgenomes_file:
        for line in topgenomes_file:
            topgenome_path = line.strip()
            topgenome_name = os.path.basename(topgenome_path).replace('.fna', '')
            compare_prefix = f"{reference_genome_name}_{topgenome_name}"
            topgenome_cmd = ["#!/bin/bash",
                             f"cd {output_dir}",
                             # Set the --maxmatch and -l 50 flags to avoid missing alignments in repetitive regions and to set the minimum length of a match to 50 bp
                             # This avoids abnormal surge in RAM usage which will cause OOM inevitably
                             f"nucmer --maxmatch -l 50 {reference_genome_path} {topgenome_path} -p {compare_prefix}",
                             f"show-coords -b -r -T -H {compare_prefix}.delta > {compare_prefix}.coords",
                             # Add the column 8 ([TAGS] The reference FastA ID and the query FastA ID) so that we can track the reference sites as well
                             "awk '{print $1, $2, $7, $8}' " + f"{compare_prefix}.coords > tmp_{compare_prefix}.coords && mv tmp_{compare_prefix}.coords {compare_prefix}.coords",
                             f"show-snps -T -C -H {compare_prefix}.delta > tmp_{compare_prefix}.snps && mv tmp_{compare_prefix}.snps {compare_prefix}.snps",
                             # Add the column 4 ([P2] SNP position in the query) so that we can track the query sites as well
                             "awk '$5 >= 20 {print $1, $2, $3, $4, $5, $9}' " + f"{compare_prefix}.snps > tmp_{compare_prefix}.snps && mv tmp_{compare_prefix}.snps {compare_prefix}.snps",
                             f"rm {compare_prefix}.delta"
                             ]
            # Write the command to a script file
            cmd_script_path = os.path.join(scripts_dir, f'mummer4_{reference_genome_name}_{topgenome_name}.sh')
            with open(cmd_script_path, 'w') as cmd_file:
                cmd_file.write('\n'.join(topgenome_cmd) + '\n')

    # Create a text with path to all scripts generated 
    with open(os.path.join(scripts_dir, f'{reference_genome_name}_cmd_list.txt'), 'w') as list_file:
        for script_file in os.listdir(scripts_dir):
            if script_file.startswith('mummer4_') and script_file.endswith('.sh'):
                list_file.write(os.path.join(scripts_dir, script_file) + '\n')
    return 0

# Call the function to generate MUMmer4 commands
get_mummer4_commands(reference_genome_name, reference_genome_path, topgenomes_filtered_list, pipeline_output_dir)