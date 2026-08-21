rule dataset_topgenomes_profiler:
    conda:
        os.path.join(BASE_PATH, "envs/datasets.yaml")
    input:
        topgenomes_result = os.path.join(config["output"], "whatsgnu", genome_name, f"{genome_name}_WhatsGNU_topgenomes.txt")
    output:
        actual_download_topgenomes = os.path.join(config["output"], "datasets_topgenomes", "actual_download_topgenomes.txt")
    threads:
        config["threads"]
    params:
        output_dir = os.path.join(config["output"], "datasets_topgenomes"),
        local_query_dir = config["local_query_dir"]
    shell:
        r"""
        mkdir -p {params.output_dir}
        topgenomes_file="{input.topgenomes_result}"
        local_query_dir="{params.local_query_dir}"

        if [ "${{local_query_dir}}" != "None" ]; then
            # LOCAL MODE. Skip WhatsGNU and NCBI Datasets entirely.
            # The input file already holds absolute paths to the user's FASTAs.
            # Link them into the same directory the public branch populates, using
            # the same <stem>.fna convention, and emit bare stems so every
            # downstream rule stays mode-agnostic.

            # These symlinks are NOT declared Snakemake outputs, so a failed or
            # interrupted run leaves them behind. Clear them first so this rule is
            # idempotent and the duplicate check below only sees the current run.
            rm -f {params.output_dir}/*.fna

            : > {output.actual_download_topgenomes}

            while IFS= read -r genome_path; do
                [ -n "${{genome_path}}" ] || continue
                base=$(basename "${{genome_path}}")
                stem="${{base%.*}}"          # strip .fasta / .fa / .fna

                # Compare against stems recorded during THIS run, not against disk
                if grep -qxF "${{stem}}" {output.actual_download_topgenomes}; then
                    echo "ERROR: two files in ${{local_query_dir}} reduce to the same" >&2
                    echo "       genome name '${{stem}}' (e.g. ${{stem}}.fasta and ${{stem}}.fna)." >&2
                    echo "       Keep one file per genome, or rename them uniquely." >&2
                    exit 1
                fi

                ln -s "${{genome_path}}" {params.output_dir}/${{stem}}.fna
                printf '%s\n' "${{stem}}" >> {output.actual_download_topgenomes}
            done < "${{topgenomes_file}}"

            n_genomes=$(wc -l < {output.actual_download_topgenomes})
            if [ "${{n_genomes}}" -eq 0 ]; then
                echo "ERROR: no genomes read from ${{topgenomes_file}}" >&2
                exit 1
            fi
            echo "[dataset_topgenomes_profiler] local mode: linked ${{n_genomes}} genome(s), no downloads performed"

        else
            # PUBLIC MODE. Run NCBI Datasets to download the top genomes determined by WhatsGNU.
            mkdir -p {params.output_dir}/scripts

            # Read the topgenomes from WhatsGNU output
            # Compared to the rules and scripts used in THRESHER
            # Instead of using a python script,
            # this rule directly uses shell commands to create the command for downloading top genomes determined by WhatsGNU

            # Skip the first row(tail -n +2) of the whatsgnu topgenomes file and extract the genome accessions at column 1, separated by tab(cut -f1)
            topgenomes=$(tail -n +2 "$topgenomes_file" | cut -f1)

            # For each genome accession, create a script to download the genome using datasets command
            while IFS= read -r genome_entry; do
                echo "#!/bin/bash
                datasets download genome accession ${{genome_entry}} --filename {params.output_dir}/${{genome_entry}}.zip
                unzip -o {params.output_dir}/${{genome_entry}}.zip -d {params.output_dir}/${{genome_entry}}/
                mv {params.output_dir}/${{genome_entry}}/ncbi_dataset/data/*/*.fna {params.output_dir}/${{genome_entry}}.fna
                rm {params.output_dir}/${{genome_entry}}.zip
                rm -rf {params.output_dir}/${{genome_entry}}" > {params.output_dir}/scripts/datasets_${{genome_entry}}.sh
            done < <(echo "$topgenomes")

            ls {params.output_dir}/scripts/datasets_*.sh > {params.output_dir}/scripts/script_list.txt
            # This is critical
            # No // in the file containing paths to the scripts otherwise there would be error!!!
            sed -i 's#//#/#g' {params.output_dir}/scripts/script_list.txt
            # Run the scripts in parallel using GNU parallel
            module load parallel
            parallel --silent --jobs {threads} bash :::: {params.output_dir}/scripts/script_list.txt
            rm -rf {params.output_dir}/scripts
            find {params.output_dir} -name "*.fna" | awk -F"/" '{{print $NF}}' | sed 's/\.fna$//g' > {output.actual_download_topgenomes}
        fi
        """