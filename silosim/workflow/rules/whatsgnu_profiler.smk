rule whatsgnu_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/whatsgnu.yaml")
    input: 
        faa = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.faa")
    params:
        species = config["species"],
        genome_name = genome_name,
        db_path = config["whatsgnu_db_path"],
        output_dir = config["output"],
        top_genome = config["top_genomes"],
        local_query_dir = config["local_query_dir"]
    output:
        whatsgnu_output = os.path.join(config["output"], "whatsgnu", genome_name, f"{genome_name}_WhatsGNU_topgenomes.txt")
    shell:
        """
        # Add defensive bash options to catch errors early
        set -euo pipefail
        
        genome_name="{params.genome_name}"
        input_faa_path="{input.faa}"
        local_query_dir="{params.local_query_dir}"
        out_file="{output.whatsgnu_output}"
        db_path=""

        mkdir -p "$(dirname "${{out_file}}")"
        if [ "${{local_query_dir}}" != "None" ]; then
            # LOCAL MODE - skip WhatsGNU and NCBI entirely.
            # The topgenomes file becomes the list of user genomes.
            local_query_dir="{params.local_query_dir}"
            : > "${{out_file}}"
            for f in "${{local_query_dir}}"/*.fasta "${{local_query_dir}}"/*.fa "${{local_query_dir}}"/*.fna; do
                [ -e "$f" ] || continue
                printf '%s\n' "$f" >> "${{out_file}}"
            done
        else
            # PUBLIC MODE - run WhatsGNU to fetch top genomes from NCBI
            
            if [ "{params.db_path}" == "None" ]; then
                mkdir -p {params.output_dir}/whatsgnu/db
                if [ "{params.species}" == "sau" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13388052/files/Saureus.tar.gz?download=1 
                    tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/Saureus/Sau_WhatsGNU_Ortholog_db.pickle"
                    whatsgnu_mode="ortholog"
                elif [ "{params.species}" == "cdiff" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13387715/files/Cdiff.tar.gz?download=1
                    tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/Cdiff/Cdiff_WhatsGNU_Ortholog_db.pickle"
                    whatsgnu_mode="ortholog"
                elif [ "{params.species}" == "kp" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.zip https://zenodo.org/record/7812697/files/Kp.zip?download=1
                    unzip -o {params.output_dir}/whatsgnu/db/whatsgnu_db.zip -d {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/Kp_Ortholog_8752.pickle"
                    whatsgnu_mode="ortholog"
                elif [ "{params.species}" == "sepi" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.gz https://zenodo.org/records/14751549/files/Sepi_WhatsGNU_basic.txt.gz?download=1
                    mkdir -p {params.output_dir}/whatsgnu/db
                    gunzip -c {params.output_dir}/whatsgnu/db/whatsgnu_db.gz > {params.output_dir}/whatsgnu/db/Sepi_WhatsGNU_basic.txt
                    db_path="{params.output_dir}/whatsgnu/db/Sepi_WhatsGNU_basic.txt"
                    whatsgnu_mode="basic"
                elif [ "{params.species}" == "ecoli" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13693568/files/Ecoli_basic.tar.gz?download=1
                    tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/Ecoli/basic_Ecoli.txt"
                    whatsgnu_mode="basic"
                elif [ "{params.species}" == "pa" ]; then
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz https://zenodo.org/records/13385380/files/PA_basic.tar.gz?download=1
                    tar -xvzf {params.output_dir}/whatsgnu/db/whatsgnu_db.tar.gz -C {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/PA/PA_basic_db.txt"
                    whatsgnu_mode="basic"
                elif [ "{params.species}" == "mtb" ]; then
                    # The Dropbox link is a shared folder, which downloads as a zip archive
                    wget -O {params.output_dir}/whatsgnu/db/whatsgnu_db.zip "https://www.dropbox.com/sh/8nqowtd4fcf7dgs/AAAdXiqcxTsEqfIAyNE9TWwRa?dl=1"
                    unzip -o {params.output_dir}/whatsgnu/db/whatsgnu_db.zip TB_Ortholog_6563.pickle -d {params.output_dir}/whatsgnu/db
                    db_path="{params.output_dir}/whatsgnu/db/TB_Ortholog_6563.pickle"
                    whatsgnu_mode="ortholog"
                fi
            elif [ "{params.db_path}" != "None" ]; then
                db_path="{params.db_path}"
                case "{params.species}" in
                    sepi|ecoli|pa)
                        whatsgnu_mode="basic"
                        ;;
                    *)
                        whatsgnu_mode="ortholog"
                        ;;
                esac
            fi
            
            cd {params.output_dir}/whatsgnu/

            # Import the top genomes number 
            top_genomes="{params.top_genome}"

            WhatsGNU_main.py \
            -d "${{db_path}}" \
            -dm "${{whatsgnu_mode}}" \
            --force \
            -t \
            -tn "${{top_genomes}}" \
            -o "${{genome_name}}" \
            "${{input_faa_path}}"

            # The kp and mtb databases name genomes as <strain>_<GenBank accession>[_<ST>_]
            # (e.g. ST_15_K087_GCA_002417375.1_ST15_), which NCBI Datasets cannot download.
            # Keep only the GenBank accession in column 1 so the downstream rules see bare accessions.
            if [ "{params.species}" == "kp" ] || [ "{params.species}" == "mtb" ]; then
                awk 'BEGIN {{FS = OFS = "\t"}}
                     NR > 1 && match($1, /GCA_[0-9]+[.][0-9]+/) {{$1 = substr($1, RSTART, RLENGTH)}}
                     {{print}}' "${{out_file}}" > "${{out_file}}.tmp"
                mv "${{out_file}}.tmp" "${{out_file}}"
            fi
        fi
        """