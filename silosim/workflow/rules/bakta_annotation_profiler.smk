rule bakta_annotation_profiler:
    conda:
        os.path.join(BASE_PATH,"envs/bakta.yaml")
    input:
        genome_path = input_genome_path,
        bakta_db = [os.path.join(config["output"],"bakta_db","bakta.db")] if config["bakta_db_path"] == "None" else []
    threads:
        config["threads"]
    output:
        gff3 = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.gff3"),
        fna = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.fna"),
        faa = os.path.join(config["output"], "bakta_annotation",genome_name,f"{genome_name}.faa")
    params:
        genome_name = genome_name,
        species = config["species"],
        output_dir = config["output"],
        db_path = [config["bakta_db_path"] if config["bakta_db_path"] != "None" else os.path.join(config["output"],"bakta_db")]
    shell:
        r"""
        mkdir -p {params.output_dir}/bakta_annotation/
        mkdir -p {params.output_dir}/bakta_annotation/tmpdir

        genome_name=({params.genome_name})
        genome_path=({input.genome_path})

        species="{params.species}"

        # Define species
        # When --species is not provided (local query mode), --genus/--species are
        # omitted and Bakta falls back to its taxon-agnostic defaults.
        if [ -z "${{species}}" ] || [ "${{species}}" == "None" ]; then
            taxon_args=""
            echo "No species provided; running Bakta without --genus/--species."
        elif [ "${{species}}" == "sau" ]; then
            taxon_args="--genus Staphylococcus --species aureus"
        elif [ "${{species}}" == "cdiff" ]; then
            taxon_args="--genus Clostridioides --species difficile"
        elif [ "${{species}}" == "kp" ]; then
            taxon_args="--genus Klebsiella --species pneumoniae"
        elif [ "${{species}}" == "sepi" ]; then
            taxon_args="--genus Staphylococcus --species epidermidis"
        elif [ "${{species}}" == "bp" ]; then
            taxon_args="--genus Bordetella --species pertussis"
        else
            echo "ERROR: Unknown species ${{species}} for Bakta annotation." >&2
            exit 1
        fi

        # Run Bakta annotation
        bakta --force --db {params.db_path} \
            --output {params.output_dir}/bakta_annotation/{params.genome_name}/ \
            --threads {threads} \
            --prefix {params.genome_name} \
            --tmp-dir {params.output_dir}/bakta_annotation/tmpdir \
            ${{taxon_args}} \
            {input.genome_path}

        rm -rf {params.output_dir}/bakta_annotation/tmpdir
        """