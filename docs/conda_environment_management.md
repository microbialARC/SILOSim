# Conda Environment Management
SILOSim uses Conda and Snakemake to automatically manage dependencies and create required environments. All four modes share the same Conda environments. After the initial installation, subsequent runs can reuse existing environments by specifying the `--conda_prefix` flag in the SILOSim command line, eliminating the need to recreate environments.

For example, run Full mode with previously created conda environments:
```
silosim profiler --input_genome /path/to/genome.fasta --species <species> --conda_prefix /path/to/conda/environments/
```
