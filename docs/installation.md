# **Installation**
SILOSim utilizes Snakemake pipelines complemented by Python scripts for input validation, configuration file generation, and execution of the workflow.

## **Bioconda Installation**  
*Coming soon.*

## **Manual Installation via Git**  
1. Clone the GitHub repository:  
   ```bash
   git clone https://github.com/microbialARC/SILOSim
   cd SILOSim
   ```

2. Install SILOSim using the installation script:
    ```bash
    bash install.sh
    ```
3. Activate SILOSim environment and display help message:
    ```bash
    conda activate silosim
    silosim -h
    ```

### What is downloaded and installed?
When you clone the GitHub repository, you receive SILOSim's core scripts and configuration files (YAML files) that specify which bioinformatics tools are needed for each function. The installation script then creates a conda environment named `silosim` and installs SILOSim along with Snakemake and its dependencies in `silosim` conda environment.

The bioinformatics tools themselves are not installed during this installation. When you first run a function (e.g., Strain Identifier or Genome Profiler), Snakemake will automatically download and install only the tools required for that specific function. This minimizes disk space usage, but requires internet connection when running the function.

Think of it like having a recipe book at home. You only buy the specific ingredients when you decide to cook a particular dish, rather than stocking your pantry with every ingredient for every recipe.