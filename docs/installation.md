# **Installation**
SILOSim utilizes Snakemake pipelines complemented by Python scripts for input validation, configuration file generation, and execution of the workflow.

There are two ways to install SILOSim:
- **[Bioconda](#bioconda-installation-recommended)** (recommended): one command, no need to download the source code.
- **[Manual installation via Git](#manual-installation-via-git)**: for developers or anyone who wants the latest code from GitHub.

## **Before You Start: Install Conda**
SILOSim is installed with **conda**, a free tool that downloads software and keeps it in separate, self-contained "environments" so different programs don't interfere with each other.

To check whether conda is already installed, open a terminal and run:
```bash
conda --version
```
If you see a version number (for example, `conda 24.9.2`), you are ready to go. If you see `command not found`, install [Miniforge](https://github.com/conda-forge/miniforge) first by following the instructions on its page, then close and reopen your terminal.

## **Bioconda Installation (Recommended)**
SILOSim is available on [Bioconda](https://bioconda.github.io/recipes/silosim/README.html#package-silosim), a collection of bioinformatics software for conda.

1. **Set up the Bioconda channels (one time only).**
   Channels are the online sources conda downloads software from. Run these three commands once per computer; you can skip this step if you have used Bioconda before.
   ```bash
   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   ```

2. **Create a conda environment and install SILOSim into it:**
   ```bash
   conda create -n silosim silosim
   ```
   This creates a new environment named `silosim` (the `-n silosim` part) and installs the `silosim` package into it. Conda will show a list of packages it is about to install and ask `Proceed ([y]/n)?`. Type `y` and press Enter. This may take a few minutes.

3. **Activate the environment and display the help message:**
   ```bash
   conda activate silosim
   silosim -h
   ```
   If the help message appears, SILOSim is installed correctly.

> **Note:** You need to run `conda activate silosim` each time you open a new terminal before using SILOSim.

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

## **What Is Downloaded and Installed?**
Both installation methods create a conda environment named `silosim` containing SILOSim's core scripts, its configuration files (YAML files that specify which bioinformatics tools each mode needs), and Snakemake with its dependencies.
- With **Bioconda**, conda downloads all of this for you in a single step.
- With **Git**, cloning the repository downloads the scripts and configuration files, and the installation script then creates the `silosim` environment and installs SILOSim and Snakemake into it.

The bioinformatics tools themselves are not installed during this installation. The first time you run a mode (profiler or simulator), Snakemake will automatically download and install only the tools required for that specific mode. This minimizes disk space usage, but requires an internet connection when running the mode.

Think of it like having a recipe book at home. You only buy the specific ingredients when you decide to cook a particular dish, rather than stocking your pantry with every ingredient for every recipe.
