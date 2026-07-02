#!/usr/bin/env bash
set -e
# Parse command line arguments for forece reinstall option
FORCE_REINSTALL=false
while [[ $# -gt 0 ]]; do
    case $1 in
        -f|--force)
            FORCE_REINSTALL=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -f, --force    Force reinstall in existing silosim environment"
            echo "  -h, --help     Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

echo "Starting SILOSim installation..."
# Detect if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed. Please install Miniconda or Anaconda first."
    exit 1
fi

# Check if GNU parallel is installed
# Try loading module first
module load parallel 2>/dev/null || true

if ! command -v parallel &> /dev/null; then
    echo "GNU parallel is not installed."
    echo "Please install it before running silosim:"
    echo ""
    echo "  Ubuntu/Debian:  sudo apt install parallel"
    echo "  Fedora/RHEL:    sudo dnf install parallel"
    echo ""
    exit 1
fi

# Initialize conda for bash shell
if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)" || . "$(conda info --base)/etc/profile.d/conda.sh"
fi

# Get the directory of the install.sh script
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Check if silosim environment exists
ENV_EXISTS=false
if conda env list | grep -q '^silosim\s'; then
    ENV_EXISTS=true
fi

# Handle installation based on environment existence and force flag
if [ "$ENV_EXISTS" = true ] && [ "$FORCE_REINSTALL" = false ]; then
    echo "Detected existing conda environment 'silosim' but --force flag not set."
    echo "Quitting installation to avoid overwriting existing environment."
    echo "If the environment is corrupted or dependencies are missing, please remove the existing silosim environment and rerun this script."
    echo "If you want to reinstall silosim by overwriting the existing environment, use:"
    echo "bash install.sh --force"
    exit 1
elif [ "$ENV_EXISTS" = true ] && [ "$FORCE_REINSTALL" = true ]; then
    echo "Force reinstalling in existing conda environment 'silosim'..."
elif [ "$ENV_EXISTS" = false ] && [ "$FORCE_REINSTALL" = true ]; then
    echo "--force flag detected but silosim environment does not exist."
    echo "Creating conda environment 'silosim'..."
    conda env create -f "$INSTALL_DIR/SILOSim.yml"
elif [ "$ENV_EXISTS" = false ] && [ "$FORCE_REINSTALL" = false ]; then
    # normal first install
    echo "Creating conda environment 'silosim'..."
    conda env create -f "$INSTALL_DIR/silosim.yml"
fi

# Activate the silosim conda environment
printf "\nActivating silosim environment...\n"
conda activate silosim

# Double check if the correct conda environment is activated
CURRENT_ENV=$(basename "$CONDA_PREFIX")
if [ "$CURRENT_ENV" != "silosim" ]; then
    echo "Wrong conda environment detected. Current environment: $CURRENT_ENV"
    echo "Please activate the silosim environment: conda activate silosim"
    exit 1
fi

# Install SILOSim with verbose info
echo "Installing SILOSim package..."
pip install "$INSTALL_DIR" -v

# Print success message
cols=$(tput cols 2>/dev/null || echo "${COLUMNS:-80}")
printf '%*s\n' "$cols" '' | tr ' ' '='

if color="$(tput setaf 3 2>/dev/null)"; then
    reset="$(tput sgr0 2>/dev/null)"
else
    color=$'\033[0;33m'
    reset=$'\033[0m'
fi

printf '%b\n' "${color}silosim installed successfully!${reset}"
printf '%b\n' "${color}Activate the environment with: conda activate silosim${reset}"
printf '%b\n' "${color}Then run 'silosim -h' to see the help page.${reset}"
printf '%*s\n' "$cols" '' | tr ' ' '='