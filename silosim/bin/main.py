#!/usr/bin/env python3
"""SILOSim Main Entry
This script creates the configuration files and executes the Snakemake workflow.
"""
VERSION = "0.1.1-beta"

# Import standard libraries and custom modules
import argparse
import os
import sys
from pathlib import Path
import traceback
# subprocess is needed to call snakemake
import subprocess
# import shutil to get terminal size for centering text
import shutil

# Import input argument validator
from silosim.bin.args_validator import (
    validate_profiler,
    validate_simulator,
    ValidationError
)

# Import config creator
from silosim.bin.config_creator import (
    profiler_config,
    simulator_config
)

# Import argument parsers for different modes
from silosim.bin.parsers.profiler_parser import add_profiler_parser
from silosim.bin.parsers.simulator_parser import add_simulator_parser

def build_parser():
    """Build the main argument parser for the command line interface"""
    description = f"""
    SILOSim: Profile per-site bacterial genome variability and simulate evolution

    Usage:
      silosim <mode> [options]
    """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS
    )

    parser.add_argument(
        "-v", "--version",
        action="version",
        version=f"SILOSim version {VERSION}",
        help="show version and exit"
    )

    subparsers = parser.add_subparsers(
        dest="command",
        metavar="<mode>",
        help="Use 'silosim <mode> -h' for more information on a specific mode."
    )

    # Add subparsers for each mode
    add_profiler_parser(subparsers)
    add_simulator_parser(subparsers)

    return parser

# Helper functions to check OS and RAM since these checks are used multiple times
def check_os(term_width):
    """Check the operating system"""
    print()
    print("=" * term_width)
    print("Checking operating system".center(term_width))
    print("=" * term_width)
    system_os = sys.platform
    if system_os != 'linux':
        print("SILOSim supports Linux only. Required dependencies are currently available only on Linux.")
        print("Please run SILOSim on a Linux system.")
        print("SILOSim quitting...")
        return 1
    else:
        print(f"Current Operating system: {system_os}. Proceeding...")
        return 0

def check_ram(term_width,
              min_ram_gb):
    """Check if the system has at least min_ram_gb of RAM"""
    print()
    print("=" * term_width)
    print("Checking system RAM".center(term_width))
    print("=" * term_width)
    ram_gb =  os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') / (1024 ** 3)
    if ram_gb < min_ram_gb:
        print("At least \033[91m40 GB\033[0m memory is required for SILOSim, as WhatsGNU requires loading large databases into memory.")
        print("Please run SILOSim on a system with sufficient RAM.")
        print("SILOSim quitting...")
        return 1
    else:
        print(f"Detected RAM: \033[92m{int(ram_gb)} GB\033[0m.")
        return 0
    

def execute_snakemake(config_path, snakefile, conda_prefix, threads, output_dir, prefix):
    """
    Execute the Snakemake workflow
    
    Args:
        config_file: Path to the config YAML file
        snakefile: Name of the Snakefile to use
        conda_prefix: Path to conda environments directory
        threads: Number of threads to use
        output_dir: Output directory for resume command scripts
        prefix: Prefix for the output files
        
    Returns:
        int: return code from Snakemake
    """
    # Locate the workflow directory
    snakefile_path = Path(__file__).parent.parent / "workflow" / snakefile
    # Print execution info
    print(f"Config file: {config_path}")
    print(f"Snakefile: {snakefile_path}")
    print(f"Threads: {threads}")
    print(f"Conda environments: {conda_prefix}")
    
    # Build and execute Snakemake command
    snakemake_cmd = [
        "snakemake",
        "--cores", str(threads),
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
        "--use-conda",
        "--conda-prefix", str(conda_prefix)
    ]
    # If accidentally interrupted, allow resuming using resume_cmd
    # --rerun-incomplete is a flag that added to the snakemake_cmd to allow re-running incomplete jobs
    resume_cmd = [
        "snakemake",
        "--cores", str(threads),
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
        "--use-conda",
        "--conda-prefix", str(conda_prefix),
        "--rerun-incomplete"
    ]
    # First unlock the working directory
    unlock_cmd = [
        "snakemake",
        "--unlock",
        "--configfile", str(config_path),
        "--snakefile", str(snakefile_path),
    ]
    # Create resume_cmd directory
    resume_cmd_dir = Path(output_dir) / "resume"
    resume_cmd_dir.mkdir(parents=True, exist_ok=True)

    # Export unlock command script
    unlock_script_path = resume_cmd_dir / f"unlock_snakemake_{prefix}.sh"
    with open(unlock_script_path, 'w') as f:
        f.write("#!/usr/bin/env bash\n\n")
        f.write("# Unlock Snakemake working directory after interruption\n")
        f.write(' '.join(unlock_cmd) + '\n')
    unlock_script_path.chmod(0o755)
    
    # Export resume command script
    resume_script_path = resume_cmd_dir / f"resume_snakemake_{prefix}.sh"
    with open(resume_script_path, 'w') as f:
        f.write("#!/usr/bin/env bash\n\n")
        f.write("# Resume Snakemake workflow after interruption\n")
        f.write(' '.join(resume_cmd) + '\n')
    resume_script_path.chmod(0o755)

    print(f"Executing Snakemake command: \033[93m{' '.join(snakemake_cmd)}\033[0m")
    # note!! 
    # If the snakemake job is interrupted, to resume, run the same command again
    print("\033[33m" + """
NOTE: To resume an interrupted workflow, simply unlock and run the executed Snakemake command with "--rerun-incomplete" flag, highlighted below.
Snakemake will skip completed steps and continue where it left off.
""" + "\033[0m")
    print(f"To unlock the working directory, run: \033[93m{' '.join(unlock_cmd)}\033[0m")
    print()
    print(f"To resume the workflow, run: \033[93m{' '.join(resume_cmd)}\033[0m")
    print("\033[33mNOTE: This is a temporary workflow resumption method. A dedicated resume feature will be implemented in a future release.\033[0m")
    print()
    print(f"\033[33mUnlock and resume scripts saved to: {unlock_script_path}, {resume_script_path}\033[0m")
    print()
    print("If you are using the same function and mode after the first finished run,")
    print("you can reuse the conda environments by specifying the same --conda-prefix path.")
    print("This avoids re-creating conda environments and saves time and disk space.")
    print(f"Conda environments are stored in: \033[93m{conda_prefix}\033[0m")
    print()
    # Flush stdout to ensure all prints appear before Snakemake output
    sys.stdout.flush()
    sys.stderr.flush()
    # Execute with subprocess.run()
    try:
        result = subprocess.run(snakemake_cmd, check=False)
        return result.returncode
    except FileNotFoundError:
        print("Error: Snakemake command not found. Please ensure Snakemake is installed and accessible in your PATH.")
        return 1
    except Exception as e:
        print(f"Error executing Snakemake: {e}")
        return 1

def main():
    """Main entry point for the SILOSim command line interface"""

    # ASCII Art for SILOSim
    ascii_art = """
███████╗██╗██╗      ██████╗ ███████╗██╗           
██╔════╝██║██║     ██╔═══██╗██╔════╝══╝           
███████╗██║██║     ██║   ██║███████╗██╗██╗███████╗
╚════██║██║██║     ██║   ██║╚════██║██║██║ ██╗ ██║
███████║██║███████╗╚██████╔╝███████║██║██║ ╚═╝ ██║
╚══════╝╚═╝╚══════╝ ╚═════╝ ╚══════╝╚═╝╚═╝     ╚═╝
   Version {VERSION}
   Center for Microbial Medicine, Children's Hospital of Philadelphia; Philadelphia, PA, USA
    """
    term_width = shutil.get_terminal_size((80, 20)).columns
    ascii_art = ascii_art.format(VERSION=VERSION)
    centered_art = "\n".join(line.center(term_width) for line in ascii_art.splitlines())
    print(centered_art)
    # Build and parse arguments
    parser = build_parser()
    args = parser.parse_args()

    print(f"SILOSim Mode: {args.command}".center(term_width))
    
    # For profiler, check OS and RAM
    if args.command == "profiler":
        # For profiler, if public query is used, check OS and RAM
        # Otherwise RAM is not checked as WhatsGNU is not used
        if not args.local_query_dir:
            print(f"\nSILOSim Mode: {args.command} with public query genomes. Checking OS and RAM...".center(term_width))
            if not args.force:
                # Check OS and abort early on failure
                if (check_os_rc := check_os(term_width)):
                    return check_os_rc
                # Check if RAM is enough (at least 40 GB for WhatsGNU database loading)
                if (check_ram_rc := check_ram(term_width, 40)):
                    return check_ram_rc
                print(f"RAM check passed for SILOSim Mode: {args.command}. Proceeding...")
            elif args.force:
                print()
                print("=" * term_width)
                print("Bypassing OS and RAM checks".center(term_width))
                print("=" * term_width)
                print(f"\033[93m--force flag set: skipping OS and RAM checks for function: {args.command}.\nProceeding may lead to dependency errors during execution.\033[0m")
        else:
            print(f"\nSILOSim Mode: {args.command} with local query genomes. Checking OS only...".center(term_width))
            if not args.force:
                # Check OS and abort early on failure
                if (check_os_rc := check_os(term_width)):
                    return check_os_rc
                print(f"OS check passed for SILOSim Mode: {args.command}. Proceeding...")
            elif args.force:
                print()
                print("=" * term_width)
                print("Bypassing OS check".center(term_width))
                print("=" * term_width)
                print(f"\033[93m--force flag set: skipping OS check for function: {args.command}.\nProceeding may lead to dependency errors during execution.\033[0m")
    elif args.command == "simulator":
        # For evo_simulator, only check OS
        if not args.force:
            if (check_os_rc := check_os(term_width)):
                return check_os_rc
        elif args.force:
            print()
            print("=" * term_width)
            print("Bypassing OS check. Proceeding...".center(term_width))
            print("=" * term_width)
    
    # Validate arguments based on the selected mode
    print()
    print("=" * term_width)
    print("Validating input arguments".center(term_width))
    print("=" * term_width)

    # Create config file based on the selected mode
    # Initialize variables
    config_path = None
    snakefile = None

    try:
        if args.command == "profiler":
            validate_profiler(args)
            config_path = profiler_config(args)
            snakefile = "Snakefile_profiler"
        elif args.command == "simulator":
            validate_simulator(args)
            config_path = simulator_config(args)
            snakefile = "Snakefile_simulator"
            # For now, if simulator mode, thread default to 1 
            # As the current simulation only supports single-threaded
            args.threads = 1
        else:
            raise ValidationError(f"Unknown SILOsim Mode: {args.command}")
            
        # Ensure config path and snakefile are set
        if config_path is None or snakefile is None:
            raise ValidationError("Configuration path or Snakefile not set properly.")
        # Print a full-width separator
        print()
        print("=" * term_width)
        print("Launching Snakemake workflow".center(term_width))
        print("=" * term_width)

        # Force flush all pending prints before Snakemake execution to ensure proper ordering
        sys.stdout.flush()
        sys.stderr.flush()

        # Execute Snakemake workflow
        # config_path, snakefile, conda_prefix, threads, output_dir, prefix
        returncode = execute_snakemake(
            config_path,
            snakefile,
            os.path.abspath(args.conda_prefix),
            args.threads,
            os.path.abspath(args.output),
            args.prefix
        )

        # Print completion message based on return code
        if returncode == 0:
            print()
            print("=" * term_width)
            print(f"SILOsim Mode: {args.command} completed successfully!".center(term_width))
            print("=" * term_width)
            print()
        else:
            print()
            print("=" * term_width)
            print(f"SILOsim Mode: {args.command} failed with return code {returncode}".center(term_width))
            print("=" * term_width)
            print()

        return returncode
    
    # Handle exceptions and print error messages for debugging purposes
    except ValidationError as e:
        print(f"Validation Error: {e}", file=sys.stderr)
        return 1
    except FileNotFoundError as e:
        print(f"File Not Found Error: {e}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as e:
        print(f"Snakemake execution failed: {e}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected Error: {e}", file=sys.stderr)
        print("\nTraceback:", file=sys.stderr)
        traceback.print_exc(file=sys.stderr)
        return 1
    
# Entry point for the script
if __name__ == "__main__":
    sys.exit(main())
    