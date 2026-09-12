"""Input validation functions for SILOSim"""
import os
import pandas as pd
import time
import re

class ValidationError(Exception):
    """Validation error"""
    pass

# Validate arguments for profiler
def validate_profiler(args):
    """Validate the arguments used in profiler function"""
    # Check if input_genome exists
    # For now I don't check if the input_genome is a valid fasta file for simplicity
    # But this will be checked during the actual analysis
    if not os.path.exists(args.input_genome):
        raise SystemExit(
            f"Input genome {args.input_genome} does not exist"
        )
    
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")

    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named silosim_profiler_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"silosim_profiler_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
    else:
        args.output = os.path.abspath(args.output)
    
    # Check ani_threshold
    if args.ani_threshold is None:
        print("ani_threshold not provided, using default 95")
        args.ani_threshold = 95
    elif args.ani_threshold < 90 or args.ani_threshold > 100:
        print("Warning: ani_threshold must be between 90 and 100 for reliable profiling.\n"
        "Using default threshold 95.")
        args.ani_threshold = 95

    # Check cov_cutoff
    if not args.cov_cutoff:
        print("cov_cutoff not provided, using default 0.7")
        args.cov_cutoff = 0.7
    elif args.cov_cutoff < 0 or args.cov_cutoff > 1:
        print("Warning: cov_cutoff must be between 0 and 1 for reliable profiling.\n"
        "Using default cutoff 0.7.")
        args.cov_cutoff = 0.7

    # Check min_ctg_len
    if not args.min_ctg_len:
        print("min_ctg_len not provided, using default 1000")
        args.min_ctg_len = 1000
    elif args.min_ctg_len <= 0:
        print("Warning: min_ctg_len must be a positive integer for reliable profiling.\n"
        "Using default value 1000.")
        args.min_ctg_len = 1000
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")

    # Check species
    if not args.species:
        args.species = "None"
    
    # Check whether local_query_dir is provided
    if args.local_query_dir:
        # args.whatsgnu_db_path and args.top_genomes will be automatically assigned to None if local_query_dir is provided
        args.whatsgnu_db_path = "None"
        args.top_genomes = 0
        # Normalize to an absolute path before any checks so error messages,
        # the config file, and the Snakemake rules all see the same resolved path
        args.local_query_dir = os.path.abspath(os.path.expanduser(args.local_query_dir))

        if not os.path.isdir(args.local_query_dir):
            raise ValidationError(f"Local query directory not found: {args.local_query_dir}")
        else:
            print(f"Local query directory found at: {args.local_query_dir}")
            # Check if the local_query_dir contains any fasta files
            fasta_files = [f for f in os.listdir(args.local_query_dir) if f.endswith(".fasta") or f.endswith(".fa") or f.endswith(".fna")]
            if not fasta_files:
                raise ValidationError(f"No FASTA files found in the local query directory: {args.local_query_dir}")
            else:
                print(f"Found {len(fasta_files)} FASTA files in the local query directory: {args.local_query_dir}")
                # Fail fast on stem collisions: the pipeline links every query genome into a
                # shared directory as <stem>.fna, so two files differing only by extension
                # would overwrite each other.
                stems = {}
                for f in fasta_files:
                    stems.setdefault(os.path.splitext(f)[0], []).append(f)

                collisions = {s: fs for s, fs in stems.items() if len(fs) > 1}
                if collisions:
                    detail = "\n".join(
                        f"  {s}: {', '.join(sorted(fs))}" for s, fs in sorted(collisions.items())
                    )
                    raise ValidationError(
                        "Multiple files in the local query directory share a genome name "
                        f"after the extension is removed:\n{detail}\n"
                        "Keep one file per genome, or rename them uniquely.")
                
                print(f"All FASTA files in the local query directory have unique stems: {args.local_query_dir}")

    else:
        # If args.local_query_dir is not provided, check args.whatsgnu_db_path and args.top_genomes
        # And assign "None" to args.local_query_dir
        args.local_query_dir = "None"
        # Check species
        # If local query not prodived, species must be specified
        if args.species == "None":
            raise ValidationError("Species must be specified when local query directory is not provided.")
        # Only check whatsgnu_db_path and top_genomes when local_query_dir is not provided
        # Check whatsgnu_db_path
        if not args.whatsgnu_db_path:
            print(f"WhatsGNU database path not provided. WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
            args.whatsgnu_db_path = "None"
        else:
            # Check if the file exists in the provided whatsgnu_db_path)
            if not os.path.exists(args.whatsgnu_db_path):
                print(f"WhatsGNU database file not found in the provided whatsgnu_db_path: {args.whatsgnu_db_path}")
                print(f"WhatsGNU database will be downloaded to {args.output}/whatsgnu/db")
                args.whatsgnu_db_path = "None"
            else:
                print(f"WhatsGNU database file found at: {args.whatsgnu_db_path}")
        # Check top_genomes
        if args.top_genomes < 100:
            raise ValidationError("top_genomes must be at least 100 for reliable profiling")

    # Check bakta_db_path and bakta_db_type
    if not args.bakta_db_path:
        print(f"Bakta database path not provided. Bakta database will be downloaded to {args.output}/bakta/db")
        args.bakta_db_path = "None"
    elif args.bakta_db_path:
        # Check if bakta.db exists in the provided bakta_db_path
        bakta_db_file = os.path.join(args.bakta_db_path, "bakta.db")
        if not os.path.exists(bakta_db_file):
            print(f"Bakta database file not found in the provided bakta_db_path: {bakta_db_file}")
            print(f"Bakta database will be downloaded to {args.output}/bakta/db")
            args.bakta_db_path = "None"
        else:
            print(f"Bakta database file found at: {bakta_db_file}")

    if not args.bakta_db_type:
        print("Bakta database type not provided, using default 'full'")
        args.bakta_db_type = "full"
    elif args.bakta_db_type not in {"full", "light"}:
        print("Unsupported Bakta database type, using default 'full'")
        args.bakta_db_type = "full"
# Validate arguments for simulator
def validate_simulator(args):
    """Validate the shared arguments used in simulator function"""
    # Check prefix
    if not args.prefix:
        print("Prefix not provided, using default prefix: current timestamp in the format of YYYY_MM_DD_HHMMSS")
        args.prefix = time.strftime("%Y_%m_%d_%H%M%S")
    
    # Check conda prefix
    if not args.conda_prefix:
        print("Conda environment path not provided, using default path: <OUTPUT>/conda_envs_<YYYY_MM_DD_HHMMSS>")
        args.conda_prefix = os.path.abspath(f"{args.output}/conda_envs_{args.prefix}")
    
    # Check seed
    if not args.seed:
        print("Seed not provided, using current date in the format of YYYYMMDD as seed")
        args.seed = int(pd.Timestamp.now().strftime("%Y%m%d"))
    elif args.seed < 0 or not isinstance(args.seed, int):
        raise ValidationError("Seed must be a non-negative integer")
    
    # Check output directory
    if not args.output:
        print(f"Output directory not provided, creating a directory named 'silosim_simulator_output_{args.prefix} under current working directory({os.getcwd()}) as output directory")
        args.output = os.path.join(os.getcwd(), f"silosim_simulator_output_{args.prefix}")
    if not os.path.exists(args.output):
        print(f"Output directory {args.output} does not exist, creating it.")
        os.makedirs(args.output)
        args.output = os.path.abspath(args.output)
    
    # Check ancestor_genome
    if not os.path.exists(args.ancestor):
        raise ValidationError(f"Ancestor genome file not found: {args.ancestor}")
    else:
        # Use absolute path
        args.ancestor = os.path.abspath(args.ancestor)
    
    # Check years
    if not args.years :
        raise ValidationError("Years must be provided")
    elif args.years <= 0 or not isinstance(args.years, int):
        raise ValidationError("Years must be a positive integer")

    if not args.taxa:
        raise ValidationError("Taxa (offspring) must be provided")
    elif args.taxa < 3 or not isinstance(args.taxa, int):
        raise ValidationError("Taxa (offspring) must be an integer no less than 3")
    
    # Check lambda_rate
    if args.lambda_rate is None:
        raise ValidationError("Lambda rate (--lambda_rate) must be provided")
    elif not isinstance(args.lambda_rate, (int, float)) or args.lambda_rate <= 0:
        raise ValidationError("Lambda rate must be a positive number")

    # Check mu_rate
    if args.mu_rate is None:
        raise ValidationError("Mu rate (--mu_rate) must be provided")
    elif not isinstance(args.mu_rate, (int, float)) or args.mu_rate < 0:
        raise ValidationError("Mu rate must be a non-negative number (0 for pure-birth Yule process)")

    # Check substitution model 
    if not args.substitution_model:
        print("Substitution model not provided, using default 'GTR'")
        args.substitution_model = "GTR"
    elif args.substitution_model not in {"JC69","K2P","K3P","GTR"}:
        print("Unsupported substitution model, using default 'GTR'")
        args.substitution_model = "GTR"

    # Check substitution model parameters for different models
    if args.substitution_model == "JC69":
        print("JC69 model selected, no additional parameters required.")
        args.model_parameters = "None"
    elif args.substitution_model == "K2P":
        print("K2P model selected, validating parameter...")
        if not args.model_parameters:
            raise ValidationError("K2P model requires transition/transversion ratio(kappa) as parameter. The parameter should be provided as a single float value greater than 0.")
        elif isinstance(float(args.model_parameters), float) or isinstance(int(args.model_parameters), int):
            if args.model_parameters <= 0:
                raise ValidationError("Transition/transversion ratio(kappa) must be a float value greater than 0 for K2P model")
            else:
                print(f"K2P model parameter (transition/transversion ratio, kappa) set to {args.model_parameters}")
        else:
            raise ValidationError("Transition/transversion ratio must be a float value greater than 0 for K2P model")
    elif args.substitution_model == "K3P":
        print("K3P model selected, validating parameters...")
        if not args.model_parameters:
            raise ValidationError("K3P model requires three parameters (alpha, gamma, beta) as parameters. The parameters should be provided as three float values greater than 0, separated by a comma (e.g., 0.4,0.6,0.8).")
        else:
            splited_params = args.model_parameters.split(",")
            if len(splited_params) != 3:
                raise ValidationError("K3P model requires three parameters (alpha, beta, gamma) as parameters. The parameters should be provided as three float values greater than 0, separated by commas (e.g., 0.4,0.6,0.8).")
            try:
                alpha = float(splited_params[0])
                beta = float(splited_params[1])
                gamma = float(splited_params[2])
                if alpha <= 0 or gamma <= 0 or beta <= 0:
                    raise ValidationError("All alpha, beta, and gamma parameters must be float values greater than 0 for K3P model")
                else:
                    print(f"K3P model parameters set to alpha: {alpha}, beta: {beta}, gamma: {gamma}")
            except ValueError:
                raise ValidationError("All alpha, beta, and gamma parameters must be float values greater than 0 for K3P model")
    elif args.substitution_model == "GTR":
        print("GTR model selected, validating parameters...")
        if not args.model_parameters:
            raise ValidationError("GTR model requires six rate parameters (a, b, c, d, e, f) as parameters. The parameters should be provided as six float values greater than 0, separated by commas (e.g., 0.1,0.2,0.3,0.15,0.15,0.1).")
        else:
            splited_params = args.model_parameters.split(",")
            if len(splited_params) != 6:
                raise ValidationError("GTR model requires six rate parameters (a, b, c, d, e, f) as parameters. The parameters should be provided as six float values greater than 0, separated by commas (e.g., 0.1,0.2,0.3,0.15,0.15,0.1).")
            try:
                rates = [float(param) for param in splited_params]
                if any(rate <= 0 for rate in rates):
                    raise ValidationError("All six rate parameters must be float values greater than 0 for GTR model")
                else:
                    print(f"GTR model parameters set to: a={rates[0]}, b={rates[1]}, c={rates[2]}, d={rates[3]}, e={rates[4]}, f={rates[5]}")
            except ValueError:
                raise ValidationError("All six rate parameters must be float values greater than 0 for GTR model")
        
    # Check mutation rate
    if not args.mutation_rate or not isinstance(args.mutation_rate, (int, float)) or args.mutation_rate <= 0:
        raise ValidationError("Mutation rate must be a positive float value")
    
    # Check if weighted mutation is on
    if not args.use_weighted_mutation:
        print("Whether using weighted mutation not provided, using default True")
        args.use_weighted_mutation = True
    elif args.use_weighted_mutation not in {True, False}:
        print("use_weighted_mutation must be either True or False. Using default True")
        args.use_weighted_mutation = True
    elif args.use_weighted_mutation:
        print("Weighted mutation enabled, validating weighted mutation input...")
        # Check weighted mutation file
        if not os.path.exists(args.weighted_mutation_file):
            raise ValidationError(f"Weighted mutation file not found: {args.weighted_mutation_file}")
    elif args.use_weighted_mutation == False:
        print("Weighted mutation disabled.")
        args.weighted_mutation_file = "None"
    
    # Check if recombination simulation is on
    if not args.use_recombination:
        print("Whether using recombination not provided, using default True")
        args.use_recombination = True
    elif args.use_recombination not in {True, False}:
        print("use_recombination must be either True or False. Using default True")
        args.use_recombination = True
    elif args.use_recombination:
        print("Recombination simulation enabled, validating recombination input...")
        # Check recombination rate 
        if not args.recombination_rate or not isinstance(args.recombination_rate, (int, float)) or args.recombination_rate <= 0:
            raise ValidationError("Recombination rate must be a positive float value. If no recombination simulation is desired, set --use_recombination to False")
        # Check mean recombination size
        if not args.mean_recombination_size or not isinstance(args.mean_recombination_size, float) or args.mean_recombination_size <= 0:
            raise ValidationError("Mean recombination size must be a positive float value")
        # Check minimal recombination size
        if not args.min_recombination_size or not isinstance(args.min_recombination_size, float) or args.min_recombination_size <= 0:
            raise ValidationError("Minimal recombination size must be a positive float value")
        elif args.min_recombination_size >= args.mean_recombination_size:
            raise ValidationError("Minimal recombination size must be smaller than mean recombination size")
        # Check nu
        if not args.nu or not isinstance(args.nu, (int, float)) or args.nu <= 0 or args.nu >=1:
            raise ValidationError("Nu must be a positive float value between 0 and 1")
        else: 
            print(f"Recombination parameters set to: recombination_rate={args.recombination_rate}, mean_recombination_size={args.mean_recombination_size}, minimal_recombination_size={args.min_recombination_size}, nu={args.nu}")
    elif args.use_recombination == False:
        print("Recombination simulation disabled.")
        args.recombination_rate = "None"
        args.mean_recombination_size = "None"
        args.min_recombination_size = "None"
        args.nu = "None"
    
    # Check if gene gain/loss simulation is on
    if not args.use_gain_loss:
        print("Whether using gene gain/loss not provided, using default True")
        args.use_gain_loss = True
    elif args.use_gain_loss not in {True, False}:
        print("use_gain_loss must be either True or False. Using default True")
        args.use_gain_loss = True
    elif args.use_gain_loss == False:
        print("Gene gain/loss simulation disabled.")
        args.gain_rate = "None"
        args.loss_rate = "None"
        args.bin = "None"
        args.mge_data = "None"
        args.position_coverage = "None"
        args.mge_fasta = "None"
        args.mge_entropy = "None"

    if args.use_gain_loss:
        print("Gene gain/loss simulation enabled, validating gene gain/loss input...")
        # Check gain rate
        if not args.gain_rate or not isinstance(args.gain_rate, (int, float)) or args.gain_rate < 0:
            raise ValidationError("Gain rate must be a non-negative float value")
        # Check loss rate
        if not args.loss_rate or not isinstance(args.loss_rate, (int, float)) or args.loss_rate < 0:
            raise ValidationError("Loss rate must be a non-negative float value")
        # Check chromosome bin file
        if not os.path.exists(args.bin):
            raise ValidationError(f"Chromosome bin file not found: {args.bin}")
        # Check MGE data file
        if not os.path.exists(args.mge_data):
            raise ValidationError(f"MGE data file not found: {args.mge_data}")
        # Check MGE fasta file
        if not os.path.exists(args.mge_fasta):
            raise ValidationError(f"MGE fasta file not found: {args.mge_fasta}")
        # Check MGE entropy file
        if not os.path.exists(args.mge_entropy):
            raise ValidationError(f"The directory of MGE entropy file not found: {args.mge_entropy}")
        # Check position coverage file
        if not os.path.exists(args.position_coverage):
            raise ValidationError(f"Position coverage file not found: {args.position_coverage}")
        else:
            print(f"Gene gain/loss parameters set to: gain_rate={args.gain_rate}, loss_rate={args.loss_rate}, bin={args.bin}, mge_data={args.mge_data}, mge_fasta={args.mge_fasta}, mge_entropy={args.mge_entropy}, position_coverage={args.position_coverage}")
        
    elif args.use_gain_loss == False:
        print("Gene gain/loss simulation disabled.")
        args.gain_rate = "None"
        args.loss_rate = "None"
        args.bin = "None"
        args.mge_fasta = "None"
        args.mge_entropy = "None"
        args.position_coverage = "None"
