import os
import yaml
# profiler function config creators
def profiler_config(args):
    """Create config file in yaml format for profiler function"""
    config = {
        'prefix': args.prefix,
        'input_genome': args.input_genome,
        'output': args.output,
        'species': args.species,
        'local_query_dir': args.local_query_dir,
        'top_genomes': args.top_genomes,
        'min_ctg_len': args.min_ctg_len,
        'min_mge_len': args.min_mge_len,
        'cov_cutoff': args.cov_cutoff,
        'ani_threshold': args.ani_threshold,
        'bakta_db_type': args.bakta_db_type,
        'bakta_db_path': args.bakta_db_path,
        'whatsgnu_db_path': args.whatsgnu_db_path,
        'threads': args.threads
    }
    # If thread is 1, give warning about long runtime for large datasets
    if args.threads == 1:
        print("Warning: Using a single thread may result in long runtimes. Consider using multiple threads for better performance.")

    if not os.path.exists(os.path.join(args.output, "config")):
        os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
# simulator function config creators
def simulator_config(args):
    """Create config file in yaml format for simulator function"""
    config = {
        'prefix': args.prefix,
        'ancestor': args.ancestor,
        'years': args.years,
        'taxa': args.taxa,
        'lambda': args.lambda_rate,
        'mu': args.mu_rate,
        'substitution_model': args.substitution_model,
        'model_parameters': args.model_parameters,
        'mutation_rate': args.mutation_rate,
        'use_weighted_mutation': args.use_weighted_mutation,
        'weighted_mutation_file': args.weighted_mutation_file,
        'use_gain_loss': args.use_gain_loss,
        'gain_rate': args.gain_rate,
        'loss_rate': args.loss_rate,
        'bin': args.bin,
        'position_coverage': args.position_coverage,
        'mge_data': args.mge_data,
        'mge_fasta': args.mge_fasta,
        'mge_entropy': args.mge_entropy,
        'use_recombination': args.use_recombination,
        'recombination_rate': args.recombination_rate,
        'min_recombination_size': args.min_recombination_size,
        'mean_recombination_size': args.mean_recombination_size,
        'nu': args.nu,
        'export_per_event_genomes': args.export_per_event_genomes,
        'output': args.output,
        'seed': args.seed,
        'threads': args.threads
    }
    
    if not os.path.exists(os.path.join(args.output, "config")):
       os.makedirs(os.path.join(args.output, "config"))

    with open(os.path.join(args.output, "config", f"config_{args.prefix}.yaml"), 'w') as f:
       yaml.dump(config, f, default_flow_style=False)

    return os.path.join(args.output, "config", f"config_{args.prefix}.yaml")
