import random
import itertools
import csv
from collections import defaultdict
import statistics
import math
from scipy.stats import norm



def run_simulation(subset_sizes):
    dataset = list(range(1, dataset_size + 1))
    subsets = {
        name: set(random.sample(dataset, size))
        for name, size in subset_sizes.items()
    }

    def compute_intersection(subset_names):
        inter_set = subsets[subset_names[0]]
        for name in subset_names[1:]:
            inter_set = inter_set.intersection(subsets[name])
        return len(inter_set)

    intersections = {}
    all_names = list(subsets.keys())
    for r in range(2, len(all_names) + 1):
        for combo in itertools.combinations(all_names, r):
            intersections[combo] = compute_intersection(combo)
    return intersections

def load_real_gene_sets(file_paths):
    gene_sets = {}
    for name, path in file_paths.items():
        with open(path, 'r') as f:
            genes = {line.strip() for line in f if line.strip()}
            gene_sets[name] = genes

    def compute_intersection(subset_names):
        inter_set = gene_sets[subset_names[0]]
        for name in subset_names[1:]:
            inter_set = inter_set.intersection(gene_sets[name])
        return len(inter_set)

    intersections = {}
    all_names = list(gene_sets.keys())
    for r in range(2, len(all_names) + 1):
        for combo in itertools.combinations(all_names, r):
            intersections[combo] = compute_intersection(combo)

    return intersections, gene_sets

def write_results_to_csv(output_path, real_intersections, all_intersections):
    """
    Writes the real vs simulated intersection data including z-scores and p-values to CSV.
    """
    with open(output_path, mode="w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Combination", "real size", "simulated mean", "z-score", "p-value"])

        for combo, real_size in sorted(real_intersections.items()):
            sorted_combo = tuple(sorted(combo))
            sim_sizes = all_intersections.get(combo) or all_intersections.get(sorted_combo)

            if sim_sizes:
                mean_val = statistics.mean(sim_sizes)
                stdev_val = statistics.stdev(sim_sizes) if len(sim_sizes) > 1 else 0
                if stdev_val > 0:
                    z_score = (real_size - mean_val) / stdev_val
                    p_value = 2 * (1 - norm.cdf(abs(z_score)))
                else:
                    z_score = float('nan')
                    p_value = float('nan')
                combo_str = ", ".join(combo)
                writer.writerow([combo_str, real_size, f"{mean_val:.2f}", f"{z_score:.5f}", f"{p_value:.10f}"])
            else:
                combo_str = ", ".join(combo)
                writer.writerow([combo_str, real_size, "N/A", "N/A", "N/A"])

# Define simulation parameters
num_simulations = 100000
dataset_size = 4191

# Define file paths
file_paths = {
    "sigma70p": "/cluster/home/zacharis/output/fimo/new_strategy/sigma70p_long_motifs_locus_ids.txt",
    "sigma54p": "/cluster/home/zacharis/output/fimo/new_strategy/sigma54p_combined_motifs_locus_ids.txt",
    "sigma38p": "/cluster/home/zacharis/output/fimo/new_strategy/sigma38p_long_motifs_locus_ids.txt",
    "sigma32p": "/cluster/home/zacharis/output/fimo/new_strategy/sigma32p_long_motifs_locus_ids.txt",
    "sigma28p": "/cluster/home/zacharis/output/fimo/new_strategy/sigma28p_combined_motifs_locus_ids.txt"
}

# Calculate real intersections
real_intersections, gene_sets = load_real_gene_sets(file_paths)

# Derive subset sizes from actual gene sets
subset_sizes = {name: len(genes) for name, genes in gene_sets.items()}

# Run simulations
all_intersections = defaultdict(list)
for _ in range(num_simulations):
    sim_result = run_simulation(subset_sizes)
    for combo, inter_size in sim_result.items():
        all_intersections[combo].append(inter_size)


# Write results
output_file = "/cluster/home/zacharis/output/fimo/new_strategy/monte_carlo_putida.csv"
write_results_to_csv(output_file, real_intersections, all_intersections)