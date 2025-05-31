import csv
import os
from collections import defaultdict

def process_sigma_csv(input_csv, output_dir):
    sigma_dict = defaultdict(list)

    os.makedirs(output_dir, exist_ok=True)

    with open(input_csv, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:


            sequence = row['sequence'][-50:-2]  # change length of promoter to include
            sigma_factors = row['sigmafactor'].strip().split(';')
            expression_levels = row['expression_level_log'].split(';')

            if not sigma_factors or len(sigma_factors) != len(expression_levels):
                continue  

            for sigma, expression_level in zip(sigma_factors, expression_levels):
                sigma = sigma.strip()
                expression_level = float(expression_level.strip())

                if 3 < expression_level <= 5: # change expression level cut-off value
                    sigma_dict[sigma].append((expression_level, sequence))

    for sigma, sequences in sigma_dict.items():
        fasta_filename = os.path.join(output_dir, f"{sigma}.fasta")
        with open(fasta_filename, 'w', buffering=1) as fasta_file:
            for idx, (expression_level, seq) in enumerate(sequences, start=1):
                fasta_file.write(f">seq_{idx}_{sigma}_{expression_level:.2f}\n")
                fasta_file.write(f"{seq}\n")


    print("FASTA files created successfully in", output_dir)

# Example usage
process_sigma_csv("/cluster/home/zacharis/data/sigma_dataset_complete_log_transformed.csv", 
                  "/cluster/home/zacharis/data/without_tss/between_3and5")
