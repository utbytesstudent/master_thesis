import os
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import logomaker

# Function to parse the XML and extract motifs
def parse_xml_motifs(xml_file_path):
    tree = ET.parse(xml_file_path)
    root = tree.getroot()
    motifs_data = []

    max_sites_elem = root.find('model/maxsites')
    max_sites = int(max_sites_elem.text) if max_sites_elem is not None else None

    motifs_elem = root.find('motifs')
    if motifs_elem is not None:
        for motif_elem in motifs_elem.findall('motif'):
            motif_id = motif_elem.get('id')
            motif_name = motif_elem.get('name')
            width = int(motif_elem.get('width'))
            sites = int(motif_elem.get('sites'))
            ic = float(motif_elem.get('ic'))
            p_value = motif_elem.get('p_value')

            probabilities_elem = motif_elem.find('probabilities')
            if probabilities_elem is not None:
                prob_matrix = []
                alphabet_matrix_elem = probabilities_elem.find('alphabet_matrix')
                for row_elem in alphabet_matrix_elem.findall('alphabet_array'):
                    row = {}
                    for value_elem in row_elem.findall('value'):
                        letter = value_elem.get('letter_id')
                        value = float(value_elem.text)
                        row[letter] = value
                    prob_matrix.append(row)
                prob_matrix_df = pd.DataFrame(prob_matrix, columns=["A", "C", "G", "T"])
            else:
                prob_matrix_df = None
            
            contributing_sites = []
            contributing_sites_elem = motif_elem.find('contributing_sites')
            if contributing_sites_elem is not None:
                for site_elem in contributing_sites_elem.findall('contributing_site'):
                    position = int(site_elem.get('position'))
                    contributing_sites.append(position)

            motif_data = {
                "motif_id": motif_id,
                "motif_name": motif_name,
                "width": width,
                "sites": sites,
                "ic": ic,
                "p_value": p_value,
                "probability_matrix": prob_matrix_df,
                "contributing_sites": contributing_sites
            }
            motifs_data.append(motif_data)

    return motifs_data, max_sites

# Function to generate sequence logo and save it as a file
def generate_sequence_logo(df, output_dir, output_file, organism_name, sigma_id, motif,
                           width, sites, p_value, percentage_sites, median_start,
                           dataset_description, ic_value, start_sites):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    # Font paths
    regular_font_path = "/cluster/home/zacharis/fonts/Montserrat-Medium.ttf"
    italic_font_path = "/cluster/home/zacharis/fonts/Montserrat-MediumItalic.ttf"
    fm.fontManager.addfont(regular_font_path)
    fm.fontManager.addfont(italic_font_path)
    regular_font = fm.FontProperties(fname=regular_font_path)
    italic_font = fm.FontProperties(fname=italic_font_path)
    plt.rcParams["font.family"] = regular_font.get_name()

    # Convert probability matrix to information content matrix
    df_bits = logomaker.transform_matrix(df, from_type="probability", to_type="information")
    fig, axes = plt.subplots(1, 2, figsize=(14, 4), gridspec_kw={'width_ratios': [2, 3]})
    ax_logo, ax_hist = axes

    # Sequence Logo Plot
    ax_logo.spines['top'].set_visible(False)
    ax_logo.spines['right'].set_visible(False)
    logo = logomaker.Logo(df_bits, ax=ax_logo, shade_below=0.5, fade_below=0.5)

    ax_logo.set_xticks([])  # Remove x-axis ticks
    ax_logo.set_xticklabels([])  # Remove x-axis tick labels
    ax_logo.set_ylabel("Bits", fontsize=18)
    ax_logo.set_yticks([1, 2])
    ax_logo.set_yticklabels(["1", "2"], fontsize=16)

    # Histogram of motif positions
    bins = np.arange(-50, 1, 1)
    ax_hist.hist(start_sites, bins=bins, edgecolor='black', alpha=0.7)
    ax_hist.set_xlabel("Motif position relative TSS", fontsize=22)
    ax_hist.set_ylabel("Occurrences", fontsize=16)
    ax_hist.set_xticks(np.arange(-50, 1, 10))
    ax_hist.set_xticklabels([str(x) for x in np.arange(-50, 1, 10)], fontsize=14)
    max_count = max(np.histogram(start_sites, bins=np.arange(-50, 1, 1))[0])
    yticks = np.linspace(0, max_count, 4, dtype=int)
    ax_hist.set_yticks(yticks)
    ax_hist.set_yticklabels([str(y) for y in yticks], fontsize=14)


    try:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved logo and histogram: {output_path}")
    except Exception as e:
        print(f"Failed to save logo and histogram: {e}")
    finally:
        plt.close(fig)

# Main function to process all XML files in the folder
def process_all_files_in_folder(folder_path, output_directory, dataset_description, sigma_label_map):
    for subdir, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith("meme.xml"):
                folder_name = os.path.basename(subdir)
                sigma_id = folder_name.split('_')[0].lower()
                organism_name, sigma_label = sigma_label_map.get(sigma_id, ("Unknown", "Unknown"))
                
                xml_file_path = os.path.join(subdir, file)
                print(f"Processing {xml_file_path} for {organism_name} sigma {sigma_label}")
                motifs, max_sites = parse_xml_motifs(xml_file_path)
                if max_sites is None or max_sites == 0:
                    print(f"Warning: max_sites not found or zero in {xml_file_path}")
                    continue

                for idx, motif_data in enumerate(motifs, start=1):
                    motif_name = motif_data["motif_name"]
                    motif_id = motif_data["motif_id"]
                    width = motif_data["width"]
                    sites = motif_data["sites"]
                    ic_value = motif_data["ic"]
                    p_value = motif_data["p_value"]
                    probability_matrix = motif_data["probability_matrix"]
                    contributing_sites = motif_data["contributing_sites"]

                    percentage_sites = (sites / max_sites * 100)
                    median_start = int(np.median(contributing_sites)) if contributing_sites else 0
                    center_sites = [start + width // 2 - 50 for start in contributing_sites]

                    output_file_name = f"{organism_name}_sigma{sigma_label}_{motif_name}_motif_{idx}_logo_bits.png"
                    print(f"Generating logo for motif {idx} ({motif_name})")

                    if probability_matrix is not None:
                        generate_sequence_logo(
                            probability_matrix, output_directory, output_file_name, organism_name, sigma_label, motif_name,
                            width, sites, p_value, percentage_sites, median_start, dataset_description, ic_value, center_sites
                        )

# Define sigma label mapping
sigma_label_map = {
    "siga": ("B. subtilis", "a"),
    "sigb": ("B. subtilis", "b"),
    "sigd": ("B. subtilis", "d"),
    "sige": ("B. subtilis", "e"),
    "sigf": ("B. subtilis", "f"),
    "sigh": ("B. subtilis", "h"),
    "sign": ("B. subtilis", "n"),
    "sigma70": ("E. coli", "70"),
    "sigma38": ("E. coli", "38"),
    "sigma32": ("E. coli", "32"),
    "sigma28": ("E. coli", "28"),
    "sigma70p": ("P. putida", "70p"),
    "sigma54p": ("P. putida", "54p"),
    "sigma38p": ("P. putida", "38p"),
    "sigma32p": ("P. putida", "32p"),
    "sigma28p": ("P. putida", "28p"),    
}

# Set paths
folder_path = "/cluster/home/zacharis/output/meme/shorter_motifs/above3/"
output_directory = "/cluster/home/zacharis/output/final_plots/clean/short_motifs/above3"
dataset_description = "Expression level: above 3"

# Run processing
process_all_files_in_folder(folder_path, output_directory, dataset_description, sigma_label_map)
