import os
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import logomaker
# Function to extract information content (ic) from the meme.html file
def extract_information_content(html_file_path):
    with open(html_file_path, 'r') as file:
        content = file.read()

    # Regular expression to extract the information content (ic)
    ic_match = re.search(r'"ic":\s*([\d.]+)', content)
    if ic_match:
        ic_value = float(ic_match.group(1))
    else:
        ic_value = None  # or set a default value if not found
    
    return ic_value

# Function to parse MEME file and extract the letter probability matrix, additional info, and motif start sites
def parse_meme_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()

    # Extract the letter probability matrix (probabilities for A, C, G, T)
    matrix_match = re.search(r"letter-probability\s+matrix:.*?nsites=.*?\n((?:[0-9.\s]+\n)+)", content, re.S)
    if not matrix_match:
        raise ValueError(f"No letter probability matrix found in the MEME file: {file_path}")

    matrix_text = matrix_match.group(1).strip()
    matrix_rows = [list(map(float, line.split())) for line in matrix_text.split("\n")]
    matrix_np = np.array(matrix_rows)

    # Convert matrix into a Pandas DataFrame with column names for A, C, G, T
    df = pd.DataFrame(matrix_np, columns=["A", "C", "G", "T"])

    # Extract motif, width, sites, p-value, and N (total sequence number)
    motif_match = re.search(r"MOTIF\s+(\S+)", content)
    width_match = re.search(r"width\s+=\s+(\d+)", content)
    sites_match = re.search(r"sites\s+=\s+(\d+)", content)
    p_value_match = re.search(r"p-value\s+=\s+([\d.e-]+)", content)
    n_match = re.search(r"N=\s+(\d+)", content)

    motif = motif_match.group(1) if motif_match else "Unknown"
    width = int(width_match.group(1)) if width_match else 0
    sites = int(sites_match.group(1)) if sites_match else 0
    p_value = p_value_match.group(1) if p_value_match else "Unknown"
    total_sequences = int(n_match.group(1)) if n_match else 0

    # Calculate the percentage of sites
    percentage_sites = (sites / total_sequences * 100) if total_sequences else 0

    # ----- New Code: Extract motif start sites and calculate median start site -----
    # Look for the section containing motif sites sorted by position
    median_start = None
    start_sites = []
    lines = content.splitlines()
    table_found = False
    for idx, line in enumerate(lines):
        if "sites sorted by position" in line:
            table_found = True
            # Advance until the header "Sequence name" is found
            while idx < len(lines) and "Sequence name" not in lines[idx]:
                idx += 1
            # Skip the header and dashed line(s)
            idx += 1
            while idx < len(lines) and set(lines[idx].strip()) == set("-"):
                idx += 1
            # Now parse the data lines until an empty line or a non-data line is reached
            while idx < len(lines) and lines[idx].strip():
                parts = lines[idx].split()
                if len(parts) >= 2:
                    try:
                        start_val = int(parts[1])
                        start_sites.append(start_val)
                    except ValueError:
                        pass
                idx += 1
            break

    if start_sites:
        median_start = int(np.median(start_sites))
    else:
        # If no sites are found, default to 0 (or handle as needed)
        median_start = 0
    # -----------------------------------------------------------------------------

    return df, motif, width, sites, p_value, total_sequences, percentage_sites, median_start

# Function to generate and save sequence logo from the probability matrix, with adjusted x-axis numbering


# Function to classify the dataset based on the folder name
def get_dataset_description(folder_name):
    # Regular expressions for different expressions
    if "above2" in folder_name:
        return "Expression level: above 2"
    elif "above3" in folder_name:
        return "Expression level: above 3"
    elif "above35" in folder_name:
        return "Expression level: above 35"
    elif "above4" in folder_name:
        return "Expression level: above 4"
    elif "above5" in folder_name:
        return "Expression level: above 5"
    elif "above6" in folder_name:
        return "Expression level: above 6"
    
    elif "between3and4" in folder_name:
        return "Expression level: between 3 and 4"
    elif "between28and3" in folder_name:
        return "Expression level: between 2.8 and 3"
    elif "between3and32" in folder_name:
        return "Expression level: between 3 and 3.2"
    elif "between32and35" in folder_name:
        return "Expression level: between 3.2 and 3.5"
    elif "between35and4" in folder_name:
        return "Expression level: between 3.5 and 4"
    
    elif "fulldataset" in folder_name:
        return "Expression level: entire dataset"
    else:
        return "Expression level of sequences: unknown"

# Modify the generate_sequence_logo function to include dataset info
# Modify the generate_sequence_logo function to include dataset info
# Modify the generate_sequence_logo function to include information content
def generate_sequence_logo(df, output_dir, output_file, organism_name, sigma_id, motif, width, sites, p_value, percentage_sites, median_start, dataset_description, ic_value):
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, output_file)

    # Font code
    regular_font_path = "/cluster/home/zacharis/fonts/Montserrat-Medium.ttf"
    italic_font_path = "/cluster/home/zacharis/fonts/Montserrat-MediumItalic.ttf"

    fm.fontManager.addfont(regular_font_path)
    fm.fontManager.addfont(italic_font_path)

    regular_font = fm.FontProperties(fname=regular_font_path)
    italic_font = fm.FontProperties(fname=italic_font_path)

    plt.rcParams["font.family"] = regular_font.get_name()

    # Convert probability matrix to information content (bits)
    df_bits = logomaker.transform_matrix(df, from_type="probability", to_type="information")

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Create sequence logo
    logo = logomaker.Logo(df_bits, ax=ax, shade_below=0.5, fade_below=0.5)

    # Adjust x-axis numbering based on the median start site
    x_ticks = list(range(width))
    x_tick_labels = [
        str(median_start - 49 + i) if (median_start - 49 + i) % 2 == 0 else ""
        for i in range(width)
    ]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.set_xlabel("Median position relative to tss", fontsize=14, labelpad=5)

    # Set y-axis label
    ax.set_ylabel("Bits", fontsize=14)
    ax.set_yticks([1, 2])
    ax.set_yticklabels(["1", "2"])

    # Set the overall title with dataset description below it
    fig.text(0.5, 1.05, f"${organism_name}$ $\\mathrm{{\\sigma}}^{{\\mathrm{{{sigma_id}}}}}$", 
             ha="center", va="bottom", fontsize=20, fontproperties=italic_font)

    # Add dataset description as a subtitle right below the title
    fig.text(0.5, 0.95, dataset_description, ha="center", va="bottom", fontsize=14, fontproperties=regular_font)

    # Prepare the information to display below the plot
    text_box = f"p-value = {p_value}\n"
    text_box += f"sites = {sites} ({percentage_sites:.2f}%)\n"
    
    if ic_value is not None:
        text_box += f"Information Content (ic) = {ic_value:.2f} bits\n"
    else:
        text_box += "Information Content (ic) = N/A\n"

    # Add the text box under the plot without the outline
    ax.text(0.5, -0.2, text_box, ha="center", va="top", fontsize=12, transform=ax.transAxes, 
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", boxstyle="round,pad=0.5"))

    # Save the sequence logo
    try:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved logo: {output_path}")
    except Exception as e:
        print(f"Failed to save logo: {e}")



# Directory containing MEME output folders
base_directory = "/cluster/home/zacharis/output/meme/larger_motifs/sigma32p/"
output_directory = "/cluster/home/zacharis/output/final_plots/large_motifs/sigma32p"

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




# Iterate over each folder in the base directory
for folder_name in os.listdir(base_directory):
    folder_path = os.path.join(base_directory, folder_name)

    if os.path.isdir(folder_path):
        meme_file_path = os.path.join(folder_path, "meme.txt")

        if os.path.exists(meme_file_path):
            # Extract sigma factor from the folder name (before the first underscore)
            sigma_factor = folder_name.split('_')[0]

            # Map the sigma factor to a human-readable label
            organism_name, sigma_id = sigma_label_map.get(sigma_factor.lower(), ("Unknown bacterium", sigma_factor.upper()))

            # Get the dataset description based on the folder name
            dataset_description = get_dataset_description(folder_name)

            # Construct output file name using the folder name
            output_file_name = f"{folder_name}_motif_logo_bits.png"

            
            # Path to the meme.html file (same folder as meme.txt)
            html_file_path = os.path.join(folder_path, "meme.html")

            # Initialize ic_value as None by default
            ic_value = None

            # Extract information content (ic) only if the HTML file exists
            if os.path.exists(html_file_path):
                ic_value = extract_information_content(html_file_path)
            else:
                print(f"Warning: meme.html not found in {folder_path}. IC will be set to N/A.")

            try:
                # Parse MEME file to extract the probability matrix, additional info, and median start site
                (motif_df, motif, width, sites, p_value, total_sequences, 
                 percentage_sites, median_start) = parse_meme_file(meme_file_path)

                # Generate sequence logo with adjusted x-axis
                print(f"Generating logo for {folder_name}")
                generate_sequence_logo(
                    motif_df, output_directory, output_file_name, organism_name, sigma_id, motif,
                    width, sites, p_value, percentage_sites, median_start, dataset_description, ic_value
                )

            except ValueError as e:
                print(f"Skipping file {meme_file_path} due to error: {e}")
        else:
            print(f"MEME file not found in {folder_path}")
