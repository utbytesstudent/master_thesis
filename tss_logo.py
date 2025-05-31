import os
from Bio import SeqIO
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm


def extract_nucleotides(fasta_file, num_nucleotides=10, position='end'):
    """
    Extracts nucleotides from each sequence in the FASTA file.
    The extraction can be from the start, end, or a specific position.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < num_nucleotides:
            continue
        if position == 'end':
            sequences.append(seq[-num_nucleotides:])
        elif position == 'start':
            sequences.append(seq[:num_nucleotides])
        elif isinstance(position, int) and 0 <= position <= len(seq) - num_nucleotides:
            sequences.append(seq[position:position + num_nucleotides])
    return sequences


def create_nucleotide_frequency_matrix(sequences):
    """
    Creates a DataFrame with the frequency of each nucleotide (A, C, G, T)
    at every position of the provided sequences.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    num_positions = len(sequences[0]) if sequences else 0
    position_counts = {i: {n: 0 for n in nucleotides} for i in range(1, num_positions + 1)}

    for seq in sequences:
        for i, nucleotide in enumerate(seq):
            if nucleotide in nucleotides:
                position_counts[i + 1][nucleotide] += 1

    num_sequences = len(sequences)
    if num_sequences:
        position_frequencies = {
            pos: {n: count / num_sequences for n, count in counts.items()}
            for pos, counts in position_counts.items()
        }
    else:
        position_frequencies = {}

    matrix = pd.DataFrame(position_frequencies).T
    return matrix


def convert_to_bits(matrix):
    """
    Converts the frequency matrix to information content (bits).
    """
    entropy_max = np.log2(4)  # Maximum entropy for 4 nucleotides

    def calc_information_content(column):
        # Add a small constant to avoid log(0)
        entropy = -np.nansum(column * np.log2(column + 1e-9))
        return entropy_max - entropy

    info_content = matrix.apply(calc_information_content, axis=1)
    return matrix.multiply(info_content, axis=0)


def extract_title(fasta_file):
    """
    Extracts a descriptive title based on the FASTA filename.
    If a match is found in the sigma_label_map, a formatted title is returned.
    Otherwise, a default title is generated.
    """
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
    basename = os.path.basename(fasta_file).replace(".fasta", "").replace(".fa", "")
    if basename in sigma_label_map:
        organism, sigma_number = sigma_label_map[basename]
        title = f"Transcription start site\n{organism} $\sigma^{{{sigma_number}}}$"
    else:
        title = f"Transcription start site\n{basename}"
    return title


def generate_sequence_logo(matrix, output_file, fasta_file):
    """
    Generates and saves a sequence logo plot using the provided matrix.
    The plot includes customized fonts, axis labels, and a title based on the FASTA file.
    """
    # Setup font paths and load fonts
    regular_font_path = "/cluster/home/zacharis/fonts/Montserrat-Medium.ttf"
    italic_font_path = "/cluster/home/zacharis/fonts/Montserrat-MediumItalic.ttf"
    fm.fontManager.addfont(regular_font_path)
    fm.fontManager.addfont(italic_font_path)
    regular_font = fm.FontProperties(fname=regular_font_path)

    plt.rcParams["font.family"] = regular_font.get_name()
    plt.rcParams["axes.unicode_minus"] = False

    font_size = 18

    # Convert frequencies to information content (bits)
    bit_matrix = convert_to_bits(matrix)

    # Create the logo plot
    fig, ax = plt.subplots(figsize=(3, 4))
    logo = logomaker.Logo(bit_matrix, ax=ax)

    # Set font properties for each letter in the logo
    for text in logo.ax.texts:
        text.set_fontproperties(regular_font)
        text.set_fontsize(font_size)

    # Set custom tick labels
    num_positions = len(bit_matrix)
    tick_positions = [i + 1 for i in range(num_positions)]
    tick_labels = [str(i - (num_positions - 1)) if i - (num_positions - 1) != 0 else "TSS" for i in range(num_positions)]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontproperties=regular_font, fontsize=font_size)
    ax.set_ylabel("Bits", fontproperties=regular_font, fontsize=font_size)
    # Example: Set y-axis ticks at 0.0, 0.5, 1.0 and label them
    y_ticks = [0.0, 0.5, 1.0]
    y_labels = ['0', '0.5', '1.0']
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_labels, fontproperties=regular_font, fontsize=font_size)

    ax.set_ylim(0, 1.2)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(output_file, dpi=300)
    plt.close(fig)


def process_single_fasta(fasta_file, output_dir, num_nucleotides, position):
    """
    Processes a single FASTA file to extract sequences, create the frequency matrix,
    and generate a sequence logo.
    """
    sequences = extract_nucleotides(fasta_file, num_nucleotides, position)
    if not sequences:
        print(f"Skipping {fasta_file}: no valid sequences found.")
        return

    matrix = create_nucleotide_frequency_matrix(sequences)
    if matrix.empty or matrix.sum().sum() == 0:
        print(f"No variation in sequences for {fasta_file}; logo not generated.")
    else:
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_file = os.path.join(output_dir, f"{base_name}_logo_bits.png")
        generate_sequence_logo(matrix, output_file, fasta_file)
        print(f"Logo saved to {output_file}")


def process_fasta_files(input_dir, output_dir, num_nucleotides=10, position='end'):
    """
    Processes all FASTA files in the given directory.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = [f for f in os.listdir(input_dir)
                   if f.endswith('.fasta') or f.endswith('.fa')]
    if not fasta_files:
        print("No FASTA files found in the provided directory.")
        return

    for fasta_file in fasta_files:
        fasta_file_path = os.path.join(input_dir, fasta_file)
        process_single_fasta(fasta_file_path, output_dir, num_nucleotides, position)


# Example usage
if __name__ == "__main__":
    output_dir = '/cluster/home/zacharis/output/logos/tss_3nt/'
    input_dir = '/cluster/home/zacharis/data/sequences_full_dataset/'  # Directory with all FASTA files
    num_nucleotides = 3  # Change as needed
    position = 'end'  # Options: 'start', 'end', or an integer for a specific position

    process_fasta_files(input_dir=input_dir, output_dir=output_dir,
                        num_nucleotides=num_nucleotides, position=position)
