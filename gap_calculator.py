import re
import numpy as np

# === CONFIGURATION ===
meme_file_path = "/cluster/home/zacharis/output/meme/shorter_motifs/above3/sigma70_zoops_4_12_above3_withouttss/meme.txt"  # ← Replace with your actual file
group1_motifs = ["CTTG", "GTTG"]  # -35 elements
group2_motifs = ["ATAT"]   # -10 elements

# === Step 1: Parse the MEME file ===
def parse_meme_sites(file_path):
    motif_name = None
    all_sites = []

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()

            # Detect the start of a new motif block
            motif_header_match = re.match(r"Motif (\S+) MEME", line)
            if motif_header_match:
                motif_name = motif_header_match.group(1).strip().upper()
                continue

            # Skip lines that are not part of the motif sites section
            if re.match(r"[-]+", line) or line.startswith("Sequence name") or not line:
                continue

            # Extract site information (sequence name, start, site match)
            parts = line.split()
            if len(parts) < 5:
                continue

            try:
                seq_name = parts[0]
                start = int(parts[1])
                site_seq = parts[4].lower()
                end = start + len(site_seq) - 1

                all_sites.append({
                    "sequence": seq_name,
                    "motif": motif_name,
                    "start": start,
                    "end": end,
                    "site": site_seq
                })
            except ValueError:
                continue  # Skip lines that don’t contain motif data

    return all_sites

# === Step 2: Organize occurrences by sequence ===
def group_sites_by_sequence(sites):
    seq_dict = {}
    for site in sites:
        seq = site["sequence"]
        if seq not in seq_dict:
            seq_dict[seq] = []
        seq_dict[seq].append(site)
    return seq_dict

# === Step 3: Compute motif gaps ===
def compute_motif_gaps(seq_dict, group1, group2):
    group1 = [m.upper() for m in group1]
    group2 = [m.upper() for m in group2]
    gaps = []

    for seq_name, occurrences in seq_dict.items():
        g1_sites = [s for s in occurrences if s["motif"] in group1]
        g2_sites = [s for s in occurrences if s["motif"] in group2]

        if not g1_sites or not g2_sites:
            continue  # Need at least one from each group

        for g1 in g1_sites:
            for g2 in g2_sites:
                if g2["start"] > g1["end"]:  # Ensure order
                    gap = g2["start"] - g1["end"]
                    gaps.append({
                        "sequence": seq_name,
                        "motif1": g1["motif"],
                        "motif2": g2["motif"],
                        "start1": g1["start"],
                        "end1": g1["end"],
                        "start2": g2["start"],
                        "gap": gap
                    })

    return gaps

# === Step 4: Print stats ===
def print_gap_stats(gap_list):
    if not gap_list:
        print("⚠️ No valid motif pairs found.")
        return

    gap_values = [g["gap"] for g in gap_list]
    gap_values.sort()

    print(f"✅ Valid motif pairs found: {len(gap_values)}")
    print(f"Mean gap: {np.mean(gap_values):.2f}")
    print(f"Median gap: {np.median(gap_values):.2f}")
    print(f"1st Quartile (25%): {np.percentile(gap_values, 25):.2f}")
    print(f"3rd Quartile (75%): {np.percentile(gap_values, 75):.2f}")
    print(f"5th Percentile: {np.percentile(gap_values, 5):.2f}")
    print(f"95th Percentile: {np.percentile(gap_values, 95):.2f}")
    print(f"Shortest gap: {min(gap_values)}")
    print(f"Longest gap: {max(gap_values)}")


# === RUN ===
all_sites = parse_meme_sites(meme_file_path)
seq_dict = group_sites_by_sequence(all_sites)
gap_list = compute_motif_gaps(seq_dict, group1_motifs, group2_motifs)
print_gap_stats(gap_list)
