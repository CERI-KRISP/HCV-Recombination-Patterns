from Bio import AlignIO
import csv

# Load alignments
stk_alignment = AlignIO.read("F5_hcv_nucleotide.stk", "stockholm")
rdp_alignment = AlignIO.read("HCV_6_aligned_wrefs.fasta", "fasta")

# Reference sequence IDs
stk_ref_id = "6_DQ278891"
rdp_ref_id = "Ref.6k.CN.x.KM45.DQ278891"

# Extract sequences
stk_ref_seq = [rec for rec in stk_alignment if stk_ref_id in rec.id][0].seq
rdp_ref_seq = [rec for rec in rdp_alignment if rdp_ref_id in rec.id][0].seq

# Extract structure (try different possible keys)
ss_cons = None
possible_keys = [
    "secondary_structure",  # Most likely candidate
    "GC:SS_cons_Alt",      # Alternative consensus
    "GC:SS_cons_Circ",     # Circular RNA structure?
    "GC:SS_cons_miR-122",  # miRNA-related structure
]

for key in possible_keys:
    if key in stk_alignment.column_annotations:
        ss_cons = stk_alignment.column_annotations[key]
        print(f"Using secondary structure from key: '{key}'")
        break

if ss_cons is None:
    raise ValueError("No valid secondary structure annotation found in Stockholm file!")

# Map stk to rdp via aligned reference sequence
map_stk_to_rdp = {}
j = 0
for i in range(len(stk_ref_seq)):
    if stk_ref_seq[i] == '-':
        continue
    while j < len(rdp_ref_seq) and rdp_ref_seq[j] == '-':
        j += 1
    if j < len(rdp_ref_seq):
        map_stk_to_rdp[i] = j
        j += 1

# Prepare data for CSV
csv_data = []
for i, symbol in enumerate(ss_cons):
    if symbol not in ".-":  # focus on structure characters
        rdp_col = map_stk_to_rdp.get(i)
        if rdp_col is not None:
            csv_data.append({
                "STK_Position": i + 1,
                "RDP_Position": rdp_col + 1,
                "Structure_Symbol": symbol
            })

# Write to CSV
csv_filename = "hcv_6_secondary_structure_mapping.csv"
with open(csv_filename, mode='w', newline='') as csvfile:
    fieldnames = ["STK_Position", "RDP_Position", "Structure_Symbol"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    writer.writerows(csv_data)

print(f"Results saved to {csv_filename}")
