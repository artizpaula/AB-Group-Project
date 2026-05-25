"""
PROJECT 4 - Conservation Analysis from MSA
Functional vs Non-Functional Regions + Coevolution
Lídia Sánchez, Alicia Mañas and Paula Artiz
"""

# 1. Obtaining and preparation of the MSA + data

# Library imports

from Bio import AlignIO
import numpy as np
import requests
import matplotlib.pyplot as plt

# Load Stockholm alignment file

alignment = AlignIO.read("PF00069.alignment.seed", "stockholm") # Protein Kinase family (Pfam PF00069)

print(f"Number of sequences : {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()} columns")

# Export alignment in FASTA format for manual inspection

AlignIO.write(alignment, "PF00069_seed.fasta", "fasta")

# Represent the MSA as a character matrix (rows = sequences, columns = positions)

msa = np.array([list(str(record.seq).upper()) for record in alignment])

# Replace any non-standard character with a gap symbol

valid_aa = set("ACDEFGHIKLMNPQRSTVWY-")

for i in range(msa.shape[0]):
    for j in range(msa.shape[1]):
        if msa[i, j] not in valid_aa:
            msa[i, j] = "-"

# Encode amino acids as integers for numerical computation (gaps encoded as -1)

aa2i = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}

m = np.full(msa.shape, -1, dtype=np.int16)

for i in range(msa.shape[0]):
    for j in range(msa.shape[1]):
        aa = msa[i, j]
        if aa in aa2i:
            m[i, j] = aa2i[aa]

# Basic quality check: compute gap fraction per column

gap_fraction = np.mean(m == -1, axis=0)
print(f"Gap fraction — mean: {gap_fraction.mean():.3f}  max: {gap_fraction.max():.3f}")

# Export prepared matrices for downstream analysis

np.save("msa_chars.npy", msa) # Character matrix
np.save("msa_int.npy", m)     # Integer-encoded matrix

# 2. Shannon entropy and functional analysis

# Load the integer-encoded MSA matrix

m = np.load("msa_int.npy")
n_seqs, n_cols = m.shape
n_aa = 20 # Number of standard amino acids

entropy = np.zeros(n_cols)

for j in range(n_cols):
    col = m[:, j] # Extract all residues at position j
    col = col[col != -1]# Exclude gap-encoded positions

    if len(col) == 0: # Skip columns that are entirely gaps
        entropy[j] = np.nan
        continue

    counts = np.bincount(col, minlength=n_aa) # Observed amino acid counts per column
    freqs = counts / counts.sum() # Relative frequencies
    freqs = freqs[freqs > 0]
    entropy[j] = -np.sum(freqs * np.log2(freqs)) # Shannon entropy formula

# Plot per-column conservation
plt.figure(figsize=(14, 4))
plt.plot(entropy, linewidth=0.8, color="steelblue")
plt.xlabel("Position in the alignment")
plt.ylabel("Shannon entropy (bits)")
plt.title("Conservation per column — PF00069 (Protein Kinase)")
plt.axhline(y=0.5, color="green", linestyle="--", linewidth=1, label="Very conserved (<0.5)")
plt.axhline(y=2.0, color="red", linestyle="--", linewidth=1, label="Very variable (>2.0)")
plt.legend()
plt.tight_layout()
plt.savefig("entropy_plot.png", dpi=150)
plt.show()

# Identify and report conserved and variable positions
conserved_positions = np.where(entropy < 0.5)[0]
variable_positions = np.where(entropy > 2.0)[0]

print(f"Positions very conserved ({len(conserved_positions)}): {conserved_positions}")
print(f"Positions very variable  ({len(variable_positions)}): {variable_positions}")
np.save("entropy.npy", entropy)

# 3. Functional Annotation

# Entropy data
entropy = np.load("entropy.npy")
conserved_set = set(np.where(entropy < 0.5)[0]) # Set of conserved column indices

# Locate AKT1_HUMAN in the alignment and retrieve its sequence
for i, record in enumerate(alignment):
    if "AKT1_HUMAN" in record.id:
        ref_seq = str(record.seq)
        break

# Build a mapping from MSA column indices to UniProt residue positions
msa_to_uniprot = {}
uniprot_pos = 150 # UniProt numbering begins at residue 150 for this fragment
for msa_pos, aa in enumerate(ref_seq):
    if aa != "-": # Only map non-gap positions
        msa_to_uniprot[msa_pos] = uniprot_pos
        uniprot_pos += 1

# Invert the mapping to enable lookup by UniProt position
uniprot_to_msa = {v: k for k, v in msa_to_uniprot.items()} # Reverse look-up

# Retrieve functional annotations for AKT1 (UniProt accession P31749) -> no manual download needed
response = requests.get("https://rest.uniprot.org/uniprotkb/P31749.json")
data = response.json()

features = []
for feature in data.get("features", []):
    feat_type = feature.get("type", "")
    if feat_type in ["Active site", "Binding site", "Region", "Motif"]: # Retain only functionally relevant annotations
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
        desc = feature.get("description", "")
        features.append({"type": feat_type, "description": desc, "start": start, "end": end})

# Cross-reference functional annotations with conserved MSA positions
print("Functional annotations in CONSERVED regions")
found = False
for feat in features:
    for up_pos in range(feat["start"], feat["end"] + 1):
        msa_pos = uniprot_to_msa.get(up_pos) # Retrieve the corresponding MSA column
        if msa_pos is not None and msa_pos in conserved_set: # Check whether the position is conserved
            print(f"  MSA pos {msa_pos} | UniProt pos {up_pos} | {feat['type']} | {feat['description']}")
            found = True

if not found:
    print("  No functional annotations found in conserved regions")

# Plot conservation profile with functional annotations overlaid

plt.figure(figsize=(14, 4))
plt.plot(entropy, linewidth=0.8, color="steelblue", label="Shannon entropy")
plt.axhline(y=0.5, color="green", linestyle="--", linewidth=1, label="Very conserved (<0.5)")
plt.axhline(y=2.0, color="red", linestyle="--", linewidth=1, label="Very variable (>2.0)")

# Mark each functional annotation as a vertical line at its MSA position
colors = {"Active site": "red", "Binding site": "orange", "Region": "purple", "Motif": "brown"}

for feat in features:
    for up_pos in range(feat["start"], feat["end"] + 1):
        msa_pos = uniprot_to_msa.get(up_pos)
        if msa_pos is not None:
            color = colors.get(feat["type"], "gray")
            plt.axvline(x=msa_pos, color=color, alpha=0.4, linewidth=1.2)

# Construct a custom legend for the annotation types
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color="red",    alpha=0.6, label="Active site"),
    Line2D([0], [0], color="orange", alpha=0.6, label="Binding site"),
    Line2D([0], [0], color="green",  linestyle="--", label="Very conserved (<0.5)"),
    Line2D([0], [0], color="red",    linestyle="--", label="Very variable (>2.0)"),
]
plt.legend(handles=legend_elements, loc="upper right")
plt.xlabel("Position in the alignment")
plt.ylabel("Shannon entropy (bits)")
plt.title("Conservation per column with functional annotations — PF00069 (AKT1_HUMAN)")
plt.tight_layout()
plt.savefig("entropy_annotated.png", dpi=150)
plt.show()

# 4. Coevolution (Mutual information MI)

def column_mi(col_i, col_j, n_aa=20):
    """Compute mutual information between two integer-encoded alignment columns (gaps excluded pairwise)."""
    mask = (col_i != -1) & (col_j != -1)
    ci, cj = col_i[mask], col_j[mask]
    n = len(ci)
    if n < 2:
        return 0.0
    # Compute joint frequency matrix
    joint = np.zeros((n_aa, n_aa), dtype=np.float64)
    np.add.at(joint, (ci, cj), 1)
    joint /= n
    # Compute marginal frequencies
    pi = joint.sum(axis=1) # marginal for column i
    pj = joint.sum(axis=0) # marginal for column J
    # Compute mutual information -> MI = Σp(a,b)·log2[ p(a,b)/(p(a)·p(b))]
    mi = 0.0
    for a in range(n_aa):
        for b in range(n_aa):
            if joint[a, b] > 0 and pi[a] > 0 and pj[b] > 0:
                mi += joint[a, b] * np.log2(joint[a, b] / (pi[a] * pj[b]))
    return mi

# Restrict analysis to columns with a gap fraction below 50% to reduce noise
gap_frac = np.mean(m == -1, axis=0)
valid_cols = np.where(gap_frac < 0.5)[0]
n_valid    = len(valid_cols)
print(f"\nValid columns:{n_cols}/{alignment.get_alignment_length()}")
print(f"Gap fraction threshold : 50%")

# Compute the full pairwise MI matrix (symmetric; only upper triangle calculated)
MI_raw = np.zeros((n_valid, n_valid), dtype=np.float64)
for idx_i, i in enumerate(valid_cols):
    for idx_j in range(idx_i + 1, n_valid):
        j = valid_cols[idx_j]
        val = column_mi(m[:, i], m[:, j])
        MI_raw[idx_i, idx_j] = val
        MI_raw[idx_j, idx_i] = val

# Apply Average Product Correction (APC) to remove background MI
# Formula: MIp(i,j) = MI(i,j) - <MI(i,*)> * <MI(*,j)> / <MI>
row_mean   = MI_raw.mean(axis=1, keepdims=True)  # Row means <MI(i,*)>
col_mean   = MI_raw.mean(axis=0, keepdims=True)  # Column means <MI(*,j)>
grand_mean = MI_raw.mean()                       # Global mean <MI>
MIp = MI_raw - (row_mean * col_mean) / (grand_mean + 1e-12)
np.fill_diagonal(MIp, 0) # a position doesn't coevolve with itself


# Identify and report the top-ranked coevolving position pairs
n_top = 20
triu  = np.triu_indices(n_valid, k=1)
scores = MIp[triu]
order  = np.argsort(scores)[::-1][:n_top]

print(f"\nTop {n_top} coevolving pairs (APC-corrected MI):")
print(f"{'MSA_i':>6}  {'MSA_j':>6}  {'MIp':>8}  {'UniProt_i':>10}  {'UniProt_j':>10}")
print("-" * 55)
top_pairs = []
for k in order:
    idx_i, idx_j = triu[0][k], triu[1][k]
    col_i = int(valid_cols[idx_i])
    col_j = int(valid_cols[idx_j])
    mip   = scores[k]
    up_i  = msa_to_uniprot.get(col_i, "—")
    up_j  = msa_to_uniprot.get(col_j, "—")
    top_pairs.append((col_i, col_j, mip))
    print(f"{col_i:6d}  {col_j:6d}  {mip:8.4f}  {str(up_i):>10}  {str(up_j):>10}")

# Plot pairwise coevolution as a heatmap
plt.figure(figsize=(8, 7))
plt.imshow(MIp, cmap="magma", aspect="auto",
           interpolation="none", vmin=0, vmax=np.percentile(MIp, 99))
plt.colorbar(label="MIp (bits)")

tick_step   = max(1, n_valid // 8)
tick_idx    = list(range(0, n_valid, tick_step))
tick_labels = [str(valid_cols[i]) for i in tick_idx]
plt.xticks(tick_idx, tick_labels, fontsize=7)
plt.yticks(tick_idx, tick_labels, fontsize=7)
plt.xlabel("Alignment column")
plt.ylabel("Alignment column")
plt.title("Pairwise coevolution heatmap (MIp) — PF00069 (Protein Kinase)")
plt.tight_layout()
plt.savefig("mi_heatmap.png", dpi=150) # Heatmap: bright spots = pairs that coevolve strongly
plt.show()


"""
Output:
Number of sequences : 37
Alignment length    : 419 columns
Gap fraction — mean: 0.361  max: 0.973
Positions very conserved (52): [  7   9  15  16  31  55  83  85  86 123 146 160 172 174 177 184 185 186
 187 188 189 190 191 192 193 209 210 211 228 229 230 231 232 233 238 253
 254 262 269 270 271 280 285 361 362 363 364 365 366 367 368 401]
Positions very variable  (243): [  1   2   3   4   5   8  10  11  13  17  18  19  20  21  22  23  24  25
  26  27  28  30  32  33  34  35  36  37  38  39  40  41  42  43  44  45
  48  49  50  51  52  53  54  56  57  58  59  60  61  62  63  66  68  69
  72  73  74  75  76  77  78  79  80  92  93  94  95  96  97 100 101 102
 104 105 106 107 109 115 116 117 118 119 120 125 126 127 128 129 130 132
 133 135 143 144 145 147 148 149 150 151 152 153 154 157 162 163 164 169
 176 180 181 182 183 194 195 196 203 204 208 212 214 215 216 217 218 219
 220 221 222 239 240 241 242 243 244 245 248 249 251 255 266 267 268 272
 273 275 276 277 278 279 281 284 286 287 288 289 290 292 293 295 296 299
 300 301 302 303 304 305 306 307 308 309 310 311 312 313 315 316 317 318
 319 320 321 322 334 335 336 337 338 339 340 341 342 343 344 345 346 347
 348 349 350 351 352 353 354 355 356 357 358 371 372 373 374 376 377 378
 379 380 381 382 384 385 386 387 388 390 391 392 395 396 399 400 402 404
 407 409 410 411 412 413 414 416 418]
Functional annotations in CONSERVED regions
  MSA pos 172 | UniProt pos 274 | Active site | Proton acceptor
  MSA pos 7 | UniProt pos 157 | Binding site | 
  MSA pos 9 | UniProt pos 159 | Binding site | 
  MSA pos 31 | UniProt pos 179 | Binding site | 
Number of sequences : 37
Alignment length    : 419 columns
Gap fraction — mean: 0.361  max: 0.973
Positions very conserved (52): [  7   9  15  16  31  55  83  85  86 123 146 160 172 174 177 184 185 186
 187 188 189 190 191 192 193 209 210 211 228 229 230 231 232 233 238 253
 254 262 269 270 271 280 285 361 362 363 364 365 366 367 368 401]
Positions very variable  (243): [  1   2   3   4   5   8  10  11  13  17  18  19  20  21  22  23  24  25
  26  27  28  30  32  33  34  35  36  37  38  39  40  41  42  43  44  45
  48  49  50  51  52  53  54  56  57  58  59  60  61  62  63  66  68  69
  72  73  74  75  76  77  78  79  80  92  93  94  95  96  97 100 101 102
 104 105 106 107 109 115 116 117 118 119 120 125 126 127 128 129 130 132
 133 135 143 144 145 147 148 149 150 151 152 153 154 157 162 163 164 169
 176 180 181 182 183 194 195 196 203 204 208 212 214 215 216 217 218 219
 220 221 222 239 240 241 242 243 244 245 248 249 251 255 266 267 268 272
 273 275 276 277 278 279 281 284 286 287 288 289 290 292 293 295 296 299
 300 301 302 303 304 305 306 307 308 309 310 311 312 313 315 316 317 318
 319 320 321 322 334 335 336 337 338 339 340 341 342 343 344 345 346 347
 348 349 350 351 352 353 354 355 356 357 358 371 372 373 374 376 377 378
 379 380 381 382 384 385 386 387 388 390 391 392 395 396 399 400 402 404
 407 409 410 411 412 413 414 416 418]
Functional annotations in CONSERVED regions
  MSA pos 172 | UniProt pos 274 | Active site | Proton acceptor
  MSA pos 7 | UniProt pos 157 | Binding site | 
  MSA pos 9 | UniProt pos 159 | Binding site | 
  MSA pos 31 | UniProt pos 179 | Binding site | 

Valid columns:419/419
Gap fraction threshold : 50%

Top 20 coevolving pairs (APC-corrected MI):
 MSA_i   MSA_j       MIp   UniProt_i   UniProt_j
-------------------------------------------------------
   300     322    0.9894         351           —
   300     316    0.9143         351           —
   129     183    0.8749         243         284
   300     319    0.8221         351           —
   316     390    0.6878           —         380
   297     318    0.6762         348           —
   307     313    0.6655         358         364
   306     322    0.6565         357           —
   118     322    0.6559         239           —
   246     318    0.6554         311           —
   300     317    0.6398         351           —
   105     306    0.6346         230         357
   118     317    0.6330         239           —
   104     290    0.6322         229         341
   302     322    0.6322         353           —
   284     318    0.6065         335           —
   318     343    0.5957           —         369
   109     251    0.5888         234         316
   298     317    0.5852         349           —
   278     309    0.5841         329         360
"""