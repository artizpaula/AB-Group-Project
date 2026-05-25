# PROJECT 4 — Conservation Analysis from MSA

## Functional vs Non-Functional Regions + Coevolution

This repository contains the work for **Project 4** of the course **Algorithms in Biology (2025–2026)** in the Bachelor's Degree in Bioinformatics.

The project analyzes conservation and coevolution in a protein family using a **Multiple Sequence Alignment (MSA)** from **Pfam (PF00069 — Protein Kinase family)**.

---

## Authors

- Alicia Mañas  
- Paula Artiz  
- Lídia Sánchez  

---

# Project Overview

The objectives of this project are:

- Compute per-column conservation using **Shannon entropy**
- Identify conserved and variable regions
- Relate conserved positions to **UniProt functional annotations**
- Detect coevolving residues using **Mutual Information (MI)** and **APC correction**
- Visualize conservation and coevolution patterns

Reference protein:

- **AKT1_HUMAN (UniProt: P31749)**

---

# Repository Contents

```bash
├── PF00069.alignment.seed
├── PF00069_seed.fasta
├── project4.py
├── entropy_plot.png
├── entropy_annotated.png
├── mi_heatmap.png
└── README.md
```

---

# Methods

## Conservation Analysis

Conservation is measured using Shannon entropy:
- Low entropy → conserved positions
- High entropy → variable positions

---

## Coevolution Analysis

Residue coevolution is computed using Mutual Information.
APC correction is applied to reduce background noise.

---

# Main Results

- 37 aligned protein sequences
- 419 alignment columns
- 52 highly conserved positions detected
- 243 highly variable positions detected
- Conserved residues overlap with:
  - active sites
  - binding sites
  - functional motifs

Top coevolving pairs were identified using APC-corrected MI scores.

---

# Tools Used

- Python
- NumPy
- Matplotlib
- Biopython
- Requests

---

# Outputs

- `entropy_plot.png` → conservation profile
- `entropy_annotated.png` → conservation + functional annotations
- `mi_heatmap.png` → coevolution heatmap

---

# Conclusion

This project demonstrates how evolutionary information from protein MSAs can identify:

- conserved functional residues
- variable regions
- coevolving amino acid networks

useful for studying protein structure and function.
