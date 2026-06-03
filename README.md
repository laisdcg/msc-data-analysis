# msc-data-analysis
🧬 R scripts developed during the master’s program for transcriptomic analysis (RNA-Seq), systems biology (WGCNA), and the reconstruction of gene regulatory networks (GRNs) to identify hubs and master regulators.

# Transcriptomic Data Analysis and Gene Regulatory Networks (MSc Thesis)

This repository contains the **R** scripts developed and used during my Master's research program at Graduate Program in Bioinformatics at the Federal University of Rio Grande do Norte.

The primary goal of this workflow is the integrative analysis of public transcriptomic datasets and biological systems modeling to identify key therapeutic targets and biological regulators.

---

## 🔬 Pipeline Overview

The scripts are organized to cover the following methodological steps:

1. **Data Mining & Processing:** Downloading and preparing public datasets from the Gene Expression Omnibus (GEO), specifically: `GSE117769`, `GSE205748`, `GSE221786`, and `GSE186063`.
2. **Differential Gene Expression (DEG):** Rigorous statistical analyses using `DESeq2`, `limma`, and `edgeR` packages, paired with data visualizations using `ggplot2`, `pheatmap`, and `ComplexHeatmap`.
3. **Weighted Gene Co-expression Network Analysis (WGCNA):** Clustering genes into functionally correlated modules.
4. **Functional Enrichment Analysis:** Gene Ontology (GO) and biological pathway analysis (KEGG/Reactome) using `clusterProfiler`, `topGO`, and `ReactomePA`.
5. **Gene Regulatory Network (GRN) Reconstruction:** Integrating Gene-Transcription Factor (TF) interactions based on transcription factor binding motif mapping.
6. **Network Topology Analysis:** Utilizing `igraph` and `ggraph` for centrality inferences to discover:
   * **Hubs:** Highly connected genes within the network.
   * **Master Regulators (MRs):** Transcription Factors (TFs) with high out-degree centrality, indicating potential control over major biological modules.

---

## 🛠️ Tech Stack & Main Packages

The project was entirely developed in **R** (version 4.x or higher). The key packages used include:

* **CRAN:** `tidyverse` (`dplyr`, `tidyr`, `ggplot2`), `igraph`, `ggraph`, `umap`, `factoextra`, `pheatmap`.
* **Bioconductor:** `GEOquery`, `DESeq2`, `limma`, `edgeR`, `WGCNA`, `clusterProfiler`, `org.Hs.eg.db`, `AnnotationDbi`, `ComplexHeatmap`, `ReactomePA`, `GENIE3`.

---

## 📂 Code Structure

* Scripts named `GSEXXXXXX_code.R` contain the individual pipeline for processing, normalization, differential expression, WGCNA, and enrichment for each specific dataset.
* The script `GRN_Hubs_MRs.R` handles the final biological network modeling, network topology metrics calculation (in-degree/out-degree), and outputs the visual representations of the main Hubs and Master Regulators.

---

## 🎓 Academic Reference

* **Thesis Title:** Beyond the Shared Inflammatory Axis: Differentiating Molecular Signatures in Psoriatic Arthritis and Ankylosing Spondylitis through Integrated Omics.
* **Author:** Laís de Carvalho Gonçalves
* **Advisor:** Prof. Dr. João Paulo Matos Santos Lima
* **Institution:** Universidade Federal do Rio Grande do Norte

If you use any part of this code or methodology in your research, please cite our work:
- Preprint: Gonçalves, L. C., Rodrigues-Neto, J. F., Gupta, S., de Souza, G. A., & Lima, J. P. M. S. (2026). Beyond the Shared Inflammatory Axis: Differentiating Molecular Signatures in Psoriatic Arthritis and Ankylosing Spondylitis through Integrated Omics. bioRxiv. (https://www.biorxiv.org/content/10.1101/2025.08.20.671331v1)
- Peer-Reviewed Article: Gonçalves, L. C., Rodrigues-Neto, J. F., Gupta, S., de Souza, G. A., & Lima, J. P. M. S. (2026). Beyond the Shared Inflammatory Axis: Differentiating Molecular Signatures in Psoriatic Arthritis and Ankylosing Spondylitis through Integrated Omics. Genes & Diseases (In Press).
---
*For questions, feedback, or collaborations, feel free to open an Issue or contact me via [lais.goncalves.034@ufrn.edu.br].*
