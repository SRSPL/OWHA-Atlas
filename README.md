
# Organ-Scale Wound Healing Atlas (OWHA)


## Overview

The Organ-Scale Wound Healing Atlas (OWHA) is a 4D multimodal reference atlas that reconstructs the complete spatial and temporal choreography of mammalian wound healing at single-cell resolution.

OWHA integrates across the full wound healing timeline.:

- Single-nucleus RNA sequencing (snRNA-seq)
- Single-cell RNA sequencing (scRNA-seq)
- Cellular Indexing of Transcriptomes and Epitopes by sequencing (CITE-seq)
- Vissium HD (High-definition) spatial transcriptomics


Preprint:
Chin Cheong et al., Organ-Scale Wound Healing Atlas, bioRxiv (2026)
https://www.biorxiv.org/content/10.64898/2026.01.15.699736v1

## Rationale

Deep skin wounds require tightly coordinated communication across epithelial, vascular, neural, immune, and stromal systems. Existing wound atlases capture partial views of this process, often limited to specific compartments, timepoints, or dissociation-compatible cell types.

OWHA overcomes these limitations by profiling:

 - Over 725,000 murine single-cell and spatial transcriptomes

 - All major skin microanatomical niches

 - Early to late healing phases

- More than 100 precisely annotated cell states

- Fragile and historically underrepresented populations, including adipocytes, Schwann cells, and transient epithelial intermediates

This atlas enables organ-scale systems reconstruction of wound repair.

## Key Biological Findings
### Central Orchestrator Populations

OWHA reveals that wound repair proceeds through sharp transcriptional and cellular inflection points driven by transient regulatory hubs termed Central Orchestrator populations. These populations synchronize repair across tissue systems through coordinated transcriptional activation and cross-tissue signaling.

### Basal IV Keratinocytes

OWHA identifies a previously uncharacterized epithelial subpopulation:

- Sox6+ Tspear+ Il20ra+ keratinocytes ("Basal IV")

After injury, Basal IV cells deviate from canonical differentiation programs and adopt a neurovasculogenic signaling state. They form a transient spatially privileged niche aligned with proliferative endothelial cells, pericytes, and Repair Schwann cells.

This niche synchronizes:

- Re-epithelialization

- Angiogenesis

- Neurite guidance

Mechanistically, this coordination is mediated through a conserved Sema3C–Nrp1/Nrp2 signaling axis.

### Human Conservation and Disease Relevance

Cross-species integration demonstrates that:

- The Basal IV / SEMA3C axis is conserved in human skin

- It is missed in conventional human scRNA-seq atlases due to dissociation-induced artifacts

- It is selectively disrupted in diabetic wounds

Topical Sema3C restores peri-wound angiogenic sprouting and accelerates re-epithelialization in diabetic ulcers in vivo.

## Repository Structure

The repository contains processed datasets, annotations, and reproducible analysis workflows corresponding to the bioRxiv release.

### 1. Figures

```
figures/
```

Contains Scripts Relatedted to Producing:

- Publication figures

- UMAP visualizations

- Spatial overlays

### 2. Analysis Pipelines

```
scripts/
```

Reproducible workflows for:

- Preprocessing and quality control

- Dataset integration (Seurat-based)

- Differential expression analysis

- Gene Ontology enrichment

- Module scoring

- Cross-species integration

Scripts correspond to the bioRxiv release.

## Example: Loading the Data in R

```
library(Seurat)
owha <- readRDS("OWHA_Integrated.rds")
DimPlot(owha, group.by = "celltype")
Citation
```

## If you use OWHA, please cite:

Chin Cheong et al. Organ-Scale Wound Healing Atlas. bioRxiv (2026).
