# Investigating the relationship between mutational spectra and antimicrobial‑resistance burden in MRSA
This repository contains a reproducible Snakemake workflow and analysis code for investigating whether mutational spectra are associated with antimicrobial‑resistance (AMR) burden and multidrug‑resistance (MDR) levels in *Staphylococcus aureus* (MRSA).

The project integrates AMRFinder outputs, recombination‑filtered phylogenies, MutTui mutational spectra, PCA and linear modelling.

## Repository structure
```
MSB_final_project/
├── .gitignore
├── 2_scripts/
│   ├── amr_analysis.py
│   ├── build_super_tables.py
│   ├── generate_labels.py
│   ├── plot_sbs_stacked.py
│   ├── spectra_PC1_linear_models.py
│   ├── spectra_pca_96_channels.py
│   └── spectra_sig_linear_models.py
├── config/
│   └── config.yaml
├── environment.yaml
├── project_plan.pdf
├── README.md
├── report.pdf
└── Snakefile
```
## Scripts overview
This project uses several Python scripts that implement individual steps of the Snakemake workflow:

- amr_analysis.py - merges AMRFinder outputs with metadata, computes AMR burden and MDR levels and generates summarising plots.

- build_super_tables.py - constructs combined “supertables” containing 6, 96‑channel mutational spectra, metadata, and study labels used for PCA and modelling.

- generate_labels.py - generates sample‑level labels (study, lineage, AMR categories) required for MutTui and PCA visualisation.

- plot_sbs_stacked.py - produces stacked barplots of 96‑channel SBS mutational spectra aggregated by study and sequence type.

- spectra_pca_96_channels.py - performs PCA on 96‑channel spectra, filters low‑mutation samples, and outputs PCA coordinates and diagnostic plots.

- spectra_PC1_linear_models.py - fits linear models testing the association between PC1 and AMR burden, with and without study adjustment.

- spectra_sig_linear_models.py - fits linear models linking individual SBS signatures (e.g., C>A, T>C) to AMR burden and study effects.


## Environment
The workflow uses:
- python 3.10
- snakemake
- AMRFinderPlus
- Parsnp
- Harvesttools
- Gubbins
- IQ-TREE
- MutTui
All dependencies are defined in environment.yaml

### Create the environment
```
conda env create -f environment.yml
conda activate mrsa_env
```

## Running the pipeline
To execute the full workflow:
```
snakemake --core 8
```
This will:
1. Merge AMRFinder results with metadata
2. Compute MDR levels,  AMR burden and perform some analysis
3. Run Parsnp, Gubbins and IQ-TREE
4. Generate MutTui spectra (for the whole tree and specified labels)
5. Perform PCA
6. Fit linear models linking PC1 and SBS signatures to AMR burden
7. Save plots and tabels to 1_data/2_processed and 1_data/3_results

## Key results
### AMR burden & MDR
- AMR burden distribution and MDR levels computed from AMRFinder metadata

- Significant differences in AMR class frequencies between studies (z‑tests)

### Mutational spectra
- 6, 96‑channel SBS spectra generated using MutTui

### Linear models
- Model 1: AMR burden alone -> no association

- Model 2: AMR burden + study -> significant positive association

- Model 3: AMR burden × study -> no consistent study‑independent effect
