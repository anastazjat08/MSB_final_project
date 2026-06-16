# Investigating the relationship between mutational spectra and antimicrobial‑resistance burden in MRSA
This repository contains a reproducible Snakemake workflow and analysis code for investigating whether mutational spectra are associated with antimicrobial‑resistance (AMR) burden and multidrug‑resistance (MDR) levels in *Staphylococcus aureus* (MRSA).

The project integrates AMRFinder outputs, recombination‑filtered phylogenies, MutTui mutational spectra, PCA and linear modelling.

## Repository structure

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
