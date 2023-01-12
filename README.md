# BGPWAS
A sparse Bayesian Genome-Phenome Wide Association Study model

## Required items
 * VCF file with variants of interest
 * BED file with gene-level information
 * Phenomic data matrix
 * Kinship matrix

 * bayesplot
 * rstanarm
 * PHENIX - see <https://github.com/jlboat/PHENIX>

## Workflow
Using the VCF file, perform LOCO PCA across all chromosomes. *But don't run with example data here except to test*

```bash
bash loco_pca.bash # will overwrite valid PCs provided here
```

Using the BED file, extract gene coordinates and then corresponding variants.

```bash
bash get_genes.bash
```

Using the Phenomic data matrix and kinship matrix, impute the missing phenotypes.

```bash
Rscript impute_phenotypes.R
```

Perform BGPWAS with the variants, PCs, and imputed phenotypes.
```bash
# Horseshoe
Rscript bgpwas.R Sobic.009G229800 hs

# Lasso
Rscript bgpwas.R Sobic.009G229800 lasso

# Ridge
Rscript bgpwas.R Sobic.009G229800 ridge
```

Generate figures for each gene/variant/model with significant associations
```bash
Rscript plot_gpwas.R Sobic.006G067700.Chr06_42806735 hs
```
