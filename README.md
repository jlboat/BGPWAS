# BGPWAS
A sparse Bayesian Genome-Phenome Wide Association Study model

## Required items
 * VCF file with variants of interest (Example.vcf.gz)
 * BED file with gene-level information (genes.bed)
 * Phenomic data matrix (unimputed\_phenotypes.tsv)
 * Kinship matrix (Kinship.csv.zip)

 * Plink - see <https://www.cog-genomics.org/plink/>
 * R - see  <https://www.r-project.org/>
 * PHENIX - R package, see <https://github.com/jlboat/PHENIX>
 * bayesplot - R package
 * rstanarm - R package

For the SAP, data are available as follows:  

Raw read data: <https://www.ebi.ac.uk/ena/browser/view/PRJEB50066>
Variants: <https://www.ebi.ac.uk/ena/browser/view/ERZ12588732>

## Workflow
Using the VCF file, perform LOCO PCA across all chromosomes. *But don't run with example data here except to test*.

```bash
bash loco_pca.bash # will overwrite valid PCs provided here
```

Using the BED file, extract gene coordinates and then corresponding variants.

```bash
bash get_genes.bash
```

Using the phenomic data matrix and kinship matrix, impute the missing phenotypes.

```bash
unzip Kinship.csv.zip
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

Generate figures for each gene/variant/model with significant associations.
```bash
Rscript plot_gpwas.R Sobic.006G067700.Chr06_42806735 hs
```

Citation (*Pending*)
```
J. Lucas Boatwright, Sirjan Sapkota, and Stephen Kresovich. Functional Genomic Effects of Indels using Bayesian Genome-Phenome Wide Association Studies in Sorghum. *Pending*
```
