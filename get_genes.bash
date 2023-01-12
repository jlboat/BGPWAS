#!/usr/bin/env bash

BED=genes.bed
VCF=Example.vcf.gz

if [ ! -e variant_data/ ]
then
    mkdir variant_data
fi

# For a list of genes ("genes.txt"), prepare the variant data
while read p
do
    # Get gene feature from BED
    grep ${p} ${BED} > ./variant_data/${p}.bed
    # Get coordinate information - CHROM START STOP
    cut -f 1-3 ./variant_data/${p}.bed > ./variant_data/${p}.txt
    # Pull variants within coordinates
    bcftools view -R ./variant_data/${p}.txt \
        -o ./variant_data/${p}.variants.vcf.gz -Oz ${VCF}
    # Generate numeric genotype data
    vcftools --gzvcf ./variant_data/${p}.variants.vcf.gz \
        --012 --out ./variant_data/${p}.variants
done < genes.txt

## Example genes.txt - must match name in BED file
# Sobic.009G229800
# Sobic.006G067700
# Sobic.007G163800
