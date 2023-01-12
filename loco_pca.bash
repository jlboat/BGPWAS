#!/usr/bin/env bash

PLINK_PATH=~/Plink/plink
BASE_NAME=SAP.400

# Convert data from VCF to Plink file format
${PLINK_PATH} \
    --allow-no-sex \
    --make-bed \
    --out ${BASE_NAME} \
    --vcf ${BASE_NAME}.vcf.gz \
    --no-pheno 

# Make directory for PCs
if [ ! -e pcs ]
then
    mkdir pcs
fi

# For each excluded chromosome, perform PCA
for i in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10
do
    ${PLINK_PATH} --bfile ${BASE_NAME} --pca --not-chr ${i} --out exclude_${i}.PCs
    mv exclude_${i}.PCs.eigenvec ./pcs/
    mv exclude_${i}.PCs.eigenval ./pcs/
    rm exclude_${i}.PCs.log exclude_${i}.PCs.nosex
done


