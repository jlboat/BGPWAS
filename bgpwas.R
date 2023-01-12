library(rstan)
library(rstanarm)
suppressMessages(library(shinystan, quietly=T, verbose=F))
suppressMessages(library(bayesplot, quietly=T, verbose=F))

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

phenomic_data = c("imputed_phenotypes.csv")
pc_count = 3

args <- commandArgs(trailingOnly=TRUE)

if (is.na(args[1])){
    print("ERROR: No file included.")
    print("Rscript bgpwas.R gene_ID model_type")
    quit("no")
}

clean_correlated <- function(df, cor_coef=0.95){
    if (dim(df)[2] <= 1){
        return(df)
    }
    df <- df[!duplicated(as.list(df))]
    cor_matrix <- cor(df)
    cor_matrix_rm <- cor_matrix
    cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
    diag(cor_matrix_rm) <- 0
    clean_df <- df[, !apply(cor_matrix_rm,
                           2,
                           function(x) any(x > cor_coef)),drop=F]
    return(clean_df)
}

plot_posteriors <- function(fit, i, model_type){
    posterior <- as.array(fit)
    interval <- mcmc_intervals(posterior,
                               prob_outer = 0.95)
    
    significant_parameters <- as.character(
                interval$data$parameter)[((interval$data$ll < 0) & 
                                         (interval$data$hh < 0)) | 
                ((interval$data$ll > 0) & (interval$data$hh > 0))]
    print(significant_parameters)
    
    write.table(as.data.frame(significant_parameters), file=paste0(i, ".", model_type, ".sig_traits.txt"), 
                quote=F, row.names=F, col.names = F)
    
    png(paste0(i, ".beta.", model_type, ".png"), width = 4, height = 14, units = 'in', res = 300)
    mcmc_intervals(posterior,
                   prob_outer = 0.95)
    dev.off()
    
    png(paste0(i, ".significant_beta.", model_type, ".png"), width = 4, height = 14, units = 'in', res = 300)
    mcmc_areas(posterior,
                   prob_outer = 0.95,
                   pars = significant_parameters)
    dev.off()
}

geno_pheno_columns = list()
# Read in pheno data
i = args[1] # args[1] = Sobic.ID
model_type = args[2] # ridge, lasso, or hs
chrom <- substr(i, 8,9)
geno <- read.table(paste0("./variant_data/", i, ".variants.012"), 
                   header=F, row.names=1)
header <- read.csv(paste0("./variant_data/", i, ".variants.012.indv"), 
                   header=F)
rownames(geno) <- header$V1
positions <- read.table(paste0("./variant_data/", i, ".variants.012.pos"), 
                        header=F)
positions$names <- paste(positions$V1, positions$V2, sep="_")
colnames(geno) <- positions$names
# remove duplicates
geno1 <- clean_correlated(geno)
geno_pheno_columns[[i]] <- colnames(geno1)

# pheno <- read.table("GCA_Hyb_males_yieldComp.txt", header=T, row.names=1)
# pheno <- read.csv("SAP_BLUPs.imputed.csv", header=T, row.names=1)
pheno <- read.csv(phenomic_data, header=T, row.names=1)
pheno <- clean_correlated(pheno)
pcs <- read.table(paste0("./pcs/exclude_Chr", 
                         chrom, 
                         ".PCs.eigenvec"), 
                  row.names=1)
pc_columns <- c()
for (j in 1:pc_count){
    pc_columns <- c(pc_columns, paste0("PC", j))
}
colnames(pcs) <- c("Taxa", pc_columns)
common <- intersect(intersect(rownames(geno), rownames(pheno)), rownames(pcs))
pheno <- pheno[common,]
pcs <- pcs[common,]
pheno_pc <- cbind(pheno[common,], pcs[,2:length(pc_columns)+1])
# pheno_pc <- scale(pheno_pc, center=T, scale=T)
geno_pheno_columns[["pheno"]] <- colnames(pheno_pc)

geno1 <- geno1[common,,drop=F]
if (dim(geno1)[2] == 0){
    print(paste0(i, " has no variants!"))
    break
}
saveRDS(t(geno1), file=paste0("./variant_data/", i, ".geno.rds"))
gpwas_data <- list(K=ncol(geno1),
                   J=ncol(pheno_pc), 
                   N=nrow(geno1), 
                   x=pheno_pc, 
                   y=geno1)
hs_data <- cbind(geno1, pheno_pc)
for (variant in 1:length(colnames(geno1))){
    formula <- eval(parse(text=paste0(paste(colnames(geno1)[variant], collapse=", "), 
                                      " ~ ", paste(colnames(pheno_pc), collapse=" + "))))
    if (model_type == "hs"){
        D <- gpwas_data$J
        n <- gpwas_data$N
        p0 <- ceiling(0.1 * D)
        tau0 <- p0/(D-p0) /sqrt(n)
        prior_coeff <- hs ( df =1 , global_df =1 , global_scale = tau0 )

        fit <- stan_glm(formula, gaussian(), prior=prior_coeff, data=hs_data, core=4)
        saveRDS(fit, file=paste0(i, ".", colnames(geno1)[variant], ".variants.hs.rds"))
    }
    if (model_type == "ridge"){
        prior_coeff <- normal(0, 5)
        fit <- stan_glm(formula, gaussian(), prior=prior_coeff, data=hs_data, core=4)
        saveRDS(fit, file=paste0(i, ".", colnames(geno1)[variant], ".variants.ridge.rds"))
    }
    if (model_type == "lasso"){
        prior_coeff <- lasso(df = 1, scale = 5)
        fit <- stan_glm(formula, gaussian(), prior=prior_coeff, data=hs_data, core=4)
        saveRDS(fit, file=paste0(i, ".", colnames(geno1)[variant], ".variants.lasso.rds"))
    }
    plot_posteriors(fit, paste0(i, ".", colnames(geno1)[variant]), model_type)
}

