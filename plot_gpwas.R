suppressMessages(library(shinystan, quietly=T, verbose=F))
suppressMessages(library(bayesplot, quietly=T, verbose=F))
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

if (is.na(args[1])){
    print("ERROR: No file included.")
    print("Rscript plot_gpwas.R gene_and_variant_id model_type")
    quit("no")
}

# Optionally, change posterior x limits
# xmin = 5
# xmax = 5

i = args[1]
model_type = args[2]
fit <- readRDS(paste0(i, ".variants.", model_type, ".rds"))
# Optionally, launch shiny app of results
# launch_shinystan(fit)
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
               regex_pars=c("_"),
               prob_outer = 0.95) # + xlim(xmin, xmax)
dev.off()

significant_parameters = significant_parameters[grep("_", significant_parameters)]

png(paste0(i, ".significant_beta.", model_type, ".png"), width = 5, height = 14, units = 'in', res = 300)
mcmc_areas(posterior,
           prob_outer = 0.95,
           pars = significant_parameters
           ) + theme(axis.text=element_text(size=20)) # + xlim(xmin, xmax)
dev.off()



