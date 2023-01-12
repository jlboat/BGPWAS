library(phenix)

df <- read.table("unimputed_phenotypes.tsv", row.names=1, head=T)
df$Tannins_F <- unname(setNames(0:1, c("Non-tannin", "Tannin"))[df$Tannins_F])
df$AnthraClassification_P <- unname(setNames(0:1, c("Suceptible", "Resistant"))[df$AnthraClassification_P])

K <- read.csv("Kinship.csv", row.names=1, head=T)
colnames(K) <- rownames(K)

common <- intersect(rownames(K), rownames(df))
K <- K[common, common]
df <- df[common,]

df <- df[,colMeans(is.na(df)) < 0.3]
df <- scale(df, center=T, scale=T)

imputed_df <- phenix(as.matrix(df), as.matrix(K))
write.csv(imputed_df$imp, file="imputed_phenotypes.csv", quote=F)
