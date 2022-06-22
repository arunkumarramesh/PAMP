setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## modify single cell file to meet input requirements of CIBERSORT. Only keeping highly variable genes
#GSE148826_intergrated_rna_count_matrix2 <- read.table(file="GSE148826_intergrated_rna_count_matrix2.csv",header = T)
#colnames(GSE148826_intergrated_rna_count_matrix2)[1] <- "Gene"
#GSE148826_intergrated_rna_count_matrix2 <- t(GSE148826_intergrated_rna_count_matrix2)
#GSE148826_intergrated_rna_count_matrix2 <- as.data.frame(GSE148826_intergrated_rna_count_matrix2)
#allnames <- GSE148826_intergrated_rna_count_matrix2[1,]
#allnames <- gsub("PLASM1","plasmatocytes",as.character(allnames[1,]))
#allnames <- gsub("PLASM2","plasmatocytes",allnames)
#allnames <- gsub("AMP","plasmatocytes",allnames)
#allnames <- gsub("MET","plasmatocytes",allnames)
#allnames <- gsub("LAM1","lamellocytes",allnames)
#allnames <- gsub("LAM2","lamellocytes",allnames)
#allnames <- gsub("LAM3","lamellocytes",allnames)
#hvg <- read.table(file = "hvg")
#GSE148826_intergrated_rna_count_matrix3 <- GSE148826_intergrated_rna_count_matrix2[rownames(GSE148826_intergrated_rna_count_matrix2) %in% hvg$V1,]
#GSE148826_intergrated_rna_count_matrix3 <- rbind(allnames,GSE148826_intergrated_rna_count_matrix3)
#write.table(GSE148826_intergrated_rna_count_matrix3, file = "GSE148826_intergrated_rna_count_matrix3.txt", quote = F, col.names = F, sep = "\t")


library(ggplot2)
library(dplyr)
library(reshape2)

## read in deconvolution results from CIBERSORT website (https://cibersortx.stanford.edu/)
CIBERSORT <- read.csv(file="CIBERSORTx_Job11_Results.csv")
CIBERSORT$group <- c("Unchallenged","Unchallenged","Unchallenged","Unchallenged","Oil","Oil","Oil","Oil","Wasp homogenate","Wasp homogenate","Wasp homogenate","Wasp homogenate")
CIBERSORT <- CIBERSORT[c(2:9,13)]
CIBERSORT <- CIBERSORT[c(9,1:8)]

## get mean proportion of each cell type for the three conditions
CIBERSORT_mean <- CIBERSORT %>%
  dplyr::group_by(group) %>%
  dplyr::summarise_all(list(mean = mean))
colnames(CIBERSORT_mean) <- colnames(CIBERSORT)

## df manipulations
CIBERSORT <- melt(CIBERSORT)
CIBERSORT_mean <- melt(CIBERSORT_mean)
CIBERSORT_mean$variable <- factor(CIBERSORT_mean$variable,levels=c("CC","AMP","MET","PLASM2","PLASM1","LAM1","LAM2","LAM3"))
CIBERSORT$variable <- factor(CIBERSORT$variable,levels=c("CC","AMP","MET","PLASM2","PLASM1","LAM1","LAM2","LAM3"))
CIBERSORT_mean$group <- factor(CIBERSORT_mean$group,levels=c("Unchallenged","Oil","Wasp homogenate"))
CIBERSORT$group <- factor(CIBERSORT$group,levels=c("Unchallenged","Oil","Wasp homogenate"))

## get total proportion of each cell type for the three conditions
CIBERSORT_mean[CIBERSORT_mean$variable %in% c("LAM1","LAM2","LAM3"),] %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(sum(value))

## plot proportion of Lamellocyte clusters in RNA-seq data
pdf(file="deconvolution.pdf", width=4, height= 2.7)
ggplot(data=CIBERSORT_mean[CIBERSORT_mean$variable %in% c("LAM1","LAM2","LAM3"),],aes(x=variable,y=value)) +
  geom_col(fill="grey") +
  geom_jitter(data=CIBERSORT[CIBERSORT$variable %in% c("LAM1","LAM2","LAM3"),],aes(x=variable,y=value),width=0.25,height=0,size=1, alpha=0.4) +
  facet_grid(~group) +
  ylab("Proportion of cells") +
  xlab("Lamellocyte cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0))
dev.off()

