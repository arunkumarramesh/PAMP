# RNA-seq anaylysis pipeline

Shift to R. Load packages and set working directory

```
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(magrittr)
library(edgeR)
library(locfit)
library(factoextra)
library(tidyr)
library(limma)
library(pheatmap)
library(GOplot)
library(car) 
library(ggplot2)
library(gridExtra)
library(GO.db)
library(org.Dm.eg.db)
library(TimeSeriesExperiment)
library(LSD)
library(grid)
library(pryr)
library(viridis)
library(ggplotify)
library(ggforce)
library(multipanelfigure)
library(cowplot)
library(ggpubr)
library(dplyr)
library(plyr)
library(reshape2)

```

# First analysing differential expression in hemocytes
Load file
```
read.counts <- read.table("gene_count_mq10.txt", header = TRUE) ## read counts from feature counts following STAR mapping
row.names(read.counts) <- read.counts$Geneid
read.counts <- read.counts[ , -c(1:6)]
read.counts <- read.counts[ , -c(1:12)]
sample_info.edger <- factor(c( rep("Hemocytes_Control", 4), rep("Hemocytes_Oil", 4), rep("Hemocytes_WaspExtract", 4))) ### treatment as grouping variables
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) ### group read counts by treatment
```
removing genes that are expressed in 3 or fewer libraries
```
keep <- rowSums( cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
dim(edgeR.DGElist)

write.table(rownames(edgeR.DGElist),file="hemocyte_universe.txt",row.names = F,quote = F, col.names = F) ## write out fat body universe genes
```

normalisation
```
edgeR.DGElist$samples$lib.size <- colSums(edgeR.DGElist$counts)
head(edgeR.DGElist$samples)
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")
edgeR.DGElist$samples
```

log fold change & relative fold change 
```
cpm_log <- cpm(edgeR.DGElist, log = TRUE)
cpm_nolog <- cpm(edgeR.DGElist, log = FALSE)
colnames(cpm_log) <- sub(".out.bam","",colnames(cpm_log))
colnames(cpm_log) <- sub("_oil","",colnames(cpm_log))

cpm_nolog_relative <- cpm_nolog/rowMeans(cpm_nolog)
colnames(cpm_nolog_relative) <- sub(".out.bam","",colnames(cpm_nolog_relative))
colnames(cpm_nolog_relative) <- sub("_oil","",colnames(cpm_nolog_relative))

colnames(cpm_nolog) <- sub(".out.bam","",colnames(cpm_nolog))
colnames(cpm_nolog) <- sub("_oil","",colnames(cpm_nolog))
```

PCA to see similarity of samples
```
Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
## modify group names to be more meaningful
Group2 <- gsub("Hemocytes_Control","Control",Group)
Group2 <- gsub("Hemocytes_Oil","Oil",Group2)
Group2 <- gsub("Hemocytes_WaspExtract","Wasp+Oil",Group2)
cpm_log_forpca <- cpm_log
## modify treatment names to be more meaningful
colnames(cpm_log_forpca) <- gsub("Hemocytes_Control","Control",colnames(cpm_log_forpca))
colnames(cpm_log_forpca) <- gsub("Hemocytes_Oil","Oil",colnames(cpm_log_forpca))
colnames(cpm_log_forpca) <- gsub("Hemocytes_WaspExtract","Wasp+Oil",colnames(cpm_log_forpca))
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE) ## do pca
## plot of pca with groups in ellipses
Hemocytespca <- fviz_pca_ind(pca,
                             col.ind = Group2, # color by groups
                             palette = c("grey",  "pink","turquoise"),
                             addEllipses = TRUE, # Concentration ellipses
                             ellipse.type = "confidence",
                             legend.title = "Groups",
                             repel = TRUE,
                             title = "Hemocytes - PCA"
)
```
screeplot to see how many informative dimensions in data
```
Hemocytesscree <- fviz_screeplot(pca, ncp=10,title = "Hemocytes - Scree plot")
```




