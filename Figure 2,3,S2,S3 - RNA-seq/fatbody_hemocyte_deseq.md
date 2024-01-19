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
                             legend.title = "Groups",
                             repel = TRUE,
                             title = "Hemocytes - PCA"
)
Hemocytespca
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/3416a3a6-f6f4-42c4-b6ce-8b591586b9fa)

screeplot to see how many informative dimensions in data
```
Hemocytesscree <- fviz_screeplot(pca, ncp=10,title = "Hemocytes - Scree plot")
Hemocytesscree
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/3b78e5ea-c12f-40bc-acc9-edae1e00d099)

model for differential expression analysis, basically a simple comparison of three groups
```
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
```
Check that the mean-variance plot looks good
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/602a553b-0992-4f10-be06-60f333e90cc1)


```
fit model
fit <- lmFit(y, mm)
head(coef(fit))
```

looking at the contrast that shows difference between wasp+oil vs oil
```
tmp.hemocyte.waspvoil <- contrasts.fit(fit, contrast=c(0,-1,1)) 
tmp.hemocyte.waspvoil <- eBayes(tmp.hemocyte.waspvoil) ## calculate DE stats
top.table <- topTable(tmp.hemocyte.waspvoil, sort.by = "P", n = Inf) ## sort by most significantly DE genes
all.top.table.hemocyte.waspvoil <- topTable(tmp.hemocyte.waspvoil, sort.by = "none", n = Inf,p.value=1,lfc=0) ## get all genes with logFC and pvalues
length(which(top.table$adj.P.Val < 0.05)) ## how many significantly DE genes

sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05) ##  significantly DE genes only
write.csv(sig_genes, file = "Hemocyte Wasp V. Oil sig genes.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes)) 
```
Heatmap of differentially expressed genes for wasp+oil versus oil
```
pheatmap(mat_DGEgenes)
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/8930c6b8-e8c2-4d32-8867-2816794bcdd3)


now doing differential expression tests for wasp+oil versus control - effect of wasp+injury, same methods as before
```
tmp.hemocyte.waspvcontrol <- contrasts.fit(fit, contrast=c(-1,0,1)) 
tmp.hemocyte.waspvcontrol  <- eBayes(tmp.hemocyte.waspvcontrol )
top.table <- topTable(tmp.hemocyte.waspvcontrol , sort.by = "P", n = Inf)
all.top.table.hemocyte.waspvcontrol <- topTable(tmp.hemocyte.waspvcontrol , sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05)) ### significantly DE genes

sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(sig_genes, file = "Hemocyte Wasp V. Control sig genes.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes))
```
Heatmap of differentially expressed genes for wasp+oil versus control
```
pheatmap(mat_DGEgenes)
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/d20a7ee5-5c2e-4d3a-b359-e74228dab66b)


now doing differential expression tests for oil versus control - effect of injury, same methods as before
```
tmp.hemocyte.oilvcontrol  <- contrasts.fit(fit, contrast=c(-1,1,0)) 
tmp.hemocyte.oilvcontrol <- eBayes(tmp.hemocyte.oilvcontrol)
top.table <- topTable(tmp.hemocyte.oilvcontrol, sort.by = "P", n = Inf)
all.top.table.hemocyte.oilvcontrol <- topTable(tmp.hemocyte.oilvcontrol, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05))

sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(sig_genes, file = "Hemocyte Oil V. Control sig genes.csv")

DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes)) 
```
Heatmap of differentially expressed genes for oil versus control
```
pheatmap(mat_DGEgenes)
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/a8b469ef-75b4-4ce6-b765-2a0d17178094)


SIGNIFICANT GENE LIST COMPARISIONS
```
Hemocyte_WaspV.Oil <- read.csv(file = "Hemocyte Wasp V. Oil sig genes.csv")
Hemocyte_WaspV.Control <- read.csv(file = "Hemocyte Wasp V. Control sig genes.csv")
Hemocyte_OilV.Control <- read.csv(file = "Hemocyte Oil V. Control sig genes.csv")
```
get venn diagrams for number of DE genes unique to wasp+oil v control and oil v control
```
GOVenn(Hemocyte_WaspV.Control,Hemocyte_WaspV.Oil,Hemocyte_OilV.Control, label = c("Hemocyte_WaspV.Control","Hemocyte_WaspV.Oil","Hemocyte_OilV.Control"))
GOVenn(Hemocyte_WaspV.Control,Hemocyte_OilV.Control, label = c("Hemocyte_WaspV.Control","Hemocyte_OilV.Control"))

A_only_vcntl <- GOVenn(Hemocyte_WaspV.Control,Hemocyte_OilV.Control, label = c("Hemocyte_WaspV.Control","Hemocyte_OilV.Control"),plot = F)$table$A_only
AB_vcntl <- GOVenn(Hemocyte_WaspV.Control,Hemocyte_OilV.Control, label = c("Hemocyte_WaspV.Control","Hemocyte_OilV.Control"),plot = F)$table$AB

cpm_nolog2 <- as.data.frame(cpm_nolog)
cpm_nolog2$FB_ID <- rownames(cpm_nolog2)

A_only_vcntl$FB_ID <- row.names(A_only_vcntl)
A_only_vcntl <- left_join(A_only_vcntl,cpm_nolog2,by="FB_ID")
write.csv(A_only_vcntl,file="Wasp_specific-Venn.csv",quote=F,row.names = F)

colnames(AB_vcntl)[1:2] <- c("logFC_Wasp","logFC_Oil")
AB_vcntl$FB_ID <- row.names(AB_vcntl)
AB_vcntl <- left_join(AB_vcntl,cpm_nolog2,by="FB_ID")

write.csv(AB_vcntl,file="Wasp_and_Injury-Venn.csv",quote=F,row.names = F)
```

make a custom venn diagram to make it easier to control aesthetic features, using numbers from before
```
df.venn <- data.frame(x = c(0.866, -0.866),y = c(-0.5, -0.5),labels = c('A', 'B'))
venn.hemocyte.ggplot  <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed(clip="off") +
  theme_void() +
  theme(legend.position = "none") +
  annotate("text", x=-1.5, y=-0.25, label= paste(intToUtf8(8593),1691)) +
  annotate("text", x=-1.5, y=-0.75, label= paste(intToUtf8(8595),1699)) +
  annotate("text", x=1.5, y=-0.25, label= paste(intToUtf8(8593),3)) +
  annotate("text", x=1.5, y=-0.75, label= paste(intToUtf8(8595),1)) +
  annotate("text", x=0, y=-0.25, label= paste(intToUtf8(8593),465)) +
  annotate("text", x=0, y=-0.75, label= paste(intToUtf8(8595),31)) +
  annotate("text", x=0, y=-1.2, label= "(1)") +
  annotate("text", x=-1.5, y=1.3, label= "Wasp homogenate", col = "turquoise",fontface =2) +
  annotate("text", x=1.5, y=1.3, label= "Oil", col = "pink",fontface =2) +
  annotate("text", x=0, y=2.1, label= "Hemocyte", fontface =2, size=4.5)
df.venn
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/634876d1-a3f4-483a-91b7-bd51d85b4cd4)

heatscatter for comparing logFC between wasp+oil vs control and oil vs control. 
```
heatscatter.hemocyte.pryr %<a-% {
  heatscatter(all.top.table.hemocyte.waspvcontrol$logFC,all.top.table.hemocyte.oilvcontrol$logFC,ylim=c(-6,6),xlim=c(-6,6),ylab =expression('Oil v. Control log'[2]*'(FC)'),xlab =expression('Wasp homogenate v. Control log'[2]*'(FC)'), main = "")
  abline(h=0,lty=2)
  abline(v=0,lty=2)
  abline(0,1,lty=2)
  text(4.5,-0.6,paste("Wasp",intToUtf8(8593)),cex=0.9)
  text(4,4.7,paste("Injury",intToUtf8(8593)),cex=0.9,srt = 45)
  text(-4.5,0.6,paste("Wasp",intToUtf8(8595)),cex=0.9)
  text(-4,-4.7,paste("Injury",intToUtf8(8595)),cex=0.9,srt = 45)
  text(0,7.2,"Hemocyte",cex=1.2)
}
heatscatter.hemocyte.pryr.grob <- as.ggplot(~heatscatter.hemocyte.pryr)

cairo_pdf("plots/heatscatter.hemocyte.pryr.pdf",width=3.8,height=4.1, family = "Arial")
heatscatter(all.top.table.hemocyte.waspvcontrol$logFC,all.top.table.hemocyte.oilvcontrol$logFC,ylim=c(-6,6),xlim=c(-6,6),ylab =expression('Oil v. Unchallenged log'[2]*'(FC)'),xlab =expression('Wasp homogenate v. Unchallenged log'[2]*'(FC)'), main = "")
abline(h=0,lty=2)
abline(v=0,lty=2)
abline(0,1,lty=2)
text(4.5,-0.6,paste("Wasp",intToUtf8(8593)),cex=0.9)
text(4,4.7,paste("Injury",intToUtf8(8593)),cex=0.9,srt = 45)
text(-4.5,0.6,paste("Wasp",intToUtf8(8595)),cex=0.9)
text(-4,-4.7,paste("Injury",intToUtf8(8595)),cex=0.9,srt = 45)
dev.off()
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/cef469a4-fa71-4678-8715-fca6e4e059dc)


FACET PLOTS HEMOCYTE
```
figure_hemocyte <- multi_panel_figure(columns = 2, rows = 1)
figure_hemocyte  %<>% 
  fill_panel(venn.hemocyte.ggplot) %<>%
  fill_panel(heatscatter.hemocyte.pryr.grob)
```
what are the gene ontology and KEGG categories for differential expressed genes

first for upregulated genes
```
genelist_hemocyte_upregulated <- mapIds(org.Dm.eg.db, as.character(subset(Hemocyte_WaspV.Control$X, Hemocyte_WaspV.Control$logFC > 0)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
universelist <- mapIds(org.Dm.eg.db, rownames(edgeR.DGElist), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

go.fisher.upregulated <-  goana(genelist_hemocyte_upregulated, universe = universelist, species = "Dm", prior.prob = NULL, covariate=NULL,plot = TRUE)
go.fisher.upregulated$ADJ<- p.adjust(go.fisher.upregulated$P.DE, method = "fdr")
go.fisher.upregulated.sig <- subset(go.fisher.upregulated, go.fisher.upregulated$ADJ < 0.05)
#write.csv(go.fisher.upregulated.sig,file = "Hemocyte Waspvcontrol Upregulated Go.csv")
go.fisher.upregulated.sig.plot <- head(go.fisher.upregulated.sig[order(go.fisher.upregulated.sig$ADJ),], n = 50)
go.fisher.upregulated.sig.plot <- go.fisher.upregulated.sig.plot[order(go.fisher.upregulated.sig.plot$N,decreasing=TRUE),]
go.fisher.upregulated.sig.plot$Term2 <- paste(go.fisher.upregulated.sig.plot$Term," (",go.fisher.upregulated.sig.plot$DE,"/",go.fisher.upregulated.sig.plot$N,")",sep="")

go.hemocyte.waspvcontrol.upregulated <- ggplot(aes(x=-log10(go.fisher.upregulated.sig.plot$P.DE), y=factor(go.fisher.upregulated.sig.plot$Term2, levels = go.fisher.upregulated.sig.plot$Term2), colour=go.fisher.upregulated.sig.plot$DE/go.fisher.upregulated.sig.plot$N, size=go.fisher.upregulated.sig.plot$N), data = go.fisher.upregulated.sig.plot) +
  geom_point(stat="identity",position = "identity") +
  expand_limits(x=18) +
  labs(x="-log 10 (p value)", y="GO term", colour="Proportion \noverrepresented", size="Count") +
  scale_color_viridis(option = "D")
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/7622f5c3-3a13-4fc7-bfa2-0655e744422f)

now for downregulated genes
```
genelist_hemocyte_downregulated <- mapIds(org.Dm.eg.db, as.character(subset(Hemocyte_WaspV.Control$X, Hemocyte_WaspV.Control$logFC < 0)), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
universelist <- mapIds(org.Dm.eg.db, rownames(edgeR.DGElist), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

go.fisher.downregulated <-  goana(genelist_hemocyte_downregulated, universe = universelist, species = "Dm", prior.prob = NULL, covariate=NULL,plot = TRUE)
go.fisher.downregulated$ADJ<- p.adjust(go.fisher.downregulated$P.DE, method = "fdr")
go.fisher.downregulated.sig <- subset(go.fisher.downregulated, go.fisher.downregulated$ADJ < 0.05)
#write.csv(go.fisher.downregulated.sig,file = "Hemocyte Waspvcontrol Downregulated Go.csv")

go.fisher.downregulated.sig.plot <- head(go.fisher.downregulated.sig[order(go.fisher.downregulated.sig$ADJ),], n = 50)
go.fisher.downregulated.sig.plot <- go.fisher.downregulated.sig.plot[order(go.fisher.downregulated.sig.plot$N,decreasing=TRUE),]
go.fisher.downregulated.sig.plot$Term2 <- paste(go.fisher.downregulated.sig.plot$Term," (",go.fisher.downregulated.sig.plot$DE,"/",go.fisher.downregulated.sig.plot$N,")",sep="")
go.hemocyte.waspvcontrol.downregulated <- ggplot(aes(x=-log10(go.fisher.downregulated.sig.plot$P.DE), y=factor(go.fisher.downregulated.sig.plot$Term2, levels = go.fisher.downregulated.sig.plot$Term2), colour=go.fisher.downregulated.sig.plot$DE/go.fisher.downregulated.sig.plot$N, size=go.fisher.downregulated.sig.plot$N), data = go.fisher.downregulated.sig.plot) +
  geom_point(stat="identity",position = "identity") +
  expand_limits(x=10) +
  labs(x="-log 10 (p value)", y="GO term", colour="Proportion \noverrepresented", size="Count") +
  scale_color_viridis(option = "D")
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/6ed06478-6663-4743-a5d2-ddefaf35f0f5)


