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

## Now analysing differential gene expression in fat body
Load file
```
read.counts <- read.table("gene_count_mq10.txt", header = TRUE) ## read counts from feature counts following STAR mapping
row.names(read.counts) <- read.counts$Geneid
read.counts <- read.counts[ , -c(1:6)]
read.counts <- read.counts[ , -c(13:24)]
sample_info.edger <- factor(c( rep("FatBody_Control", 4), rep("FatBody_Oil", 4), rep("FatBody_WaspExtract", 4))) ### treatment as grouping variables
edgeR.DGElist <- DGEList(counts = read.counts, group = sample_info.edger) ### group read counts by treatment
```
removing genes that are expressed in 3 or fewer libraries
```
keep <- rowSums(cpm(edgeR.DGElist) >= 2) >= 4
edgeR.DGElist <- edgeR.DGElist[keep,]
dim(edgeR.DGElist)
```

add tissue filtering, remove larval salivary gland and male testes genes

Data from flyatlas2 data, courtesy of David Leader and Darren Obbard , values are FPKM
```
FlyAtlas2_Alltissues_Allgenes <- read.csv(file="FlyAtlas2_Alltissues_Allgenes.csv")
## tissue specific index
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-x/max(x))/(length(x)-1)
}
```
get larval tissues
```
FlyAtlas2_larvae <- FlyAtlas2_Alltissues_Allgenes[grep("^larval",colnames(FlyAtlas2_Alltissues_Allgenes))]
```
calculate tau
```
FlyAtlas2_larvae$tau <- apply(FlyAtlas2_larvae,1,tau) 
rownames(FlyAtlas2_larvae) <- FlyAtlas2_Alltissues_Allgenes$FlyBaseID
```
Only keep tissue specific genes, tau > 0.8
```
FlyAtlas2_larvae <- FlyAtlas2_larvae[FlyAtlas2_larvae$tau > 0.8,]
```
remove those without tau
```
FlyAtlas2_larvae <- FlyAtlas2_larvae[!is.na(FlyAtlas2_larvae$tau),]
```
 identify which tissue has maximal expression
```
max_tissue <- ""
for (i in 1:nrow(FlyAtlas2_larvae)){
  max_tissue[i] <- colnames(FlyAtlas2_larvae)[which(FlyAtlas2_larvae[i,] == max(FlyAtlas2_larvae[i,]))]
}
```
only keep tissue specific keep  where Salivary_gland has maximal expression
```
FlyAtlas2_larvae <- FlyAtlas2_larvae[max_tissue %in% "larval_Salivary_gland",]
```
only keep genes with FPKM above 1
```
FlyAtlas2_larvae <- FlyAtlas2_larvae[FlyAtlas2_larvae$larval_Salivary_gland > 1,] 
```
Correctly identifying genes with highest expression in larvae
```
pheatmap(FlyAtlas2_larvae[1:(ncol(FlyAtlas2_larvae)-1)]/(rowMeans(FlyAtlas2_larvae[1:(ncol(FlyAtlas2_larvae)-1)])))
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/4224e98b-5ba2-4197-8cf8-ed6d33cf7a69)

get male adult tissues
```
FlyAtlas2_adult_male <- FlyAtlas2_Alltissues_Allgenes[grep("^male",colnames(FlyAtlas2_Alltissues_Allgenes))] 
```
calculate tau
```
FlyAtlas2_adult_male$tau <- apply(FlyAtlas2_adult_male,1,tau)
rownames(FlyAtlas2_adult_male) <- FlyAtlas2_Alltissues_Allgenes$FlyBaseID
```
only keep tissue specific genes, tau > 0.8
```
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[FlyAtlas2_adult_male$tau > 0.8,]
```
remove those without tau
```
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[!is.na(FlyAtlas2_adult_male$tau),]
```
identify which tissue has maximal expression
```
max_tissue <- ""
for (i in 1:nrow(FlyAtlas2_adult_male)){
  max_tissue[i] <- colnames(FlyAtlas2_adult_male)[which(FlyAtlas2_adult_male[i,] == max(FlyAtlas2_adult_male[i,]))]
}
```
only keep tissue specific keep  where male testis has maximal expression
```
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[max_tissue %in% "male_Testis",]
```
only keep genes with FPKM above 1
```
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[FlyAtlas2_adult_male$male_Testis > 1,] 
```
Correctly identifying genes with highest expression in male testes
```
pheatmap(FlyAtlas2_adult_male[1:(ncol(FlyAtlas2_adult_male)-1)]/(rowMeans(FlyAtlas2_adult_male[1:(ncol(FlyAtlas2_adult_male)-1)])))
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/af0a7fb3-5635-4d5d-8959-8efc3ea79972)


now filter fatbody genes
```
dim(edgeR.DGElist[rownames(edgeR.DGElist) %in% rownames(FlyAtlas2_adult_male),])
edgeR.DGElist <- edgeR.DGElist[!rownames(edgeR.DGElist) %in% rownames(FlyAtlas2_adult_male),]
dim(edgeR.DGElist[rownames(edgeR.DGElist) %in% rownames(FlyAtlas2_larvae),])
edgeR.DGElist <- edgeR.DGElist[!rownames(edgeR.DGElist) %in% rownames(FlyAtlas2_larvae),]
write.table(rownames(edgeR.DGElist),file="fatbody_universe.txt",row.names = F,quote = F, col.names = F) ## write out fat body universe genes
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
```

PCA to see similarity of samples
```
Group <- edgeR.DGElist$samples[1]
Group <- as.factor(unlist(Group))
Group2 <- gsub("FatBody_Control","Control",Group)
Group2 <- gsub("FatBody_Oil","Oil",Group2)
Group2 <- gsub("FatBody_WaspExtract","Wasp+Oil",Group2)
cpm_log_forpca <- cpm_log
colnames(cpm_log_forpca) <- gsub("FatBody_Control","Control",colnames(cpm_log_forpca))
colnames(cpm_log_forpca) <- gsub("FatBody_Oil","Oil",colnames(cpm_log_forpca))
colnames(cpm_log_forpca) <- gsub("FatBody_WaspExtract","Wasp+Oil",colnames(cpm_log_forpca))
pca <- prcomp(t(cpm_log_forpca), scale. = TRUE)
```
plot of pca with groups in ellipses
```
fatbodypca <- fviz_pca_ind(pca,
                           col.ind = Group2, # color by groups
                           palette = c("grey",  "pink","turquoise"),
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "Groups",
                           repel = TRUE,
                           title = "Fat body - PCA"
)
fatbodypca
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/c57ef963-46ad-488d-9632-8d4d43868e56)
screeplot to see how many informative dimensions in data
```
fatbodyscree <- fviz_screeplot(pca, ncp=10, title = "Fat body - Scree plot")
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/b71b3b71-555d-42b7-b300-3c7cb9947384)

model for differential expression analysis, basically a simple comparison of three groups
```
mm <- model.matrix(~0+edgeR.DGElist$samples$group, data = edgeR.DGElist$samples)
colnames(mm) <- levels(edgeR.DGElist$samples$group)
y <- voom(edgeR.DGElist, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))
```

Looking at the contrast that shows difference between wasp+oil vs oil
```
tmp.fatbody.waspvoil <- contrasts.fit(fit, contrast=c(0,-1,1)) # Directly test second coefficient
tmp.fatbody.waspvoil <- eBayes(tmp.fatbody.waspvoil)
top.table <- topTable(tmp.fatbody.waspvoil, sort.by = "P", n = Inf)
all.top.table.fatbody.waspvoil <- topTable(tmp.fatbody.waspvoil, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05))
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(sig_genes, file = "Fatbody Wasp V. Oil sig genes.csv")
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes)) 
```

looking at the contrast that shows difference between wasp+oil vs control
```
tmp.fatbody.waspvcontrol <- contrasts.fit(fit, contrast=c(-1,0,1)) # Directly test second coefficient
tmp.fatbody.waspvcontrol <- eBayes(tmp.fatbody.waspvcontrol)
top.table <- topTable(tmp.fatbody.waspvcontrol, sort.by = "P", n = Inf)
all.top.table.fatbody.waspvcontrol <- topTable(tmp.fatbody.waspvcontrol, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05))
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(sig_genes, file = "Fatbody Wasp V. Control sig genes.csv")
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes)) 
colnames(mat_DGEgenes) <- sub("FatBody_","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("WaspExtract","Wasp",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_"," ",colnames(mat_DGEgenes))
symbols <- mapIds(org.Dm.eg.db, as.character(rownames(mat_DGEgenes)), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
#symbols["FBgn0265375"] <- "lncRNA:CR44316"
symbols["FBgn0265577"] <- "lincRNA-IBIN"
#Change names of serine proteases to nomeculature proposed by Cao 2018 "Building a platform for predicting functions of serine protease-related proteins in Drosophila melanogaster and other insects". 
Serine_Proteases_Nomenculature <- read.csv(file="Serine_Proteases_Nomenculature.csv")
for (i in 1:length(symbols)){
  if (names(symbols[i]) %in% Serine_Proteases_Nomenculature$FB_ID){
    symbols[i] <- Serine_Proteases_Nomenculature[Serine_Proteases_Nomenculature$FB_ID %in% names(symbols[i]),]$Gene
  }
}
Serine_Proteases_Nomenculature <- Serine_Proteases_Nomenculature[order(Serine_Proteases_Nomenculature$FB_ID),]
rownames(mat_DGEgenes) <- symbols
df <- data.frame(c(rep("Control",4),rep("Oil",4),rep("Wasp",4)))
rownames(df) <- colnames(mat_DGEgenes)
colnames(df) <- "group"
mycolors <- c("grey","pink","turquoise") ## define color of annotation bars
names(mycolors) <- c("Control","Oil","Wasp")
mycolors <- list(group = mycolors)
heatmap_fatbody_waspvcontrol <- pheatmap(mat_DGEgenes,treeheight_row = 0, clustering_method = "average", show_colnames = F,legend_labels = "Log",annotation_col = df,annotation_legend = F,annotation_names_col = F,border_color=NA,annotation_colors=mycolors,labels_row = as.expression(lapply(rownames(mat_DGEgenes),function(x) bquote(italic(.(x))))))
heatmap_fatbody_waspvcontrol <- as.ggplot(heatmap_fatbody_waspvcontrol) 
heatmap_fatbody_waspvcontrol <- heatmap_fatbody_waspvcontrol +
  #theme(plot.margin = unit(c(0, 0, 5, 0), "lines")) +
  scale_x_continuous(labels=c("", "Unchallenged", "Oil", "Wasp homogenate", ""), breaks=c(0.00001,0.14,0.39,0.64,0.99999)) +
  theme(axis.text.x = element_text(size=12))+
  coord_cartesian(clip = 'off') +
  theme(plot.margin = margin(0, 0, 0.1, 0, "cm"))+
  theme(legend.position = "bottom") +
  geom_text(x=0.91, y=0.9, label="Multiples of mean \nCPM per gene")
heatmap_fatbody_waspvcontrol
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/bd2749d4-da2d-4c41-a488-dd6db0707c83)

Looking at the contrast that shows difference between oil vs control
```
tmp.fatbody.oilvcontrol <- contrasts.fit(fit, contrast=c(-1,1,0)) # Directly test second coefficient
tmp.fatbody.oilvcontrol <- eBayes(tmp.fatbody.oilvcontrol)
top.table <- topTable(tmp.fatbody.oilvcontrol, sort.by = "P", n = Inf)
all.top.table.fatbody.oilvcontrol <- topTable(tmp.fatbody.oilvcontrol, sort.by = "none", n = Inf,p.value=1,lfc=0)
length(which(top.table$adj.P.Val < 0.05))
sig_genes <- subset(top.table, top.table$adj.P.Val < 0.05)
write.csv(sig_genes, file = "Fatbody Oil V. Control sig genes.csv")
DGEgenes <- rownames(subset(top.table, top.table$adj.P.Val < 0.05))
mat_DGEgenes <- cpm_nolog_relative[DGEgenes, ]
colnames(mat_DGEgenes) <- sub(".out.bam","",colnames(mat_DGEgenes))
colnames(mat_DGEgenes) <- sub("_oil","",colnames(mat_DGEgenes)) 
```

SIGNIFICANT GENE LIST COMPARISIONS
```
Fatbody_WaspV.Oil <- read.csv(file = "Fatbody Wasp V. Oil sig genes.csv")
Fatbody_WaspV.Control <- read.csv(file = "Fatbody Wasp V. Control sig genes.csv")
Fatbody_OilV.Control <- read.csv(file = "Fatbody Oil V. Control sig genes.csv")
```
get venn diagrams for number of DE genes unique to wasp+oil v control and oil v control
```
df.venn <- data.frame(x = c(0.866, -0.866),y = c(-0.5, -0.5),labels = c('A', 'B'))
venn.fatbody.ggplot  <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed(clip="off") +
  theme_void() +
  theme(legend.position = "none") +
  annotate("text", x=-1.5, y=-0.5, label= paste(intToUtf8(8593),29)) +
  annotate("text", x=1.5, y=-0.5, label= 0) +
  annotate("text", x=-0.1, y=-0.5, label= 0) +
  annotate("text", x=-1.5, y=1.3, label= "Wasp homogenate", col = "turquoise",fontface =2) +
  annotate("text", x=1.5, y=1.3, label= "Oil", col = "pink",fontface =2) +
  annotate("text", x=0, y=2.1, label= "Fat body", fontface =2, size=4.5)
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/e8033ae3-5d5b-458e-a49a-cf112728f864)

heatscatter for comparing logFC between wasp+oil vs control and oil vs control. 
```
heatscatter.fatbody.pryr %<a-% {
  heatscatter(all.top.table.fatbody.waspvcontrol$logFC,all.top.table.fatbody.oilvcontrol$logFC,ylim=c(-6,6),xlim=c(-6,6),ylab =expression('Oil v. Control log'[2]*'(FC)'),xlab =expression('Wasp homogenate v. Control log'[2]*'(FC)'), main = "")
  abline(h=0,lty=2)
  abline(v=0,lty=2)
  abline(0,1,lty=2)
  text(4.5,-0.6,paste("Wasp",intToUtf8(8593)),cex=0.9)
  text(4,4.7,paste("Injury",intToUtf8(8593)),cex=0.9,srt = 45)
  text(-4.5,0.6,paste("Wasp",intToUtf8(8595)),cex=0.9)
  text(-4,-4.7,paste("Injury",intToUtf8(8595)),cex=0.9,srt = 45)
  text(0,7.2,"Fat body",cex=1.2)
}
heatscatter.fatbody.pryr.grob <- as.ggplot(~heatscatter.fatbody.pryr)
heatscatter.fatbody.pryr.grob
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/5e2b7d75-1072-442a-ba60-867b25eef3fe)

FACET PLOTS FAT BODY
```
figure_fatbody <- multi_panel_figure(columns = 2, rows = 1)
figure_fatbody  %<>% 
  fill_panel(venn.fatbody.ggplot) %<>%
  fill_panel(heatmap_fatbody_waspvcontrol)
```

FACET PLOTS FAT BODY and HEMOCYTE
```
figure_hemocyte_fatbody <- multi_panel_figure(columns = 20, rows = 4,height = 210,width = 260)
figure_hemocyte_fatbody  %<>% 
  fill_panel(venn.hemocyte.ggplot, row = 1, column = 1:3) %<>%
  fill_panel(venn.fatbody.ggplot, row = 1, column = 4:6) %<>%
  fill_panel(heatmap_fatbody_waspvcontrol, row = 2:4, column = 1:6)  %<>%
  fill_panel(heatscatter.hemocyte.pryr.grob, row = 1:2, column = 13:20) #%<>%
#fill_panel(heatscatter.fatbody.pryr.grob, row = 3:4, column = 7:10)
figure_hemocyte_fatbody

#ggsave(file="plots/hemocyte_fatbody_grid_plot.pdf", width = 260, height = 210, units = "mm", plot = figure_hemocyte_fatbody,device=cairo_pdf, family = "Arial")

#cairo_pdf("plots/venn.hemocyte.ggplot.pdf",width=3,height=1.9, family = "Arial")
venn.hemocyte.ggplot
#dev.off()

#cairo_pdf("plots/venn.fatbody.ggplot.pdf",width=3,height=1.9, family = "Arial")
venn.fatbody.ggplot
#dev.off()

#cairo_pdf("plots/heatmap_fatbody_waspvcontrol.pdf",width=3,height=1.9, family = "Arial")
heatmap_fatbody_waspvcontrol
#dev.off()

#cairo_pdf("plots/heatmap_fatbody_waspvcontrol.pdf",width=3,height=1.9, family = "Arial")
heatmap_fatbody_waspvcontrol
#dev.off()
```

volcano plots for differentially expressed genes in each category
```
#pdf(file="plots/Volcano.pdf",width=8.27, height=11.69)
par(mfrow=c(2,2))
##logfc vs -log10(pvalue)all genes sig genes

#with(all.top.table.hemocyte.waspvoil, plot(all.top.table.hemocyte.waspvoil$logFC, -log10(all.top.table.hemocyte.waspvoil$P.Value), pch=20, main="hemocyte waspvoil", xlim=c(-2.5,2),col=ifelse(all.top.table.hemocyte.waspvoil$adj.P.Val < 0.05, "red", "black")))
with(all.top.table.hemocyte.waspvcontrol, plot(all.top.table.hemocyte.waspvcontrol$logFC, -log10(all.top.table.hemocyte.waspvcontrol$P.Value), pch=20, xlab="Log(2)FC",ylab="Log(10)P-value",main="Hemocyte - Wasp homogenate v. Unchallenged", xlim=c(-2.5,2),col=ifelse(all.top.table.hemocyte.waspvcontrol$adj.P.Val < 0.05, "red", "black")))
with(all.top.table.hemocyte.oilvcontrol, plot(all.top.table.hemocyte.oilvcontrol$logFC, -log10(all.top.table.hemocyte.oilvcontrol$P.Value), pch=20,xlab="Log(2)FC",ylab="Log(10)P-value", main="Hemocyte - Oil v. Unchallenged", xlim=c(-2.5,2),col=ifelse(all.top.table.hemocyte.oilvcontrol$adj.P.Val < 0.05, "red", "black")))

#with(all.top.table.fatbody.waspvoil, plot(all.top.table.fatbody.waspvoil$logFC, -log10(all.top.table.fatbody.waspvoil$P.Value), pch=20, main="fatbody waspvoil", xlim=c(-2.5,2),col=ifelse(all.top.table.fatbody.waspvoil$adj.P.Val < 0.05, "red", "black")))
with(all.top.table.fatbody.waspvcontrol, plot(all.top.table.fatbody.waspvcontrol$logFC, -log10(all.top.table.fatbody.waspvcontrol$P.Value), pch=20,xlab="Log(2)FC",ylab="Log(10)P-value", main="Fat body - Wasp homogenate v. Unchallenged", xlim=c(-2.5,2),col=ifelse(all.top.table.fatbody.waspvcontrol$adj.P.Val < 0.05, "red", "black")))
with(all.top.table.fatbody.oilvcontrol, plot(all.top.table.fatbody.oilvcontrol$logFC, -log10(all.top.table.fatbody.oilvcontrol$P.Value), pch=20, xlab="Log(2)FC",ylab="Log(10)P-value", main="Fat body - Oil v. Unchallenged", xlim=c(-2.5,2),col=ifelse(all.top.table.fatbody.oilvcontrol$adj.P.Val < 0.05, "red", "black")))
#dev.off()
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/ad46f0d3-f7e3-44f2-b8a8-d39435132480)

library metrics - number of uniquely mapped exonic reads
```
allmetrics <- read.csv(file = "All_metrics_with_labels.csv")
allmetrics$ExonicReads <- as.numeric(gsub(",","",gsub(" ","",gsub("\\(.*","", allmetrics$Exonic))))
allmetrics <- allmetrics[order(allmetrics$Tissue,allmetrics$Treatment),]
allmetrics$Treatment <- gsub("Wasp Extract","Wasp homogenate",allmetrics$Treatment)
plot.p <- ggplot(data = allmetrics, aes(x=factor(Library_name,levels = Library_name),y=ExonicReads,fill=Treatment))+
  geom_bar(stat="identity")+
  facet_wrap(~Tissue, scales="free")+ 
  scale_fill_manual(values = c("grey","pink","turquoise"))+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  xlab("Samples")+
  ylab("Number of uniquely mapped exonic reads")
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/888d0668-92ca-45b7-b64c-7ae8c442cd6f)

test to see if treatment of tissue had an impact on number of reads
```
anova(lm(data = allmetrics, ExonicReads ~ Treatment * Tissue))
table.p <- ggtexttable(anova(lm(data = allmetrics, ExonicReads ~ Treatment * Tissue)))
table.p <- table.p %>% tab_add_title("anova(lm(Exonic Reads ~ Treatment * Tissue))", hjust = -0.8)
nummappedreads <- ggarrange(plot.p, table.p, ncol = 1, heights = c(1, 0.5))

#pdf(file="plots/exonic reads.pdf",height=5.83, width=8.27)
nummappedreads
#dev.off()
```
![image](https://github.com/arunkumarramesh/PAMP/assets/23363383/029fc420-7d18-4215-a396-0f13e3f8f7b4)


