library(ggplot2)
library(plyr)
library(tidyr)
library(reshape2)
library(ggrepel)
library(cowplot)

#set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This changes the path director to the folder where the file is saved
list.files()
data = read.table("Fat_body_qPCR_Confirmation.csv", sep = ",", dec = ".", header = TRUE)
head(data)
str(data)


#Create data summaries
data_summary = ddply(data,.(Plate,Gene,Treatment, Biological_Replica, Target_Name, Target_Type),  summarise,
                     mean_Ct = mean(Ct))

data_summary = ddply(data,.(Plate,Target_Name, Treatment, Biological_Replica),  summarise,
                     mean_Ct = mean(Ct))

data_wide <- spread(data_summary, Target_Name, mean_Ct)

data_Delta_Ct = data_wide

data_Delta_Ct$CG33461 = (data_wide$Rpl-data_wide$CG33461)
data_Delta_Ct$CG6788 = (data_wide$Rpl-data_wide$CG6788)
data_Delta_Ct$IM1 = (data_wide$Rpl-data_wide$IM1)
data_Delta_Ct$IM2 = (data_wide$Rpl-data_wide$IM1)
data_Delta_Ct$IM3 = (data_wide$Rpl-data_wide$IM3)
data_Delta_Ct$'lectin-24A' = (data_wide$Rpl-data_wide$'lectin-24A')
data_Delta_Ct$CG33462 = (data_wide$Rpl-data_wide$CG33462)
data_Delta_Ct$Tep1 = (data_wide$Rpl-data_wide$Tep1)
data_Delta_Ct$CG10764 = (data_wide$Rpl-data_wide$CG10764)
data_Delta_Ct$CG11313 = (data_wide$Rpl-data_wide$CG11313)
data_Delta_Ct$CG34436 = (data_wide$Rpl-data_wide$CG34436)
data_Delta_Ct$CG43124 = (data_wide$Rpl-data_wide$CG43124)
data_Delta_Ct$CG30002 = (data_wide$Rpl-data_wide$CG30002)

data_Delta_Ct = data_Delta_Ct[,!names(data_Delta_Ct) %in% c('Rpl', 'Plate')]

mdata_Delta_Ct = melt(data_Delta_Ct, id.vars = c('Treatment', 'Biological_Replica'))
colnames(mdata_Delta_Ct) = c('Treatment', 'Biological_Replica', 'Gene', 'DeltaCt')

mdata_Delta_Ct = mdata_Delta_Ct[!mdata_Delta_Ct$DeltaCt=='NA',]
mdata_Delta_Ct = na.omit(mdata_Delta_Ct)

se = function(x) sd(x)/sqrt(length(x))
mdata_Delta_Ct_summary = ddply(mdata_Delta_Ct,.(Gene,Treatment),  summarise,
                          mean_DeltaCt = mean(DeltaCt),
                          CI_Lower = -se(DeltaCt)+mean(DeltaCt),
                          CI_Upper = se(DeltaCt)+mean(DeltaCt))
mdata_Delta_Ct_summary$Treatment <- gsub("Wasp_Extract","Wasp + Oil",mdata_Delta_Ct_summary$Treatment)

#Plot Delta Ct
p = ggplot(mdata_Delta_Ct_summary, aes(x=Treatment, y=mean_DeltaCt))+
  geom_bar(stat='identity', aes(fill=Treatment))+
  geom_errorbar(aes(ymin=CI_Lower, ymax=CI_Upper), width=.4)+
  scale_fill_manual(values = c('grey20', 'grey40', 'grey60'))+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  xlab(NULL)+
  ylab(expression(paste(Delta, 'Ct')))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0, color='grey30',size=12),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.text.y = element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"))+
  facet_grid(.~Gene) +
  theme(strip.text = element_text(face = "italic"))
p

#Plot DeltaDelta Ct - just averages
DeltaDelta_Ct_summary = spread(mdata_Delta_Ct_summary[,1:3], Treatment, mean_DeltaCt)
colnames(DeltaDelta_Ct_summary)[4] <- "Wasp_Extract"
DeltaDelta_Ct_summary$Oil_FC = 2^-(DeltaDelta_Ct_summary$Control-DeltaDelta_Ct_summary$Oil)
DeltaDelta_Ct_summary$Wasp_Extract_FC = 2^-(DeltaDelta_Ct_summary$Control-DeltaDelta_Ct_summary$Wasp_Extract)
mDeltaDelta_Ct_summary = melt(DeltaDelta_Ct_summary[,c(1,5,6)], id.vars = c('Gene'))
colnames(mDeltaDelta_Ct_summary) = c('Gene', 'Treatment', 'FC')

p1 = ggplot(mDeltaDelta_Ct_summary, aes(x=Treatment, y=log2(FC)))+
  geom_bar(stat='identity', aes(fill=Treatment))+
  scale_fill_manual(values = c('grey40', 'grey20'))+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  xlab(NULL)+
#  ylab(expression(paste(Delta, 'Ct')))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0, color='grey30',size=12),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.text.y = element_text(size=15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"))+
  facet_grid(.~Gene)
p1


p2 = ggplot(DeltaDelta_Ct_summary, aes(x=log2(Wasp_Extract_FC), y=log2(Oil_FC), label = Gene))+
  geom_point()+
  geom_text_repel(fontface = "italic")+
  theme_bw()+
  scale_y_continuous(limits = c(-1, 9))+
  scale_x_continuous(limits = c(0, 9))+
  geom_abline(intercept = 0, slope = 1)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  ggtitle("qPCR") +
  xlab('log2(Wasp + Oil FC)')+
  ylab('log2(Oil FC)')+
  theme(axis.line=element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.y = element_text(angle = 90,vjust = 0, color='grey30',size=15),
        axis.text.x = element_text(angle = 90,vjust = 0, color='grey30',size=15),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.title.x = element_text(colour = "grey30", face='bold',size=15),
        panel.border = element_blank())
p2

# Plot RNAseq
list.files()
RNAseq_data = read.table("logfc_fatbody_RNAseq_Data.csv", sep = ",", dec = ".", header = TRUE)
FB_IDs = read.table("genes_tested_by_qPCR_IDs.csv", sep = ",", dec = ".", header = TRUE)
head(RNAseq_data)
str(RNAseq_data)

#Filter data
RNAseq_data = RNAseq_data[RNAseq_data$FB_ID %in% FB_IDs$FB_ID,]
RNAseq_data = merge(RNAseq_data, FB_IDs, by='FB_ID')


p3 = ggplot(RNAseq_data, aes(x=waspvcontrol, y=waspvoil, label = Gene))+
  geom_point()+
  geom_text_repel(fontface = "italic")+
  theme_bw()+
  scale_y_continuous(limits = c(-1, 9))+
  scale_x_continuous(limits = c(0, 9))+
  geom_abline(intercept = 0, slope = 1)+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)+
  xlab('log2(Wasp + Oil FC)')+
  ylab('log2(Oil FC)')+
  ggtitle("RNA-seq")+
  theme(axis.line=element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text.y = element_text(angle = 90,vjust = 0, color='grey30',size=15),
        axis.text.x = element_text(angle = 90,vjust = 0, color='grey30',size=15),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.title.x = element_text(colour = "grey30", face='bold',size=15),
        panel.border = element_blank())
p3

p23 <- plot_grid(p2,p3,labels = c("b","c"))
plot_grid(p,p23,ncol=1,labels = c("a",""))

pdf(file="Fat_Body_qPCR.pdf",height=6,width=12)
p
dev.off()
