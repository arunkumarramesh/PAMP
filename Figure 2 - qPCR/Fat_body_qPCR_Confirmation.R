library(ggplot2)
library(plyr)
library(tidyr)
library(reshape2)

#set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #This changes the path director to the folder where the file is saved
data2 = read.table("Fat_body_qPCR_Confirmation.csv", sep = ",", dec = ".", header = TRUE)
gene_names=read.csv("gene_names.csv")
data=merge(data2,gene_names)
dim(data)
dim(data2)
ind=which(data$Gene==data$Target_Name)
data$Gene=data$Gene2
data$Target_Name[ind]=data$Gene2[ind]

data$ID = paste(data$Sample_Name, data$Target_Name)
data$Target_Type=as.factor(data$Target_Type)

data$Gene=as.factor(data$Gene)

#Check is technical replicas are ok
data_wide = spread(data, key=c(Technical_Replica), value=c(Ct))
data_wide$DeltaCt = abs(data_wide$'1'-data_wide$'2')
hist(data_wide$DeltaCt)#A few replicas have Delta Cts between replicas>1.
data_wide[data_wide$DeltaCt>1,]

#Remove samples where DeltaCt between technical replicas >1
data=data[!data$ID %in% data_wide$ID[data_wide$DeltaCt>1],]

#Check genes that, after filtering, have <3 biological replicas per treatment
table(data$Treatment, data$Target_Name, data$Technical_Replica)
genes_to_drop=c('CG10764', 'CG6788')
data=data[!data$Gene %in% genes_to_drop,]
data=droplevels(data)

#Create data summaries to calculate average of Ct between technical replicas
data_summary = ddply(data,.(Plate,Target_Name, Treatment, Biological_Replica),  summarise,
                     mean_Ct = mean(Ct))

data_wide <- spread(data_summary, Target_Name, mean_Ct) #The transformation of data_summary into wide format leaves many NAs because there were a total of 4 qPCR plates and in each there was a Rpl32 control, to each gene expression is normalised.

#Create data frame with DeltaCt between target gene and Rpl
data_Delta_Ct = data_wide
for(i in 4:ncol(data_Delta_Ct)){
  data_Delta_Ct[,i] = data_wide$Rpl-data_wide[,i]
}


data_Delta_Ct = data_Delta_Ct[,!names(data_Delta_Ct) %in% c('Rpl', 'Plate')] # Removes columns for Rpl and Plate

mdata_Delta_Ct = melt(data_Delta_Ct, id.vars = c('Treatment', 'Biological_Replica'))
colnames(mdata_Delta_Ct) = c('Treatment', 'Biological_Replica', 'Gene', 'DeltaCt')

mdata_Delta_Ct = mdata_Delta_Ct[!mdata_Delta_Ct$DeltaCt=='NA',]
mdata_Delta_Ct = na.omit(mdata_Delta_Ct)

se = function(x) sd(x)/sqrt(length(x))
mdata_Delta_Ct_summary = ddply(mdata_Delta_Ct,.(Gene,Treatment),  summarise,
                          mean_DeltaCt = mean(DeltaCt),
                          CI_Lower = -se(DeltaCt)+mean(DeltaCt),
                          CI_Upper = se(DeltaCt)+mean(DeltaCt))
colours=rep(c("red","green","blue"),11)

#Plot Delta Ct
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


ylims <- mdata_Delta_Ct %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(floor(min(DeltaCt)),ceiling(max(DeltaCt)))
colnames(ylims)[2:3] <- c("ymin","ymax")

p = ggplot(data=mdata_Delta_Ct,aes(x=Treatment,y=DeltaCt,fill=Treatment,color=Treatment))+
  #next line stopes error bars overlapping border
  geom_errorbar(data=mdata_Delta_Ct_summary,aes(x=Treatment, y=mean_DeltaCt,ymin=CI_Lower-0.5, ymax=CI_Upper+0.5), width=.4,colour="white")+
  geom_errorbar(data=mdata_Delta_Ct_summary,aes(x=Treatment, y=mean_DeltaCt,ymin=CI_Lower, ymax=CI_Upper), width=.4,color="black",)+
  geom_jitter(width=0.25,height=0,size=1, alpha=0.6) +
  geom_point(data=mdata_Delta_Ct_summary,stat='identity', aes(x=Treatment, y=mean_DeltaCt),size=1,shape=2,color="black",fill="black")+
  
  scale_color_manual(values = c('grey40', 'pink','cyan'))+
  scale_y_continuous(expand = c(0, 0))+
  theme_bw()+
  xlab(NULL)+
  ylab(expression(log[2]*"(Relative Expression)"))+
  theme(panel.grid = element_blank(),
        axis.text.x =element_blank(),# element_text(angle = 90,vjust = 0, color='grey30',size=12),
        axis.title.y = element_text(colour = "grey30", face='bold',size=12),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(size = 10,face = "italic"),
        legend.position = "none",
        #       panel.border = element_blank(),strip.text = element_text(face = "italic")
        axis.line = element_line(colour = "grey30"))+
  facet_wrap_custom(~Gene,scales="free_y",nrow=3,scale_overrides = list(
    scale_override(1, scale_y_continuous(limits = c(as.integer(ylims[1,2]), as.integer(ylims[1,3])))),
    scale_override(2, scale_y_continuous(limits = c(as.integer(ylims[2,2]), as.integer(ylims[2,3])))),
    scale_override(3, scale_y_continuous(limits = c(as.integer(ylims[3,2]), as.integer(ylims[3,3])))),
    scale_override(4, scale_y_continuous(limits = c(as.integer(ylims[4,2]), as.integer(ylims[4,3])))),
    scale_override(5, scale_y_continuous(limits = c(as.integer(ylims[5,2]), as.integer(ylims[5,3])))),
    scale_override(6, scale_y_continuous(limits = c(as.integer(ylims[6,2]), as.integer(ylims[6,3])))),
    scale_override(7, scale_y_continuous(limits = c(as.integer(ylims[7,2]), as.integer(ylims[7,3])))),
    scale_override(8, scale_y_continuous(limits = c(as.integer(ylims[8,2]), as.integer(ylims[8,3])))),
    scale_override(9, scale_y_continuous(limits = c(as.integer(ylims[9,2]), as.integer(ylims[9,3])))),
    scale_override(10, scale_y_continuous(limits = c(as.integer(ylims[10,2]), as.integer(ylims[10,3])))),
    scale_override(11, scale_y_continuous(limits = c(as.integer(ylims[11,2]), as.integer(ylims[11,3]))))
  ))
p

pdf(file="qPCR.pdf",height=4.5,width=5)
p
dev.off()


#legend
pdf(file="legend.pdf",height=3,width=6.5)

plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
text(0,7,"Unchallenged",pos=4,col="grey40",font=2)
text(0,5,"Oil",pos=4,col="cyan",font=2)
text(0,3,"Wasp Homogenate",pos=4,col="pink",font=2)

dev.off()


#Plot DeltaDelta Ct - just averages
DeltaDelta_Ct_summary = spread(mdata_Delta_Ct_summary[,1:3], Treatment, mean_DeltaCt)
DeltaDelta_Ct_summary$Oil_FC = 2^-(DeltaDelta_Ct_summary$Control-DeltaDelta_Ct_summary$Oil)
DeltaDelta_Ct_summary$Wasp_Extract_FC = 2^-(DeltaDelta_Ct_summary$Control-DeltaDelta_Ct_summary$Wasp_Extract)
mDeltaDelta_Ct_summary = melt(DeltaDelta_Ct_summary[,c(1,5,6)], id.vars = c('Gene'))
colnames(mDeltaDelta_Ct_summary) = c('Gene', 'Treatment', 'FC')

#Change Gene IDs. Serine proteases will  match Cao 2018 "Building a platform for predicting functions of serine protease-related proteins in Drosophila melanogaster and other insects"
Gene_IDs = read.table("genes_tested_by_qPCR_IDs.csv", sep = ",", dec = ".", header = TRUE)
Gene_IDs
colnames(DeltaDelta_Ct_summary)[1]='Symbol'
DeltaDelta_Ct_summary = merge(DeltaDelta_Ct_summary, Gene_IDs, by='Symbol')

p2 = ggplot(DeltaDelta_Ct_summary, aes(x=log2(Wasp_Extract_FC), y=log2(Oil_FC)))+
  geom_point()+
  geom_text(aes(y=log2(Oil_FC)+0.2,label=Gene))+
  scale_y_continuous(limits = c(-1, 8))+
  scale_x_continuous(limits = c(-1, 8))+
  geom_abline(intercept = 0, slope = 1, linetype='dashed')+
  geom_vline(xintercept = 0, linetype='dashed')+
  geom_hline(yintercept = 0, linetype='dashed')+
  ylab(expression(paste('Oil vs. Control log'[2], '(FC)')))+
  xlab(expression(paste('Wasp homogenate vs. Control log'[2], '(FC)')))+
  theme_bw()+
  theme(axis.line=element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(angle = 90,vjust = 0, color='grey30',size=15),
        axis.text.x = element_text(angle = 0,vjust = 0, color='grey30',size=15),
        axis.title.y = element_text(colour = "grey30", face='bold',size=15),
        axis.title.x = element_text(colour = "grey30", face='bold',size=15, margin=margin(t=20)))
p2



####Stats
#Create dataset to add p values for all comparisions: Control Vs. Oil and Control Vs. Wasp homogenate
pvalues = data.frame(Symbol=rep(NA,11), pvalue_Oil=rep(NA,11), pvalue_Wasp=rep(NA,11))
pvalues

for(i in 1:11){
  m = lm(data=mdata_Delta_Ct, DeltaCt~Treatment, subset=Gene==levels(mdata_Delta_Ct$Gene)[i])
  pvalues$Symbol[i]=levels(mdata_Delta_Ct$Gene)[i]
  pvalues$pvalue_Oil[i]=coef(summary(m))[,4][2]
  pvalues$pvalue_Wasp[i]=coef(summary(m))[,4][3]
}

pvalues
