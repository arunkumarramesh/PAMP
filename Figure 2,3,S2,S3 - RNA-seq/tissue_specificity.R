setwd("/Users/ramesha/Downloads/")

## flyatlas2 data, courtesy of David Leader and Darren Obbard , values are FPKM
FlyAtlas2_Alltissues_Allgenes <- read.csv(file="FlyAtlas2_Alltissues_Allgenes.csv")

## tissue specific index
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-x/max(x))/(length(x)-1)
}


FlyAtlas2_larvae <- FlyAtlas2_Alltissues_Allgenes[grep("^larval",colnames(FlyAtlas2_Alltissues_Allgenes))] ### get larval tissues
FlyAtlas2_larvae$tau <- apply(FlyAtlas2_larvae,1,tau) # calculate tau
rownames(FlyAtlas2_larvae) <- FlyAtlas2_Alltissues_Allgenes$FlyBaseID
FlyAtlas2_larvae <- FlyAtlas2_larvae[FlyAtlas2_larvae$tau > 0.8,] ## only keep tissue specific genes, tau > 0.8
FlyAtlas2_larvae <- FlyAtlas2_larvae[!is.na(FlyAtlas2_larvae$tau),] ## remove those without tau
## identify which tissue has maximal expression
max_tissue <- ""
for (i in 1:nrow(FlyAtlas2_larvae)){
  max_tissue[i] <- colnames(FlyAtlas2_larvae)[which(FlyAtlas2_larvae[i,] == max(FlyAtlas2_larvae[i,]))]
}
FlyAtlas2_larvae <- FlyAtlas2_larvae[max_tissue %in% "larval_Salivary_gland",] ## only keep tissue specific keep  where Salivary_gland has maximal expression
FlyAtlas2_larvae <- FlyAtlas2_larvae[FlyAtlas2_larvae$larval_Salivary_gland > 1,] ## only keep genes with FPKM above 1

#pheatmap(FlyAtlas2_larvae[1:(ncol(FlyAtlas2_larvae)-1)]/(rowMeans(FlyAtlas2_larvae[1:(ncol(FlyAtlas2_larvae)-1)])))


FlyAtlas2_adult_male <- FlyAtlas2_Alltissues_Allgenes[grep("^male",colnames(FlyAtlas2_Alltissues_Allgenes))] ### get male adult tissues
FlyAtlas2_adult_male$tau <- apply(FlyAtlas2_adult_male,1,tau) # calculate tau
rownames(FlyAtlas2_adult_male) <- FlyAtlas2_Alltissues_Allgenes$FlyBaseID
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[FlyAtlas2_adult_male$tau > 0.8,] ## only keep tissue specific genes, tau > 0.8
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[!is.na(FlyAtlas2_adult_male$tau),] ## remove those without tau
## identify which tissue has maximal expression
max_tissue <- ""
for (i in 1:nrow(FlyAtlas2_adult_male)){
  max_tissue[i] <- colnames(FlyAtlas2_adult_male)[which(FlyAtlas2_adult_male[i,] == max(FlyAtlas2_adult_male[i,]))]
}
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[max_tissue %in% "male_Testis",] ## only keep tissue specific keep  where male testis has maximal expression
FlyAtlas2_adult_male <- FlyAtlas2_adult_male[FlyAtlas2_adult_male$male_Testis > 1,] ## only keep genes with FPKM above 1

#pheatmap(FlyAtlas2_adult_male[1:(ncol(FlyAtlas2_adult_male)-1)]/(rowMeans(FlyAtlas2_adult_male[1:(ncol(FlyAtlas2_adult_male)-1)])))
