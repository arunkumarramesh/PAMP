library(ape)
library(MCMCglmm)

setwd("/scratch/Arun/Projects/Oil_injection")

Tree <- read.tree("ins.phy")
Data <- read.csv("Oil_Injection_Insect_Species_Phenotypes_V3.csv", row.names = 1) 
Data_subset2 <- as.data.frame(subset(Data, select = c("Melanized","Non_Melanized", "Parasitoid","Divergence_Time_from_Drosophila_melanogaster", "Specimen_part")))
Data_taxon2 <- data.frame(Data_subset2, animal=row.names(Data_subset2))
prior_MCMCglmm_disc_cont <- list (R=list(V=1,nu=0.002), G=list(G1=list(V=1, nu=0.002)))

##model 3 used
model <-MCMCglmm(cbind(Melanized, Non_Melanized)~ 1,random=~animal, family="multinomial2", data=Data_taxon2,prior=prior_MCMCglmm_disc_cont,pedigree=Tree, pl = TRUE, nitt=100000000, burnin=10000000, thin=100000, verbose=TRUE)
model2 <-MCMCglmm(cbind(Melanized, Non_Melanized)~ Divergence_Time_from_Drosophila_melanogaster,random=~animal, family="multinomial2", data=Data_taxon2,prior=prior_MCMCglmm_disc_cont,pedigree=Tree, pl = TRUE, nitt=100000000, burnin=10000000, thin=100000, verbose=TRUE)
model3 <-MCMCglmm(cbind(Melanized, Non_Melanized)~ Divergence_Time_from_Drosophila_melanogaster + Parasitoid,random=~animal, family="multinomial2", data=Data_taxon2,prior=prior_MCMCglmm_disc_cont,pedigree=Tree, pl = TRUE, nitt=100000000, burnin=10000000, thin=100000, verbose=TRUE)
model4 <-MCMCglmm(cbind(Melanized, Non_Melanized)~ Divergence_Time_from_Drosophila_melanogaster + Parasitoid + Specimen_part,random=~animal, family="multinomial2", data=Data_taxon2,prior=prior_MCMCglmm_disc_cont,pedigree=Tree, pl = TRUE, nitt=100000000, burnin=10000000, thin=100000, verbose=TRUE)
model5 <-MCMCglmm(cbind(Melanized, Non_Melanized)~ Divergence_Time_from_Drosophila_melanogaster +  Specimen_part,random=~animal, family="multinomial2", data=Data_taxon2,prior=prior_MCMCglmm_disc_cont,pedigree=Tree, pl = TRUE, nitt=100000000, burnin=10000000, thin=100000, verbose=TRUE)
