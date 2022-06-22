#load required packages
library(ggplot2)
library(plyr)
library(Hmisc)

#load dataset 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data = read.table('Oil_Melanization_Phenotype_Characterization.csv', sep=',',dec='.', head=T)

#Reorder Treatment factor levels
data$Treatment = factor(data$Treatment, levels = c('Oil', 'Oil + Wasp Infection',
                                                   'Oil + Female Wasp Extract', 'Oil + Male Wasp Extract',
                                                   'Oil + Fly Extract'))

#Create data summary table
data_summary = ddply(data,.(Treatment),summarise, 
                     Melanized = sum(Melanized), 
                     Total = sum(Total))

data_summary_2 = ddply(data_summary,.(Treatment),summarise, 
                       Proportion_Melanization = binconf(x=Melanized, n=Total, alpha=0.05, method='wilson')[1],
                       Lower_Confidence_Interval_Value = binconf(x=Melanized, n=Total, alpha=0.05, method='wilson')[2],
                       Upper_Confidence_Interval_Value = binconf(x=Melanized, n=Total, alpha=0.05, method='wilson')[3],
                       Total = sum(Total))

### Plot data
p=ggplot(data_summary_2, aes(x=Treatment, y=Proportion_Melanization))+
  geom_bar(stat='identity',width=0.8,fill="grey50")+
  geom_errorbar(data=data_summary_2, aes(ymin=Lower_Confidence_Interval_Value, ymax=Upper_Confidence_Interval_Value), width=.4)+ 
  geom_text(aes(y=Upper_Confidence_Interval_Value+0.025, label=paste('n=',Total, sep='')), size=4,colour = "grey30")+
  geom_text(aes(x=Treatment, y=0.99, label=c('a','b','b','b','a')), size=4,colour = "grey30")+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_discrete(labels=c("Oil + Wasp Infection" = 'Oil +\nWasp infection',
                            "Oil + Female Wasp Extract" = "Oil +\n\u2640 Wasp homogenate",
                            "Oil + Male Wasp Extract" = "Oil +\n\u2642 Wasp homogenate",
                            "Oil + Fly Extract"="Oil +\nFly homogenate"))+
  ylab(expression('Proportion Oil Droplets Melanized'))+
  xlab('\nImmune Challenge')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,  hjust = 1, vjust = 1, color='grey30', size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(colour = "grey30",  size=14),
        axis.title.x = element_text(colour = "grey30",  size=12),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
p

#use this to embed male and female symbols in pdf. Think only works on mac (?)
#quartz.save(type = 'pdf', file="Melanisation.pdf",height=5,width=5)
#p
#dev.off()


svg(file="Oil_Melanization_Characterization.svg",height=5,width=5)
p
dev.off()

###Stats
library(car)
library(multcomp)
m1= glm(data=data, cbind(data$Melanized,data$Non_Melanized)~Treatment,family="quasibinomial")
summary(m1)
Anova(m1)
summary(glht(model = m1, linfct = mcp(Treatment = "Tukey")),test = adjusted('bonferroni'))

