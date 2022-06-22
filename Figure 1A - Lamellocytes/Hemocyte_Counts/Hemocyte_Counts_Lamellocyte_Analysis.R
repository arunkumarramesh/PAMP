#load required packages
library(ggplot2)
library(plyr)

#load dataset 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data = read.table('Hemocyte_Counts.csv', sep=',',dec='.', head=T)
data$Treatment=as.factor(data$Treatment)

#Define standard error function
se <- function(x) sd(x)/sqrt(length(x))

#Create data summary table with plasmatocytes and lamellocyte mean and standard error
data_summary = ddply(data,.(Treatment),summarise,
                     mean_Plasmatocytes = mean(Plasmatocytes, na.rm =TRUE),
                     se_Plasmatocytes = se(Plasmatocytes),
                     mean_Lamellocytes = mean(Lamellocytes,na.rm=TRUE),
                     se_Lamellocytes = se(Lamellocytes))

### Plot data
p=ggplot(data_summary, aes(x=Treatment, y=mean_Lamellocytes))+
  geom_bar(stat='identity',fill="grey50",width=0.8)+
  geom_jitter(data=data, aes(x=Treatment, y=Lamellocytes),color='grey10', width = 0.1)+
  geom_errorbar(data=data_summary, aes(ymin=mean_Lamellocytes-se_Lamellocytes, ymax=mean_Lamellocytes+se_Lamellocytes), width=.4)+ 
  geom_text(aes(x=Treatment, y=28, label=c('a','b','c')), , size=4,colour = "grey30")+
  scale_y_continuous(limits = c(0, 30),expand = c(0, 0))+
  scale_x_discrete(labels=c("Control" = 'Unchallanged',
                            "Oil" = "Oil",
                            "Oil + Wasp Extract" = "Oil +\n\u2642 Wasp homogenate"))+
  theme_bw()+
  xlab("\nImmune Challenge")+
  ylab(expression(paste('Lamellocytes (x100/ ',mu*L,')')))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,  hjust = 1, vjust = 1, color='grey30', size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(colour = "grey30",  size=14),
        axis.title.x = element_text(colour = "grey30",  size=12),
        panel.border = element_blank(),
        axis.line = element_line(colour = "grey30"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
p


svg(file="Lamellocytes.svg",height=5,width=3)
p
dev.off()


###Stats
m1=aov(data=data, Lamellocytes~Treatment)
summary(m1)
TukeyHSD(x=m1, 'Treatment', conf.level=0.95)

###Stats
library(car)
library(multcomp)
m1=lm(data=data, Lamellocytes~Treatment)
summary(m1)
Anova(m1)
summary(glht(model = m1, linfct = mcp(Treatment = "Tukey")),test = adjusted('bonferroni'))
