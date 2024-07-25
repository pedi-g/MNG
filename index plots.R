#Pedram Ghahramani
#pedram.ghrm@gmail.com
#visualization of indices
#6-27-2024
#libraries----------
library(tidyverse)
library(ggthemes)
library(ggfortify)
library(ggpubr)
library(psych)

theme_set(theme_bw())
#input-----------
inds <- read.csv("output/indices.csv", row.names = 1)
# inds<- inds[inds[1]>1,]
indsL <- read.csv("output/indices(log).csv", row.names = 1)
inds<- indsL
# indsL<- indsL[indsL[1]>1,]
inds1<- inds%>%
  select(-Status,-habitat, -area,- station, -season)
indsL1<- indsL%>%
  select(-Status,-habitat, -area,- station, -season)
#defining factor variables --------
source('R-scripts/factors.R')

inds<- cbind(inds1,factors)%>%
  filter(S !=1)
inds1<- inds1%>%
  filter(S !=1)
factors<- factors[rownames(factors) %in% rownames(inds1),]
####################### mean and standard error of all indices ##############################

meanSAH<-list()
meanAH<-list()
meanSA<- list()
meanST<- list()
for (i in seq(ncol(inds1))){
  dum<-data.frame(index = inds1[,i],season = factors$season,
                  area = factors$area, habitat = factors$habitat)
  
  ST<- data.frame(index = inds1[,i],station = factors$station,
                 habitat = factors$habitat) %>%
    group_by(station, .add = T)%>%
    summarise(mean = mean(index),SD = sd(index))
meanST[[length(meanST)+1]]<- ST
names(meanST)[length(meanST)]<-names(inds1[i])

  SAH<- dum %>%
    group_by(season, area, habitat)%>%
    summarise(mean = mean(index),SD = sd(index))
meanSAH[[length(meanSAH)+1]]<- SAH
names(meanSAH)[length(meanSAH)]<-names(inds1[i])

SA<- dum %>%
  group_by(season, area)%>%
  summarise(mean = mean(index),SD = sd(index))
meanSA[[length(meanSA)+1]]<- SA
names(meanSA)[length(meanSA)]<-names(inds1[i])

AH<- dum %>%
  group_by(area, habitat)%>%
  summarise(mean = mean(index),SD = sd(index))
meanAH[[length(meanAH)+1]]<- AH
names(meanAH)[length(meanAH)]<-names(inds1[i])

}

################ bar plots and box plots for all 3 factors ##################

plotSAH<- function(index= "H", f1 = "season", f2 = "area", f3 = "habitat"){
   colind <- which(colnames(inds) == index)
  ############### barplot #######################.
  ind<-as.data.frame(meanSAH[[as.character(index)]])
  if(f1 == "season" & f2 == "area"){
    factor<- interaction(ind$area,ind$season)
    fill<- ind$habitat
    name<- "season.area"
  }
  if(f1 == "season" & f2 == "habitat"){
    factor<- interaction(ind$season,ind$habitat)
    fill<- ind$area
    name<- "season.habitat"
  }
  if(f1 == "area" & f2 == "habitat"){
    factor<- interaction(ind$area,ind$habitat)
    fill<- ind$season
    name<- "area.habitat"
  }
  ind<- data.frame(fill,factor, ind[4:5])
  
  bar<- ggplot(ind, aes(factor, mean,
                        fill = fill))+
    geom_col(position = "dodge",col = "black")+ #by dodge the fill factor will be side to side instead of stacked
    labs(x = name, y = as.character(index))+
    geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = c("#E7298A", "#66A61E", "#E6AB02"))
    # scale_fill_brewer(palette = "Dark2")
    #scale_color_manual(values = c("darkblue", "#E6AB02"))
    #geom_line(aes(group = area))
  
  ################ Boxplot #########################.
  if(f1 == "season" & f2 == "area"){
    factor<- interaction(inds$season,inds$area)
    fill<- inds$habitat
  }
  if(f1 == "season" & f2 == "habitat"){
    factor<- interaction(inds$season,inds$habitat)
    fill<- inds$area
  }
  if(f1 == "area" & f2 == "habitat"){
    factor<- interaction(inds$area,inds$habitat)
    fill<- inds$season
  }
  
  data<- data.frame(fill,factor, value = inds[,colind])
  box<- ggplot(data, aes(factor, value, fill = fill))+
    geom_boxplot(position = "dodge")+
    labs(x = name, y = as.character(index))+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"))+
    # scale_fill_brewer(palette = "Dark2", direction = -1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(bar,box))
  
}

AllSAH<-list()
for (i in names(inds1)){
  AllSAH[[length(AllSAH)+1]]<- plotSAH(i)
  names(AllSAH)[length(AllSAH)]<-names(inds1[i])
}

AllSAH[["S"]][[1]]

SAHplots<- ggarrange(AllSAH[["S"]][[1]], AllSAH[["J"]][[1]], AllSAH[["H"]][[1]], 
          AllSAH[["FRic"]][[1]], AllSAH[["FEve"]][[1]], AllSAH[["Rao"]][[1]],
          common.legend = TRUE)
 ggsave("figs/indices_3factors_6inds.png", width = 12, height = 8)

###################### EcoQS#########################################

EcoQplot<- ggarrange(AllSAH[["AMBI"]][[2]], AllSAH[["BENTIX"]][[2]], AllSAH[["M.AMBI"]][[2]], 
                     common.legend = TRUE)
################ plots for season and area only#####################
plotSA<- function(index= "H", f1 = "season", f2 = "area"){
  colind <- which(colnames(inds) == index)
  ############### barplot #######################.
  ind<-as.data.frame(meanSA[[as.character(index)]])
 
  factor<- ind$area
  fill <- ind$season
  ind<- data.frame(fill,factor, ind[3:4])
  
  bar<- ggplot(ind, aes(factor, mean,
                        fill = fill))+
    geom_col(position = "dodge",col = "black")+ #by dodge the fill factor will be side to side instead of stacked
    labs(x = "", y = as.character(index))+
    geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = c("#66A61E", "#E6AB02"))

  ################ Boxplot #########################.
  factor<- inds$area
  fill<- inds$season
  
  data<- data.frame(fill,factor, value = inds[,colind])
  box<- ggplot(data, aes(factor, value, fill = fill))+
    geom_boxplot(position = "dodge")+
    labs(x = "", y = as.character(index))+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"))+
    # scale_fill_brewer(palette = "Dark2", direction = -1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(bar,box))
}


AllSA<-list()
for (i in names(inds1)){
  AllSA[[length(AllSA)+1]]<- plotSA(i)
  names(AllSA)[length(AllSA)]<-names(inds1[i])
}

AllSA[["S"]][[1]]
ggarrange(AllSA[["S"]][[1]], AllSA[["J"]][[1]], AllSA[["H"]][[1]], 
                     AllSA[["FRic"]][[1]], AllSA[["FEve"]][[1]], AllSA[["Rao"]][[1]],
                     common.legend = TRUE)
  ggsave("figs/indices_area.season_6inds.png", width = 7, height = 5)


################ plots for area and habitat only ########################
plotAH<- function(index= "H", f1 = "area", f2 = "habitat"){
  colind <- which(colnames(inds) == index)
  ############### barplot #######################.
  ind<-as.data.frame(meanAH[[as.character(index)]])
  
  factor<- ind$area
  fill <- ind$habitat
  ind<- data.frame(fill,factor, ind[3:4])
  
  bar<- ggplot(ind, aes(factor, mean,
                        fill = fill))+
    geom_col(position = "dodge",col = "black")+ #by dodge the fill factor will be side to side instead of stacked
    labs(x = "", y = as.character(index))+
    geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = c("#E7298A", "#66A61E", "#E6AB02"))
  
  ################ Boxplot #########################.
  factor<- inds$area
  fill<- inds$habitat
  
  data<- data.frame(fill,factor, value = inds[,colind])
  box<- ggplot(data, aes(factor, value, fill = fill))+
    geom_boxplot(position = "dodge")+
    labs(x = "", y = as.character(index))+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"))+
    # scale_fill_brewer(palette = "Dark2", direction = -1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(bar,box))
}


AllAH<-list()
for (i in names(inds1)){
  AllAH[[length(AllAH)+1]]<- plotAH(i)
  names(AllAH)[length(AllAH)]<-names(inds1[i])
}

AllAH[["S"]][[1]]
ggarrange(AllAH[["S"]][[1]], AllAH[["J"]][[1]], AllAH[["H"]][[1]], 
          AllAH[["FRic"]][[1]], AllAH[["FEve"]][[1]], AllAH[["Rao"]][[1]],
          common.legend = TRUE)
 ggsave("figs/indices_area.habitat_6inds.png", width = 8, height = 5)

################ plots divided by stations ########################

plotST<- function(index= "H", f1 = "transect"){
  colind <- which(colnames(inds) == index)
  ############### barplot #######################.
  ind<-as.data.frame(meanST[[as.character(index)]])
  
  factor<- ind$station
  ind<- data.frame(factor, ind[2:3])
  
  bar<- ggplot(ind, aes(factor, mean, fill = factor))+
    geom_col(position = "dodge",col = "black")+ #by dodge the fill factor will be side to side instead of stacked
    labs(x = "", y = as.character(index))+
    geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    #scale_fill_manual(values = c("#66A61E", "#E6AB02"))
  
  ################ Boxplot #########################.
  factor<- inds$station

  data<- data.frame(factor, value = inds[,colind], fill = factor)
  box<- ggplot(data, aes(factor, value))+
    geom_boxplot(position = "dodge")+
    labs(x = "", y = as.character(index))+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"))+
    # scale_fill_brewer(palette = "Dark2", direction = -1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(list(bar,box))
}

AllST<-list()
for (i in names(inds1)){
  AllST[[length(AllST)+1]]<- plotST(i)
  names(AllST)[length(AllST)]<-names(inds1[i])
}

AllST[["S"]][[1]]
ggarrange(AllST[["S"]][[1]], AllST[["J"]][[1]], AllST[["H"]][[1]], 
          AllST[["FRic"]][[1]], AllST[["FEve"]][[1]], AllST[["Rao"]][[1]],
          common.legend = TRUE)
ggsave("figs/indices_area.habitat_6inds.png", width = 8, height = 5)

# scale_fill_brewer(palette = "Dark2")

#boxplot(S ~ habitat:season:area, data = inds)

############################### correlation#############################

pairs.panels(inds1)
pairs.panels(indsL1)
png("figs/corr.inds.pairs1.png", width = 800, height = 1000)
pairs.panels(indsL1, gap=0,pch=20,method="spearman", 
             ellipses = F, lm= T, density = F, hist.col ="#236573", 
             jiggle = F,rug=F,smoother=F,stars=T,ci=T,alpha=.05,)
dev.off()




