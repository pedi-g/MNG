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
inds <- read.csv("output/indices.csv", row.names = 1)%>%
  filter(S>1) %>%
  select(-Dm, -D ,
         -AMBI, - BENTIX, -M.AMBI, -Status)
indsL <- read.csv("output/indices(log).csv", row.names = 1)%>%
  select(-Dm, -D, 
         # -AMBI, - BENTIX, -M.AMBI,
         -Status)
ab<-rowSums(log(1+(read.csv('raw data/ab.csv', row.names = 1)%>%
                               select(-Habitat, -Site,- station, -Season))))
mass<- rowSums(log(1+(read.csv('raw data/mass.csv', row.names = 1))))
indsL<- cbind(indsL,ab, mass)
source('R-scripts/factors.R') # factors
inds1<- inds%>%
  select(-Habitat, -Site,- station, -Season)
indsL1<- indsL%>%
  select(-Habitat, -Site,- station, -Season)
inds<- cbind(inds1,factors)
indsL<- cbind(indsL1,factors)


### deleting Bahu


# ### threshold 2.5
# indsL[row.names(indsL) == "C513_s",]$S <- NA
# indsL[row.names(indsL) == "C312_w",]$J <- NA
# indsL[row.names(indsL) == "C423_w",]$J <- NA
# indsL[row.names(indsL) == "C723_w",]$FEve <- NA
# 
# # threshhold 2
# indsL[row.names(indsL) == "C512_s",]$S <- NA
# indsL[row.names(indsL) == "C422_w",]$H <- NA
# indsL[row.names(indsL) == "C513_s",]$H <- NA
# indsL[row.names(indsL) == "C311_s",]$J <- NA
# indsL[row.names(indsL) == "C821_w",]$FRic <- NA
# indsL[row.names(indsL) == "C512_s",]$FRic <- NA
# indsL[row.names(indsL) == "C222_w",]$FEve <- NA
# indsL[row.names(indsL) == "C412_w",]$Rao <- NA
# indsL[row.names(indsL) == "C421_w",]$Rao <- NA
# indsL[row.names(indsL) == "C423_w",]$Rao <- NA
# indsL[row.names(indsL) == "C522_w",]$Rao <- NA

### personal choice deletion
# indsL[row.names(indsL) == "C513_s",]$H <- NA


# eliminating bch Habitat:
inds<- inds[inds$Habitat!="bch",]
indsL<- indsL[indsL$Habitat!="bch",]

inds<- inds[inds$Site!="Bah",]
indsL<- indsL[indsL$Site!="Bah",]
# inds<- inds%>%
#   filter(transect != "1", transect != "2", transect!= "7", transect!= "8")
# indsL<- indsL%>%
#   filter(transect != "1", transect != "2", transect!= "7", transect!= "8")
inds1<- inds1[row.names(inds1) %in% row.names(inds),]
indsL1<- indsL1[row.names(indsL1) %in% row.names(inds),]
 
# eliminating stations with richness of one (all the species)
# inds<- inds[inds$S2>1,]
# indsL<- indsL[indsL$S2>1,]
inds1<- inds1[row.names(inds1) %in% row.names(inds),]
indsL1<- indsL1[row.names(indsL1) %in% row.names(inds),]

#eliminating stations with richness of one (trait species)
# inds<- inds[inds$S2>2,]
# indsL<- indsL[indsL$S2>2,]
# inds1<- inds1[row.names(inds1) %in% row.names(inds),]
# indsL1<- indsL1[row.names(indsL1) %in% row.names(inds),]

##### turning NAs to 0s
# inds[is.na(inds)]<- 0
# indsL[is.na(inds)]<- 0

########################## histograms #########################################
########## not transformed 
png('figs/ind_hist_s2.png',width=900, height=700) 
par(mfrow = c(2,4))
for(i in c(1,3:8)){
  hist(inds1[,i], ylab = '',
       main = names(inds1)[i])
}
dev.off()

########## log transformed 
png('figs/ind_hist_log_s2.png',width=900, height=700)
par(mfrow = c(2,4))
for(i in c(1,3:8)){
  hist(indsL1[,i], ylab = '',
       main = names(inds1)[i])
}
dev.off()

################################# bar and box #################################

################## a function for bar and box plots #
plotSAH<- function(index= "H", f1 = "Season", f2 = "Site", f3 = "Habitat"){
  colind <- which(colnames(indsL) == index)
  
  ############### barplot #####.
  ind<-as.data.frame(meanSAH[[as.character(index)]])
  
  if(f1 == "Season" & f2 == "Site"){
    factor<- interaction(ind$Site,ind$Season)
    fill<- ind$Habitat
    name<- "Site.Season"
  }
  if(f1 == "Season" & f2 == "Habitat"){
    factor<- interaction(ind$Season,ind$Habitat)
    fill<- ind$Site
    name<- "Season.Habitat"
  }
  if(f1 == "Site" & f2 == "Habitat"){
    factor<- interaction(ind$Site,ind$Habitat)
    fill<- ind$Season
    name<- "Site.Habitat"
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
  #geom_line(aes(group = Site))
  
  ################ Boxplot #####.
  if(f1 == "Season" & f2 == "Site"){
    factor<- interaction(inds$Site,inds$Season)
    fill<- inds$Habitat
  }
  if(f1 == "Season" & f2 == "Habitat"){
    factor<- interaction(inds$Season,inds$Habitat)
    fill<- inds$Site
  }
  if(f1 == "Site" & f2 == "Habitat"){
    factor<- interaction(inds$Site,inds$Habitat)
    fill<- inds$Season
  }
  
  data<- data.frame(fill,factor, value = df[,colind])
  data <- data[!is.na(data$value) & is.finite(data$value), ]
  box<- ggplot(data, aes(factor, value, fill = fill))+
    geom_boxplot(position = "dodge")+
    labs(x = name, y = as.character(index), fill = "")+
    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"),
                      labels = c("Mudflat", "Vegetated"))+
    # scale_y_continuous(
    #   breaks = round(seq(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE), length.out = 4), 2),
    #   expand = expansion(mult = c(0.1, 0.2))
    # ) +
    # scale_fill_brewer(palette = "Dark2", direction = -1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 9))
  return(list(bar, box, data))
  
}

################ mean and standard error of all indices #######################

df<- indsL #or indsL

meanSAH<-list()
for (i in seq(ncol(indsL1))){
  dum<-data.frame(index = df[,i],Season = inds$Season,
                  Site = inds$Site, Habitat = inds$Habitat)

  SAH<- dum %>%
    group_by(Season, Site, Habitat)%>%
    summarise(mean = mean(index, na.rm =T),SD = sd(index, na.rm =T))
meanSAH[[length(meanSAH)+1]]<- SAH
names(meanSAH)[length(meanSAH)]<-names(indsL1[i])

}

################ bar plots and box plots for all 3 factors ##################

AllSAH<-list()
for (i in names(indsL1)){
  AllSAH[[length(AllSAH)+1]]<- plotSAH(i)
  names(AllSAH)[length(AllSAH)]<-names(indsL1[i])
}

######## editing x axis lable 
for (i in 1:10) {
  AllSAH[[i]][[2]] <- AllSAH[[i]][[2]] + theme(axis.title.x = element_blank())
  AllSAH[[i]][[1]] <- AllSAH[[i]][[1]] + theme(axis.title.x = element_blank())
}

bar<- list()
for( i in c(1,3:10)){
  bar[[length(bar)+1]]<-  AllSAH[[i]][[1]]
}

box<- list()
for( i in c(1,3:10)){
 box[[length(box)+1]]<-  AllSAH[[i]][[2]]
}


 p1<- ggarrange(bar[[1]],bar[[3]],bar[[2]],bar[[4]],bar[[5]],bar[[6]],ggplot() + theme_void(),
          bar[[7]], common.legend = T,
          labels = c("a",'b','c','d','e','f','g'),
          font.label = list(size = 14, color = "black", face = "plain", family = NULL),hjust = -1.4)

############### Mehtod two for box plot 2 #############################

 box[[1]] <- box[[1]] + annotate("text", x = 4.3, y = 1.5, label = "a", size = 5, color = "gray30")
 box[[3]] <- box[[3]] + annotate("text", x = 4.3, y = 0.805, label = "b", size = 5, color = "gray30")
 box[[2]] <- box[[2]] + annotate("text", x = 4.3, y = 0.8, label = "c", size = 5, color = "gray30")
 box[[4]] <- box[[4]] + annotate("text", x = 4.3, y = 0.05, label = "d", size = 5, color = "gray30")
 box[[5]] <- box[[5]] + annotate("text", x = 4.3, y = 0.5, label = "e", size = 5, color = "gray30")
 box[[6]] <- box[[6]] + annotate("text", x = 4.3, y = 0.17, label = "f", size = 5, color = "gray30")

 # Combine plots with consistent labels
 p2 <- ggarrange(box[[1]], box[[3]], box[[2]], box[[4]], box[[5]], box[[6]],
                 common.legend = TRUE)
 p2

 ############ method one for box plot 1 ###############
 # top_margin_theme <- theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 5))
 # Update the top three plots with the increased top margin
 # box[[1]] <- box[[1]] + top_margin_theme
 # box[[3]] <- box[[3]] + top_margin_theme
 # box[[2]] <- box[[2]] + top_margin_theme
 # box[[4]] <- box[[4]] + top_margin_theme
 # box[[5]] <- box[[5]] + top_margin_theme
 # box[[6]] <- box[[6]] + top_margin_theme
 # 
 # p2<-
 #  ggarrange(box[[1]],box[[3]],box[[2]],box[[4]],box[[5]],box[[6]] , common.legend = T,
 #               labels = c('a','b','c','d','e',' f'),
 #               font.label = list(size = 13, color = "gray40", face = "bold", family = NULL),
 #               hjust = -3.2, vjust = 1.1)
 # p2
 # 
 # 
 # # SAHplots<- ggarrange(AllSAH[["S"]][[1]], AllSAH[["J"]][[1]], AllSAH[["H"]][[1]], 
 # #           AllSAH[["FRic"]][[1]], AllSAH[["FEve"]][[1]], AllSAH[["Rao"]][[1]],
 # #           common.legend = TRUE)
 
 ############# box inds manually ########################
 
#  data<- AllSAH$S[[3]]
#  index<- 'S'
#  name<- "Site.Season"
#  (box[[1]] <- ggplot(data, aes(factor, value, fill = fill))+
#    geom_boxplot(position = "dodge")+
#    labs(x = name, y = as.character(index), fill = "")+
#    scale_fill_manual(values = c("#7570b3", "#1b9e77", "#d95f02"),
#                      labels = c("Mudflat", "Vegetated"))+
#    scale_y_continuous(
#    # breaks = round(seq(min(data$value, na.rm = TRUE), max(data$value, na.rm = TRUE), length.out = 5), 2),
#    expand = expansion(mult = c(0.1, 0.1))
#    ) +
#    # scale_fill_brewer(palette = "Dark2", direction = -1)+
#    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#          axis.text.y = element_text(size = 9),
#          axis.title.x = element_blank()))
#  
#  

 
#################### Export ###############################################
# ggsave("figs/indices_bar.png", plot = p1, width = 10, height = 10)
ggsave("figs/indices_box2.png", plot = p2, width = 8, height = 6)
 








###################### EcoQS#########################################

EcoQplot<- ggarrange(AllSAH[["AMBI"]][[2]], AllSAH[["BENTIX"]][[2]], AllSAH[["M.AMBI"]][[2]], 
                     common.legend = TRUE)


############################### correlation#############################

png("figs/corr.inds.png", width = 800, height = 1000)
pairs.panels(inds1[names(inds1) != "S2"], gap=0,pch=20,method="spearman", 
             ellipses = F, lm= T, density = F, hist.col ="#236573", 
             jiggle = F,rug=F,smoother=F,stars=F,ci=F,alpha=.05,)
dev.off()

png("figs/corr.inds.Log.png", width = 800, height = 1000)
pairs.panels(indsL1[names(indsL1) != "S2"], gap=0,pch=20,method="spearman", 
             ellipses = F, lm= T, density = F, hist.col ="#236573", 
             jiggle = F,rug=F,smoother=F,stars=F,ci=F,alpha=.05,)
dev.off()


