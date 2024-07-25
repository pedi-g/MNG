#RDA for species and functions compared with environmental variables
#Pedram Ghahramani
#pedram.ghrm@gmail.com
#last edit 6-9-2024

#libraries
library(tidyverse)
library(vegan)
library(ggpubr)
#input ----
cwm <- read.csv("output/CWM.csv", row.names = 1)
#cwm.n <- read.csv("output/CWM(new).csv", row.names = 1)
env_raw<- read.csv("raw data/env.csv", row.names = 1)
env<- env_raw[rownames(cwm),]
#ab <- read.csv("raw data/ab(for traits).csv",row.names = 1)
ab_raw <- read.csv("raw data/ab.csv",row.names = 1)
#traits <- read.csv("output/comXtrait(LOGab).csv",row.names = 1)

###### factors 
source('R-scripts/factors.R')
cwm2 <- cbind(cwm,factors)

# color palettes ----
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with grey:
cbPalette <- c(  "#7570b3", "#1b9e77", "#d95f02", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())
# scatter plot of modalities compared with env data----
# par(mfrow = c(2, 2))
# par(mar = c(3, 4, 2, 1))
# plot(env$Sal, cwm$Weight, xlab = "salinity",
#      ylab = "CWM wight", pch = 20)
# plot(env$Sal, cwm$LD_direct, xlab = "salinity",
#      ylab = "CWM direct dev", pch = 20)
# plot(env$Sal, cwm$LD_lecitotrophic.larvae, xlab = "salinity",
#      ylab = "CWM lecitotrophic", pch = 20)
# plot(env$Sal, cwm$LD_planktotrophic.larvae, xlab = "salinity",
#      ylab = "CWM planktotrophik", pch = 20)


#boxplot----
#cwm2<- cwm%>% group_by(season, area, habitat)
# ggplot(cwm2, aes(Weight,habitat) )+
#   geom_boxplot()


################bar charts for CWM ##########################

 ######### nun fuzzy traits ##########

nonfuzz<- function(trait,a){
t<- cwm2%>% 
  group_by(season, area, habitat)%>%
  select(starts_with(trait)) %>%
  summarize(across(everything(), list(mean), .names = "{col}"))

t$factors <- interaction(t$habitat, t$area, t$season)
t$season.area <- interaction(t$season, t$area)

ggplot(t, aes(x = season.area, y = t[[4]], fill = habitat ))+
  geom_bar(position = "dodge",stat = "identity", col = "gray20")+
  scale_fill_manual(values=cbPalette)+
  labs(x = "", y = trait)+
  theme_bw(base_size = 15)
}

                ########### fuzzy traits ##############

# function for selecting each function and summarizing them according to selected factors
fdbar<- function(trait, season = F, area = F, habitat = F){
  
  # factors of season, area, and habitat are Boolean (TRUE or FALSE)
  # for trait insert the first initials of intended trait as a character string
    # for example if your trait is feeding group and in your table its modalities
    # all start with FD (e.g. FDP, FDG, etc.) but "FD" as input for the trait

  if (season == T & area == T & habitat == T){
    fd <- cwm2 %>% group_by(season.area.habitat)
  } else if (season == T & area == F & habitat == T){
    fd <- cwm2 %>% group_by(season.habitat)
  } else if (season == T & area == T & habitat == F){
    fd <- cwm2 %>% group_by(season.area)
  } else if (season == F & area == T & habitat == T){
    fd <- cwm2 %>% group_by(area.habitat)
  } else if (season == T & area == F & habitat == F){
    fd <- cwm2 %>% group_by(season)
  } else if (season == F & area == T & habitat == F){
    fd <- cwm2 %>% group_by(area)
  } else if (season == F & area == F & habitat == T){
    fd <- cwm2 %>% group_by(habitat)
  } else {
    fd <- cwm2
  }
  
  fd<- fd %>% select(starts_with(trait))%>%
    summarise(across(everything(), list(mean), .names = "{col}"))
 if (nrow(fd) != 1){
   names(fd)[1] <- "factor" 
 }
  
 if(ncol(fd) == 1){
    names(fd)[1] <- "factor"
    fd_l <- fd
    fd_l$n <- 1
    plot <- ggplot(fd_l, aes( x = 1, y = factor)) +
      geom_bar(stat = "identity", col = "gray20") +
      theme_bw(base_size = 16)+
      labs(x = "", y = "", fill = "") +
      theme_update(legend.position = "top")+
      theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 10)
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else{
    fd_l<- fd %>%
    pivot_longer(cols = starts_with(trait), names_to = "traits", values_to = "value")%>%
    mutate(traits = factor(traits))
    
    if(nrow(fd) == 1 ) {
      plot<- ggplot(fd_l, aes(x = "all", y = value, fill = traits)) +
        geom_bar(stat = "identity", col = "gray20") +
        theme_bw(base_size = 16)+
        theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 10)
      )+
        labs(x = "", y = "", fill = "") +
        theme_update(legend.position = "top") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        #scale_fill_manual(values=cbPalette)
        scale_fill_brewer(palette="Spectral")
      
    } else {
      plot<- ggplot(fd_l, aes(x = factor, y = value, fill = traits)) +
        geom_bar(stat = "identity", col = "gray20") +
        theme_bw(base_size = 13)+
        labs(x = "", y = "", fill = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")+
        #scale_fill_manual(values=cbPalette)
        scale_fill_brewer(palette="Spectral")+
        theme(
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 10)
        )
      }
  }
  
  return(list(fd_l, plot))
}

############ output ###############################################
w<- nonfuzz("w")
fl<- nonfuzz("FL")
tn<- nonfuzz("TN")

FD<- fdbar("FD",T,T,T)[[2]]
HB<- fdbar("HB",T,T,T)[[2]]
BT<- fdbar("BT",T,T,T)[[2]]
# LD<- fdbar("LD",T,T,T)[[2]]
RT<- fdbar("RT",T,T,T)[[2]]
FR<- fdbar("FR",T,T,T)[[2]]
MV<- fdbar("MV",T,T,T)[[2]]
# RS<- fdbar("RS",T,T,T)[[2]]
# SD<- fdbar("SD",T,T,T)[[2]]
g1<- ggarrange(w,fl,tn, nrow = 1, common.legend = T)
g2<- ggarrange(FD,HB,BT,RT,FR,MV)
ggarrange(g1,g2,
          # LD,RS,SD,
          ncol = 1, nrow = 2, heights = c(1,2.5))

ggsave("figs/cwm_all.png", width = 11, height = 13)

# ########## non fuzzy 
# tiff("figs/CWM-TN-TTT.tiff",  width=1000, height=800, res=180)
#  tn
#  dev.off()
# 
#  tiff("figs/CWM-fl-TTT.tiff",  width=1500, height=1000, res=180)
#  fl
#  dev.off()
# 
#  tiff("figs/CWM-W-TTT.tiff",  width=1000, height=800, res=180)
#  w
#  dev.off()


# ##################### fuzzy 

#  tiff("figs/CWM-FD-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("FD",T,T,T)
#  dev.off()
 # tiff("figs/CWM-MV-TTT.tiff",  width=1000, height=800, res=180)
 # pm
 # dev.off()
#  tiff("figs/CWM-HB-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("HB",T,T,T)
#  dev.off()
#  tiff("figs/CWM-BT-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("BT",T,T,T)
#  dev.off()
#  tiff("figs/CWM-LD-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("LD",T,T,T)
#  dev.off()
#  tiff("figs/CWM-RT-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("RT",T,T,T)
#  dev.off()
#  tiff("figs/CWM-FR-TTT.tiff",  width=1000, height=800, res=180)
#  fdbar("FR",T,T,T)
#  dev.off()


# tiff("figs/CWM-FR-TFF.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",T,F,F)
# dev.off()
# tiff("figs/CWM-FR-FTF.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",F,T,F)
# dev.off()
# tiff("figs/CWM-FR-FFT.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",F,F,T)
# dev.off()
# tiff("figs/CWM-FR-TTF.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",T,T,F)
# dev.off()
# tiff("figs/CWM-FR-TFT.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",T,F,T)
# dev.off()
# tiff("figs/CWM-FR-FTT.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",F,T,T)
# dev.off()
# tiff("figs/CWM-FR-FFF.tiff",  width=1500, height=1000, res=200)
# fdbar("FR",F,F,F)
# dev.off()


