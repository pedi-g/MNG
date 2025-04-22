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
# env_raw<- read.csv("raw data/env.csv", row.names = 1)
# env<- env_raw[rownames(cwm),]
#ab <- read.csv("raw data/ab(for traits).csv",row.names = 1)
# ab_raw <- read.csv("raw data/ab.csv",row.names = 1)
 traits <- read.csv("output/comXtrait(LOGab).csv",row.names = 1)
# cwm<- traits
###### factors and data filtering
source('R-scripts/factors.R')

cwm1 <- cbind(cwm,factors)%>%
  filter(Habitat != "bch")
cwm<- cwm[row.names(cwm) %in% row.names(cwm1),]
factors<- factors[row.names(factors) %in% row.names(cwm1),]


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
#cwm1<- cwm%>% group_by(Season, Site, Habitat)
# ggplot(cwm1, aes(Weight,Habitat) )+
#   geom_boxplot()


################bar charts for CWM ##########################

nname<- c('Weight(gr)', 'Size', 'Feeding', 'Movement', 'Bioturbation',
          'Reporductive strategy','Life span', 'Toughness', 'Form', 'Flexibility')
trait<- c('W', 'S.', 'FD', 'MV', 'BT', 'RS', 'Ls', 'TN', 'FR', 'FL')
 ######### non fuzzy traits ##########

nonfuzz <- function(i){
t <- cwm1%>%
  group_by(Season, Site, Habitat) %>%
  select(a = starts_with(trait[i])) %>%
  summarize(mean = mean(a, na.rm =T),SD = sd(a, na.rm =T))

t$factors <- interaction(t$Habitat, t$Site, t$Season)
t$Season.Site <- interaction(t$Site, t$Season)
if(trait[i] != 'W'){
t$mean<- t$mean/3
t$SD<- t$SD/3
}


p <- ggplot(t, aes(x = Season.Site, y = t[[4]], fill = Habitat ))+
  geom_bar(position = "dodge",stat = "identity", col = "gray20")+
  geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                position = position_dodge2(width = 0.5, padding = 0.5)) +
  scale_fill_manual(values=cbPalette, labels = c("Mudflat", "Vegetated"))+
  labs(x = "", y = nname[i])+
  theme_bw(base_size = 11)+
 
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank(),#element_text(size = 12, face = "bold"),
        legend.key.size = unit(0.5, "cm"),
        axis.text = element_text(size = 10))
  
if(trait[i] != 'W'){
 p <-  p +
   scale_y_continuous(labels = scales::percent,
                      n.breaks = 6,
                      limits = (0:1))
}
  return(p)

}

                ########### fuzzy traits ##############

# function for selecting each function and summarizing them according to selected factors
fdbar<- function(i, Season = F, Site = F, Habitat = F){
  
  # factors of Season, Site, and Habitat are Boolean (TRUE or FALSE)
  # for trait insert the first initials of intended trait as a character string
    # for example if your trait is feeding group and in your table its modalities
    # all start with FD (e.g. FDP, FDG, etc.) but "FD" as input for the trait

  if (Season == T & Site == T & Habitat == T){
    fd <- cwm1 %>% group_by(Season.Site.Habitat)
  } else if (Season == T & Site == F & Habitat == T){
    fd <- cwm1 %>% group_by(Season.Habitat)
  } else if (Season == T & Site == T & Habitat == F){
    fd <- cwm1 %>% group_by(Season.Site)
  } else if (Season == F & Site == T & Habitat == T){
    fd <- cwm1 %>% group_by(Site.Habitat)
  } else if (Season == T & Site == F & Habitat == F){
    fd <- cwm1 %>% group_by(Season)
  } else if (Season == F & Site == T & Habitat == F){
    fd <- cwm1 %>% group_by(Site)
  } else if (Season == F & Site == F & Habitat == T){
    fd <- cwm1 %>% group_by(Habitat)
  } else {
    fd <- cwm1
  }
  
  fd<- fd %>% select(starts_with(trait[i]))
  # mods <- sum(!sapply(fd, is.factor))
  # 
  # fd[2:ncol(fd)]<- fd[2:ncol(fd)]/mods
  
  fd<- fd%>%
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
      labs(x = "", y = nname[i], fill = "") +
      theme_update(legend.position = "top")+
      theme(
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 10)
      )+
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
  } else{
    fd_l<- fd %>%
    pivot_longer(cols = starts_with(trait[i]), names_to = "traits", values_to = "value")%>%
    mutate(traits = factor(traits))
    fd_l$value<- fd_l$value/5
    
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
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))+
        #scale_fill_manual(values=cbPalette)
        scale_fill_brewer(palette="Spectral")
      
    } else {
    #### the real plot ###############
      plot<- ggplot(fd_l, aes(x = factor, y = value, fill = traits)) +
        geom_bar(stat = "identity", col = "gray20") +
        theme_bw(base_size = 11)+
        labs(x = "", y = nname[i], fill = "") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), legend.position = "top")+
        #scale_fill_manual(values=cbPalette)
        scale_fill_brewer(palette="Spectral")+
        theme(
          legend.text = element_text(size =10), #face = "bold"
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.5, "cm"))+
        scale_y_continuous(labels = scales::percent,
                           n.breaks = 6)
      }
  }
  
  return(list(fd_l, plot))
}

############ output ###############################################
w<- nonfuzz(1)
fl<- nonfuzz(10)
tn<- nonfuzz(8)

FD<- fdbar(3,T,T,T)[[2]]
# HB<- fdbar("HB",T,T,T)[[2]]
BT<- fdbar(5,T,T,T)[[2]]
# LD<- fdbar(7,T,T,T)[[2]]
RS<- fdbar(6,T,T,T)[[2]]
FR<- fdbar(9,T,T,T)[[2]]
MV<- fdbar(4,T,T,T)[[2]]
S<- fdbar(2,T,T,T)[[2]]
LS<- fdbar(7,T,T,T)[[2]]
# SD<- fdbar("SD",T,T,T)[[2]]

g1<- ggarrange(
                ggplot() + theme_void(),FD,
                ggplot() + theme_void(),
                nrow = 1, common.legend = F,
                widths = c(1,2,1))
g2<- ggarrange(MV,BT,RS,LS, nrow = 2, ncol = 2)

g3<- ggarrange( w,
                tn,
                # fl,
                # ggplot() + theme_void(),
                nrow = 1, common.legend = T
                # widths = c(1,2,1))
)
ggarrange(g3,g2,g1,
          # LD,RS,SD,
          ncol = 1, nrow = 3, heights = c(1,2,1))


# ggsave("figs/cwm_new.png", width = 11, height = 13)

ggsave("figs/cwm_new_-bch.png", width = 8, height = 13)
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


