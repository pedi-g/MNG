#nMDS ordination diagrams
#Pedram Ghahramani
#pedram.ghrm@gmail.com
#last mod: 2-28-2024

#####libraries ----
library(tidyverse)
library(ggthemes)
library(vegan)
library(FD)
library(ggfortify)
library(ggpubr)
library(patchwork)

#####input and data prepration----

source('R-scripts/factors.R') #recalling factors 
#species-trait distance 
sp.gow<- read.csv("output/dist.spxtrait(gowdis).csv", row.names = 1) #gower
sp.bray<- read.csv("output/dist.spxtrait(bray).csv", row.names = 1) #bray curtis
#community-trait matrix and distance
com.trait<- read.csv("output/comxtrait(LOGab).csv", row.names = 1)
com.gow<- read.csv("output/dist.comXtrait(gowdis).csv", row.names = 1) #gower
com.bray<- read.csv("output/dist.comxtrait(bray).csv", row.names = 1) #bray curtis
#community abundance matrix
ab_raw<- read.csv("raw data/ab.csv", row.names = 1)
ab<- ab_raw[,-((ncol(ab_raw)-3):ncol(ab_raw))]
ab_raw<- cbind(ab, factors)
#environmental data
env<- read.csv("raw data/env.csv", row.names = 1)
#community Weighted Mean
cwm <- read.csv("output/CWM.csv", row.names = 1)



#selecting habitats
factors<- factors%>% filter(habitat != "bch")
ab_raw<- ab_raw[row.names(ab_raw) %in% row.names(factors),]
ab<- ab_raw%>%
  select(-names(factors))
com.bray<- com.bray[rownames(com.bray) %in% rownames(factors),
                    colnames(com.bray) %in% rownames(factors)]
com.gow<- com.gow[rownames(com.gow) %in% rownames(factors),
                    colnames(com.gow) %in% rownames(factors)]
com.trait<- com.trait[rownames(com.trait) %in% rownames(factors),]
env<- env[rownames(env) %in% rownames(ab_raw),]
cwm<- cwm[rownames(cwm) %in% rownames(ab_raw),]

############### we define plot functions to prevent repetition#####
  # a plot theme for all three factors in use

my_theme<-
  theme_classic()+#base_size = 12)+ #another way to detemine the text size
  theme(panel.grid.major = element_blank(), # remove the major lines
          panel.grid.minor = element_blank())+ # remove the minor lines
  theme( legend.title = element_text(size=8), legend.position = "right",
           legend.text = element_text(size = 9),
           axis.title = element_text(size= 12),
           axis.text = element_text(size = 10))+
  theme(legend.title = element_text(face = "bold"))+
  theme(axis.title.y = element_text(vjust = 1))+ # increase distance from the y-axis
  theme(axis.title.x = element_text(vjust = -1)) # increase distance from the x-axis
  # theme(legend.position = c(0.93, 0.3))

plot1<- function(){
  ggplot(mds, aes(MDS1,MDS2))+
    geom_point(aes(shape = season.area , color = habitat),stroke = 1.7, alpha = 0.8, size = 2)+
    scale_color_manual( values=c("#264653", "#2a9d8f", "#e9c46a"))+
#                        labels = c("Creek","Beach","Mangrove"))+
    scale_shape_manual(values = c(17, 2, 16, 1))+
    geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",label = paste("stress:",stress))+
    labs(x= "MDS1", y= "MDS2", shape = "", colour = "")+
    my_theme
    #geom_text(x= max(mds$MDS1)-0.25,y = max(mds$MDS2) ,label = paste("stress:",stress))
  #stat_ellipse(aes(shape = season.area, color = season.area), color= "gray50",type = "t", linetype= 1)+
  #stat_ellipse(aes(shape = season.area, color = season.area), geom="polygon",alpha= 0.08 , color= 'gray',type = "t")
}

  # a plot theme for two factors of season and habitat (use for areas separatly)
plot2<- function(){
  ggplot(mds, aes(MDS1,MDS2))+
    geom_point(aes(shape = season , color = habitat), size = 1.5)+
    scale_color_manual( values=c("#264653", "#2a9d8f", "#e9c46a"))+
    scale_shape_manual(values = c(16, 17))+
    geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",label = paste("stress:",stress))+
    labs(x= "nMDS1", y= "nMDS2")+
    my_theme
  #geom_text(x= max(mds$MDS1)-0.25,y = max(mds$MDS2) ,label = paste("stress:",stress))
}


###############################  NMDS ########################################

############## nmds of community-trait composition ########

          ###### based on **bray curtis** dissimilarity ####

    # all plots
# nmds
# nmds<- metaMDS(com.bray)
# mds<- cbind(as.data.frame(nmds$points), factors)
# stress<- round(nmds$stress,3)
# #mds<- mds %>% filter(area == "Tis")
# #plotting
# (b1 <- ggplot(mds, aes(MDS1, MDS2, color = season))+
#   geom_point())
# (b2 <- ggplot(mds, aes(MDS1, MDS2, color = area))+
#   geom_point())
# (b3 <- ggplot(mds, aes(MDS1, MDS2, color = habitat))+
#   geom_point())
# (b4<- plot1())
# #(p4 <- plot2())#+ggtitle("all"))
# (bray<- ggarrange(b1, b2, b3, b4))
# 


        ###### based on gower distance (gowdis function) #####
  nmds<- metaMDS(com.gow)
  mds<- cbind(as.data.frame(nmds$points), factors)
  stress<- round(nmds$stress,3)
  
  g1 <- ggplot(mds, aes(MDS1, MDS2, color = season))+
    geom_point()
  g2 <- ggplot(mds, aes(MDS1, MDS2, color = area))+
    geom_point()
  g3 <- ggplot(mds, aes(MDS1, MDS2, color = habitat))+
    geom_point()
  
  (g4<- plot1()+
      stat_ellipse(aes(fill = area), color = "gray0", alpha= 0.2, type = "t", level = 0.9))
  
  (gower<- ggarrange(g1, g2, g3, g4))
  
############## nmds of community-species composition ########

# nmds <- metaMDS(ab, distance = "bray")
# mds<- cbind(as.data.frame(nmds$points), factors)
# stress<- round(nmds$stress,3)
# (sp2 <- plot1()+
#     stat_ellipse(aes(fill = area), color = "gray0", alpha= 0.2, type = "t", level = 0.9))


#another way
sp_dis <- vegdist(log(1+ab))
nmds<- metaMDS(sp_dis)
mds<- cbind(as.data.frame(nmds$points), factors)
#scores(nmds , "species")
stress<- round(nmds$stress,3)
(sp <- plot1()+
    stat_ellipse(aes(fill = area), color = "gray0", alpha= 0.2, type = "t", level = 0.9))


  #ggtitle("species composition")
 
  #stat_ellipse(aes(shape = season.area, color = season.area), color= "gray50",type = "t", linetype= 1)+
  #stat_ellipse(aes(shape = season.area, color = season.area), geom="polygon",alpha= 0.08 , color= 'gray',type = "t")

  ###############sp and trait compositon together##########################
# (two<- g4+labs(tag = "a")| sp+ labs(tag = "b"))
space<- ggplot()+theme_tufte() # a third, blanck plot to increase the space
(two1<- ggarrange(sp,space,g4, common.legend = T, labels = c('a','','b'),
                  ncol = 3,
                  nrow = 1,
                  label.x = 0.12,
                  label.y = 1,
                  hjust = -1.5,
                  vjust = 1.5,
                  font.label = list(size = 15, color = "gray30",
                                    face = "bold", family = NULL),
                  align = c("none", "h", "v", "hv"),
                  widths = c(1,0.09,1),
                  heights = 1,
                  legend = "right",
                  legend.grob = NULL))


#-------------------------nmds of species-trait distance ---------------------
# 
# # bray
# nmds<- metaMDS(sp.bray)
# mds<- cbind(as.data.frame(nmds$points))
# 
# ggplot(mds, aes(-MDS1,MDS2))+
#   theme_classic()+
#   labs(x= "nMDS1", y= "nMDS2")+
#   geom_point(size = 2)+
#   ggtitle("sp bray")+
#   theme( legend.title = element_text(size=7) ,legend.position = "none",
#          #legend.text = element_text(size = 6),
#          axis.title = element_text(size= 14),
#          axis.text = element_text(size = 12))+
#   geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",label = paste("stress:",stress))
# 
# 
# 
# #gower (gowdis)
# nmds<- metaMDS(sp.gow)
# mds<- cbind(as.data.frame(nmds$points))
# 
# ggplot(mds, aes(-MDS1,MDS2))+
#   theme_classic()+
#   labs(x= "nMDS1", y= "nMDS2")+
#   geom_point(size = 2)+
#   ggtitle("sp gowdis")+
#   theme( legend.title = element_text(size=7) ,legend.position = "none",
#          #legend.text = element_text(size = 6),
#          axis.title = element_text(size= 14),
#          axis.text = element_text(size = 12))+
#   geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",label = paste("stress:",stress))
# 

#----------------from imported coordinationd (PRIMER) --------------------

nmds2<- read.csv("raw data/nmds.traits.raw.csv", row.names = 1)%>%
  filter(habitat != "bch")
nmds1<- read.csv("raw data/nmds.raw.csv", row.names = 1)%>%
  filter(habitat != "bch")

#sp comp 
MDS1<- ggplot(nmds1, aes(x,y, color= habitat, shape= area))+
  geom_point(size= 2)+
  theme_classic()+
  labs(x= "nMDS1", y= "nMDS2")+
scale_color_colorblind()+
  theme( legend.title = element_text(size=10) ,#legend.position = "top" ,
        legend.text = element_text(size = 8),
        axis.title = element_text(size= 9))

  p1<- 
    MDS1+# ggtitle("Species Composition")+
  stat_ellipse(geom="polygon",alpha= 0.07 , color= 'gray',type = "t")+
  stat_ellipse(color= "gray50",type = "t", linetype= 1)

 p1
  
 #trait comp
 MDS2<- ggplot(nmds2, aes(x,y, color= habitat, shape= area))+
   geom_point(size= 2)+
   #geom_point(data= nmds2[nmds2$season=="s",],color= "white",size=1.5)+
   theme_classic()+
   labs(x= "nMDS1", y= "nMDS2")+
   scale_color_colorblind()+
   theme( legend.title = element_text(size=10) ,#legend.position = "top" ,
          legend.text = element_text(size = 8),
          axis.title = element_text(size= 9))

 p2<- 
   MDS2+
   stat_ellipse(geom="polygon",alpha= 0.08 , color= 'gray',type = "t")+
   stat_ellipse(color= "gray50",type = "t", linetype= 1)
 
 p2
  #stat_ellipse(data= nmds[nmds$habitat!="beach"|nmds$area!= "Gwatr",], aes(x, y, shape= area),color= "gray", inherit.aes = F,type = "t", linetype= 1)+ #excluding gwtr beach stations for ellipse
  #stat_ellipse(aes(x,y,fill = habitat),geom= "polygon",alpha=0.07,type = "t", linetype= 1, inherit.aes = F)+
  #stat_ellipse(aes(x,y,color= habitat),type = "t", level= 0.8,linetype= 3, inherit.aes = F)
#stat_ellipse(aes(shape=season),type = "norm", linetype= 2)
 
################################ PCoA ######################################### 
  pco <- dudi.pco(as.dist(com.gow), scannf = F, nf = 4)
 scatter(pco) ### same info as nMDS --> no need for this
 sum(pco$eig[1:2]) / sum(pco$eig)
 
 cor(pco$li[, 1:2], com.trait)
 
 scatter(dudi.pca(com.trait, scannf = F, nf = 4))
 
 ######################## envfit ###########################
 # work on this later
 env1 <- env[!rownames(env)=="C611_w"&!is.na(env$Temp),]
 env1 <- scale(dplyr::select(env1, -DO))
 com <- com.gow %>% filter( rownames(com.gow) %in% rownames(env1))
 com <- com %>% dplyr::select(rownames(env1))
 
 nmds <- metaMDS(com)
 fit <-  envfit(nmds, env1, perm = 999) #correlation with env data
 spp.scrs <- as.data.frame(scores(fit, display = "vectors"))
 scrs <- as.data.frame(scores(nmds, display = "sites"))
 
 
 plot(nmds)
 plot(fit, add = T)
 plot(fit, p.max = 0.05, col = "red")
 
 
 ############################## RDA #############################################
 env1 <- env[!is.na(env$Temp),]%>%
   select(-DO)
 env1 <- scale(env1)
 ab <- ab[rownames(env1),]
 env1 <- as.data.frame(env1)
 rda1 <- rda(ab ~ ., data = env1, scale = TRUE)
 plot(rda1, display = c("bp", "sp"))
 
 #RDA for cwm relation with env data----
 #omiting missing values 
 env3<- na.omit(env)
 env2<- na.omit(filter(env ,grepl("w", rownames(env))))
 
 cwm2<- filter(cwm, rownames(cwm) %in% rownames(env2))
 traits2 <- filter(com.trait, rownames(com.trait) %in% rownames(env2))
 ab2 <- filter(ab, rownames(ab) %in% rownames(env2))
 
 par(mfrow = c(1, 1))
 rda.cwm <- rda(cwm2 ~ ., data = env2)
 
 png("figs/RDA.png", width = 1000, height = 900,pointsize = 25)
 plot(rda.cwm, type = "n", scaling = "sites")
 text(rda.cwm, dis = "cn", scaling = "sites")
 text(rda.cwm, display =  "sp", scaling = "sites", col = "darkblue", cex = 0.7 ) #modalities are shown as 
 dev.off()
 
 
 ##output####
 ggsave("figs/nMDS.triat(gower)&sp.png", plot = two1, width = 12, height = 5) 
 # ggsave("figs/nMDS.triat.bray.4.png", plot = bray , width = 12, height = 8)  
 # ggsave("figs/nMDS.triat.bray.png", plot = b4 , width = 8, height = 5) 
 ggsave("figs/nMDS.trait.gower.4.png", plot = gower , width = 12, height = 8)  
 ggsave("figs/nMDS.trait.gower.png", plot = g4, width = 8, height = 5)
 # ggsave("figs/nMDS.sp(nolog).png", plot = sp2, width = 8, height = 5)
 ggsave("figs/nMDS.spLog.png", plot = sp, width = 8, height = 5)
 #ggsave("figs/sp.nMDS.area.png", plot = p1, width = 4.5, height = 3)
 #ggsave("figs/trait.nMDS.area.png", plot = p2, width = 4.5, height = 3)
 
 #ggsave("penguin_plot.pdf", dpi = 600, width = 100, height = 60, unit = "mm")
 
 #pdf("ggplot-cmyk.pdf", width = 12 / 2.54, height = 8 / 2.54, colormodel = "cmyk")
 
 #print(penguin_plot)
 
 #dev.off()


