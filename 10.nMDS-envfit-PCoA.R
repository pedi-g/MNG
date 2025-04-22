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
# sp.bray<- read.csv("output/dist.spxtrait(bray).csv", row.names = 1) #bray curtis
#community-trait matrix and distance
com.trait<- read.csv("output/comxtrait(LOGab).csv", row.names = 1)
com.gow<- read.csv("output/dist.comXtrait(gowdis).csv", row.names = 1) #gower
# com.bray<- read.csv("output/dist.comxtrait(bray).csv", row.names = 1) #bray curtis
#community abundance matrix
ab_raw<- read.csv("raw data/ab.csv", row.names = 1)
ab<- ab_raw[,-((ncol(ab_raw)-3):ncol(ab_raw))]
ab_raw<- cbind(ab, factors)
#func group abundance matrix
abg<- read.csv("output/ab_grouped.csv", row.names = 1)

#environmental data
env<- read.csv("raw data/env.csv", row.names = 1)
#community Weighted Mean
cwm <- read.csv("output/CWM.csv", row.names = 1)
#biomass
mass<- read.csv("raw data/mass.csv", row.names = 1)


#selecting Habitats
factors<- factors%>% filter(Habitat != "bch")
ab_raw<- ab_raw[row.names(ab_raw) %in% row.names(factors),]
abg<- abg[row.names(abg) %in% row.names(factors),]
ab<- ab_raw%>%
  select(-names(factors))
sp.dis <- vegdist(log(1+ab))

mass<- mass[row.names(mass) %in% row.names(factors),]
mass_raw<- cbind(mass, factors)


# com.bray<- com.bray[rownames(com.bray) %in% rownames(factors),
                    # colnames(com.bray) %in% rownames(factors)]
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
  theme( legend.title = element_text(size=9), legend.position = "right",
           legend.text = element_text(size = 10),
           axis.title = element_text(size= 12),
           axis.text = element_text(size = 10))+
  theme(legend.title = element_text(face = "bold"))+
  theme(axis.title.y = element_text(vjust = 1))+ # increase distance from the y-axis
  theme(axis.title.x = element_text(vjust = -1)) # increase distance from the x-axis
  # theme(legend.position = c(0.93, 0.3))

plot1<- function(){
  ggplot(mds, aes(MDS1,MDS2))+
    geom_point(aes(shape = Season.Site , color = Habitat),stroke = 1.7, alpha = 0.8, size = 2)+
    scale_color_manual( values=c("#264653", "#2a9d8f", "#e9c46a"))+
#                        labels = c("Creek","Beach","Mangrove"))+
    scale_shape_manual(values = c(17, 16, 2, 1))+
    geom_text(x= Inf, y = Inf,vjust="inward",
              hjust="inward",label = paste("stress:",stress))+
    labs(x= "MDS1", y= "MDS2", shape = "Site.Season", colour = "Habitat")+
    my_theme
    #geom_text(x= max(mds$MDS1)-0.25,y = max(mds$MDS2) ,label = paste("stress:",stress))
  #stat_ellipse(aes(shape = Season.Site, color = Season.Site), color= "gray50",type = "t", linetype= 1)+
  #stat_ellipse(aes(shape = Season.Site, color = Season.Site), geom="polygon",alpha= 0.08 , color= 'gray',type = "t")
}

  # a plot theme for two factors of Season and Habitat (use for Sites separatly)
plot2<- function(){
  ggplot(mds, aes(MDS1,MDS2))+
    geom_point(aes(shape = Season , color = Habitat), size = 1.5)+
    scale_color_manual( values=c("#264653", "#2a9d8f", "#e9c46a"))+
    scale_shape_manual(values = c(16, 17))+
    geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",
              label = paste("stress:",stress))+
    labs(x= "nMDS1", y= "nMDS2")+
    my_theme
  #geom_text(x= max(mds$MDS1)-0.25,y = max(mds$MDS2) ,label = paste("stress:",stress))
}


###############################  NMDS ########################################

############## nmds of community-trait composition ########


        ###### based on gower distance (gowdis function) #
  nmds<- metaMDS(com.gow)
  mds<- cbind(as.data.frame(nmds$points), factors)
  stress<- round(nmds$stress,3)
  
  # g1 <- ggplot(mds, aes(MDS1, MDS2, color = Season))+
  #   geom_point()
  # g2 <- ggplot(mds, aes(MDS1, MDS2, color = Site))+
  #   geom_point()
  # g3 <- ggplot(mds, aes(MDS1, MDS2, color = Habitat))+
  #   geom_point()
  # 
  (g4<- plot1()+
      stat_ellipse(aes(fill = Site), color = "gray0", alpha= 0.2, type = "t", level = 0.9))
  
  # (gower<- ggarrange(g1, g2, g3, g4))
  
############## nmds of community-species composition ########

# nmds <- metaMDS(ab, distance = "bray")
# mds<- cbind(as.data.frame(nmds$points), factors)
# stress<- round(nmds$stress,3)
# (sp2 <- plot1()+
#     stat_ellipse(aes(fill = Site), color = "gray0", alpha= 0.2, type = "t", level = 0.9))


#another way
nmds<- metaMDS(sp.dis)
mds<- cbind(as.data.frame(nmds$points), factors)
#scores(nmds , "species")
stress<- round(nmds$stress,3)
(sp <- plot1()+
    stat_ellipse(aes(fill = Site), color = "gray0", alpha= 0.2, type = "t", level = 0.9))

  #ggtitle("species composition")
 
  #stat_ellipse(aes(shape = Season.Site, color = Season.Site), color= "gray50",type = "t", linetype= 1)+
  #stat_ellipse(aes(shape = Season.Site, color = Season.Site), geom="polygon",alpha= 0.08 , color= 'gray',type = "t")

################### nmds of biomass composition ##############
# mass_dis <- vegdist(log(1+mass))
# nmds<- metaMDS(mass_dis)
# mds<- cbind(as.data.frame(nmds$points), factors)
# #scores(nmds , "species")
# stress<- round(nmds$stress,3)
# (p_mass <- plot1()+
#     stat_ellipse(aes(fill = Site), color = "gray0", alpha= 0.2, type = "t", level = 0.9))

  




################# sp and mass together #############

# space<- ggplot()+theme_tufte() # a third, blanck plot to increase the space
# (two1<- ggarrange(sp,space,p_mass, common.legend = T, labels = c('a','','b'),
#                   ncol = 3,
#                   nrow = 1,
#                   label.x = 0.12,
#                   label.y = 1,
#                   hjust = -1.5,
#                   vjust = 1.5,
#                   font.label = list(size = 15, color = "gray30",
#                                     face = "bold", family = NULL),
#                   align = c("none", "h", "v", "hv"),
#                   widths = c(1,0.09,1),
#                   heights = 1,
#                   legend = "right",
#                   legend.grob = NULL))



#-------------------------nmds of species-trait distance ---------------------

# #gower (gowdis)
# nmds<- metaMDS(sp.gow)
# mds<- cbind(as.data.frame(nmds$points))
# stress<- round(nmds$stress,3)
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



################################ PCoA ######################################### 
 #  pco <- dudi.pco(as.dist(com.gow), scannf = F, nf = 4)
 # scatter(pco) ### same info as nMDS --> no need for this
 # sum(pco$eig[1:2]) / sum(pco$eig)
 # 
 # cor(pco$li[, 1:2], com.trait)
 # 
 # scatter(dudi.pca(com.trait, scannf = F, nf = 4))
 
 ######################## envfit ###########################

env1 <- select(env, -c(DO, sand, clay, silt, TOM))
env1 <- na.omit(env1[!rownames(env1)=="C611_w",])
env1 <- scale(env1)
env2 <- select(env, -c(DO, sand, Temp, Sal, pH))
env2 <- na.omit(env2[!rownames(env2)=="C611_w",])
env2 <- scale(env2)
 # function for the envfit plot
 envfitplot<-function(df, size = 0.5){
   
   nmds0<- metaMDS(df)
   stress<- round(nmds0$stress,3)
   mds0<- cbind(as.data.frame(nmds0$points), factors)
   
   p<- ggplot(mds0, aes(MDS1,MDS2))+
     geom_point(aes(shape = Season.Site , color = Habitat),
                stroke = 1.7, alpha = 0.8, size = 2)+
     scale_color_manual( values=c("#264653", "#2a9d8f", "#e9c46a"))+
     #                        labels = c("Creek","Beach","Mangrove"))+
     scale_shape_manual(values = c(17, 16, 2, 1))+
     geom_text(x= Inf, y = Inf,vjust="inward",hjust="inward",
               label = paste("stress:",stress))+
     labs(x= "MDS1", y= "MDS2", shape = "Site.Season (Shape)", colour = "Habitat (color)")+
     my_theme+
     stat_ellipse(aes(fill = Site), color = "gray0", alpha= 0.2,
                  type = "t", level = 0.9)
   
   p
   
   
   # Define a function to process each subset
   process_env <- function(vars) {
     env1 <- select(env, all_of(vars))              # Select specific columns
     env1 <- na.omit(env1[rownames(env1) != "C611_w", , drop = FALSE]) # Remove NA and specific row
     env1 <- scale(env1)                         # Scale the data
     
     com1 <- df[rownames(df) %in% rownames(env1), colnames(df) %in% rownames(env1)]
     nmds1 <- metaMDS(com1)
     fit1 <-  envfit(nmds1, env1, perm = 9999) #correlation with env data
     return( fit1)
   }
   
   TOM<- process_env( "TOM")
   water<- process_env(c("Temp", "Sal", "pH"))
   sed<- process_env(c('silt','clay'))
   
   spp.scrs <- rbind(size*as.data.frame(scores(TOM, display = "vectors")),
                     size*as.data.frame(scores(water, display = "vectors")),
                     size*as.data.frame(scores(sed, display = "vectors")))
   # Define the column groups for each subset
   
   fit<-rbind(cbind(r2=TOM[['vectors']][["r"]],Pr = TOM[["vectors"]][["pvals"]]),
              cbind(r2=water[['vectors']][["r"]],Pr = water[["vectors"]][["pvals"]]),
              cbind(r2=sed[['vectors']][["r"]],Pr = sed[["vectors"]][["pvals"]]))
   
   # Add environmental vectors
   p<- p + 
     geom_segment(data = spp.scrs, aes(x = 0, y = 0,
                                       xend = NMDS1, yend = NMDS2),
                  arrow = arrow(length = unit(0.2, "cm")),
                  color = "gray10",size = 0.6) +
     # Highlight significant vectors (assuming p.max = 0.05 means significant vectors)
     geom_segment(data = subset(spp.scrs, fit[,2] <= 0.05),
                  aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                  arrow = arrow(length = unit(0.2, "cm")),
                  color = "red2", size = 0.6) 
     # Add variable names at a fixed distance from the tip of the arrows
     # geom_text(data = spp.scrs,
     #           aes(x = NMDS1 * 1.1, y = NMDS2 * 1,1, label = rownames(spp.scrs)),
     #           vjust = 0.5, hjust = 0.5, color = "gray20")
     # annotate("text", x = (spp.scrs$NMDS1*X),
     #          y = (spp.scrs$NMDS2*Y),
     #          label = rownames(spp.scrs), color = "gray20")
   
   return(list(fit, p, spp.scrs))
  }
 

 
 (sp<- envfitplot(as.data.frame(as.matrix(sp.dis)), 0.5))
 (g4<- envfitplot(com.gow, 0.4))
 
 
 (sp[[2]]<- sp[[2]]+ 
   annotate("text",
            x = (sp[[3]]$NMDS1*c(1.3,1.4,2.0,1.0,1.2,1.2)),
            y = (sp[[3]]$NMDS2*c(1.2,1.0,1.2,1.2,1.1,0.8)),
            label = rownames(sp[[3]]), color = "gray40", size = 4,fontface = "bold"))

(g4[[2]]<- g4[[2]]+ 
     annotate("text",
              x = (g4[[3]]$NMDS1*c(1.3,1.1,1.1,1.1,1.2,1.8)),
              y = (g4[[3]]$NMDS2*c(1.3,1.15,1.05,1.1,1.2,1.1)),
              label = rownames(g4[[3]]), color = "gray40",size = 4,fontface = "bold"))
 
g4[[1]]
sp[[1]]
 ###############sp and trait compositon together#############
 # (two<- g4+labs(tag = "a")| sp+ labs(tag = "b"))
 space<- ggplot()+theme_tufte() # a third, blanck plot to increase the space
 (two<- ggarrange(sp[[2]],space,g4[[2]], common.legend = T, labels = c('a','','b'),
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
 
 
 ############################## RDA #############################################
 # env1 <- na.omit(env[!rownames(env)=="C611_w",])%>%
 #   select(-DO)
 # env1 <- scale(env1)
 # abg2 <- abg[rownames(env1),]%>%
 #   select(starts_with("G"))
 # env1 <- as.data.frame(env1)
 # rda1 <- rda(abg2 ~ ., data = env1, scale = TRUE)
 # plot(rda1, display = c("bp", "sp"))
 # 
 # #RDA for cwm relation with env data
 # 
 # cwm2<- filter(cwm, rownames(cwm) %in% rownames(env1))%>%
 #   select(starts_with(c('S.','FD','MV','BT','TN','FL','RT','LS')))
 # traits2 <- filter(com.trait, rownames(com.trait) %in% rownames(env1))
 # ab2 <- filter(ab, rownames(ab) %in% rownames(env1))
 # 
 # par(mfrow = c(1, 1))
 # rda.cwm <- rda(cwm2 ~ ., data = env1)
 # 
 # png("figs/RDA.png", width = 1000, height = 900,pointsize = 25)
 # plot(rda.cwm, type = "n", scaling = "sites")
 # text(rda.cwm, dis = "cn", scaling = "sites")
 # text(rda.cwm, display =  "sp", scaling = "sites", col = "darkblue", cex = 1 ) #modalities are shown as 
 # dev.off()
 
 
 ##output####
 ggsave("figs/nMDS.triat(gower)&sp.png", plot = two, width = 12, height = 5) 
 # ggsave("figs/nMDS.triat.bray.4.png", plot = bray , width = 12, height = 8)  
 # ggsave("figs/nMDS.triat.bray.png", plot = b4 , width = 8, height = 5) 
 # ggsave("figs/nMDS.trait.gower.4.png", plot = gower , width = 12, height = 8)  
 ggsave("figs/nMDS.trait.gower.png", plot = g4[[2]], width = 8, height = 5)
 # ggsave("figs/nMDS.sp(nolog).png", plot = sp2, width = 8, height = 5)
 ggsave("figs/nMDS.spLog.png", plot = sp[[2]], width = 8, height = 5)
 #ggsave("figs/sp.nMDS.Site.png", plot = p1, width = 4.5, height = 3)
 #ggsave("figs/trait.nMDS.Site.png", plot = p2, width = 4.5, height = 3)
 
 #ggsave("penguin_plot.pdf", dpi = 600, width = 100, height = 60, unit = "mm")
 
 #pdf("ggplot-cmyk.pdf", width = 12 / 2.54, height = 8 / 2.54, colormodel = "cmyk")
 
 #print(penguin_plot)
 
 #dev.off()


