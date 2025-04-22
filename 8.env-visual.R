#Pedram Ghahramani
#pedram.ghrm@gmail.com
#env data analysis and ploting
#3-6-2024
#libraries----------
 library(tidyverse)
 library(ggfortify) #for autoplot of PCA
 library(ggpubr)
 library(psych)     #for correlation pair panels
 library(extrafont) #for fonts

font_import()
loadfonts(device="win")     #Register fonts for Windows bitmap output
# fonts()           
theme_set(theme_bw())
#input-----------
env <- read.csv("raw data/env2.csv", row.names = 1)%>%
   filter(Habitat!= "bch")
env1<- env%>%
  select(-Habitat, -Site, -Season)
#defining factor variables--------

# env$station <- factor(env$station)
env$Season <- factor(env$Season)
env$Habitat<- factor(env$Habitat, levels= c("Mud", "Veg", "bch"))
env$Site <- factor (env$Site, levels = c("Plt", "Nat"))
env$Site.Habitat<- interaction(env$Site, env$Habitat)
env$Season.Habitat<- interaction(env$Season, env$Habitat)
env$Season.Site <- interaction(env$Season, env$Site)
env$factors <- interaction(env$Season, env$Site, env$Habitat)
# env$transect <- factor(env$transect)

facts<- env[9:15]


#correlations####
png("figs/corr.env.pairs.png", width = 900, height = 900, pointsize = 20)
pairs.panels(env1[-which(colnames(env1) == "DO")], gap=0,pch=20,method="spearman", 
             ellipses = F, lm= T, density = F, hist.col ="#236573", 
             jiggle = F,rug=F, smoother=F,stars=F,ci=F,alpha=.05)
dev.off()


############################# mean and SD ###############################
ave <- env %>% group_by(Season, Site, Habitat)%>%
  summarise(Temp = mean(Temp, na.rm = T),Sal =  mean(Sal, na.rm = T), DO = mean(DO, na.rm = T),
            pH = mean(pH, na.rm = T), Clay = mean(clay, na.rm = T), Silt =  mean(silt, na.rm = T),
            Sand = mean(sand, na.rm = T), TOM = mean(TOM, na.rm = T))
std <- env %>% group_by(Season, Site, Habitat) %>%
  summarise(Temp_sd = sd(Temp, na.rm = T),Sal_sd = sd(Sal, na.rm = T),#/ sqrt(sum(!is.na(Sal))), 
            DO_sd = sd(DO, na.rm = T), pH_sd = sd(pH, na.rm = T),
            Clay_sd = sd(clay, na.rm = T), Silt_sd =  sd(silt, na.rm = T),
            Sand_sd = sd(sand, na.rm = T), TOM_sd = sd(TOM, na.rm = T))

av_sd <- merge(ave,std)

factors<- av_sd[, 1:3]
factors$Site.Habitat<- interaction(factors$Site, factors$Habitat)

factors$Season.Site <- interaction(factors$Season, factors$Site)
av_sd <- cbind(av_sd, factors)
av_sd$factors <- paste(factors$Season, factors$Site.Habitat, sep = "_")
av_sd <- av_sd[,c(4:25)]

########################### outliars #################################

# boxplot by variable and factor
box<- function(factor, variable, col = "Habitat"){
  ggboxplot(env, factor, variable, color = col)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

box("factors", "Temp")
box("factors", "Sal")
box("factors", "pH")
box("factors", "TOM")
box("factors", "DO")

    ## Function to detect outliars using Z-scores
out <- function(factor, threshold = 3) {
  groups<- unique(env[[factor]])
  outliars<- list()
  
  for (i in groups){
    outliar <- env[env[[factor]] == i,]%>%
      select(1:8)%>%
      mutate(across(everything(), ~ ifelse(abs(scale(.)) > threshold, ., NA)))
    
    outliars[[as.character(i)]] <- outliar
    
    # Get non-NA values and their indices
    non_na_indices <- which(!is.na(outliar), arr.ind = TRUE)
    non_na_values <- data.frame(
      row = rownames(outliar)[non_na_indices[, 1]],
      column = colnames(outliar)[non_na_indices[, 2]],
      value = outliar[non_na_indices])
    
    outliars[[paste(as.character(i),"1")]] <- non_na_values
    
  }
  
  # all of the variables across all stations
  all<- env%>%
    select(1:8)%>%
    mutate(across(everything(), ~ ifelse(abs(scale(.)) > threshold, ., NA)))
  
  outliars[["all"]] <- all
  
  non_na_indices <- which(!is.na(all), arr.ind = TRUE)
  non_na_values <- data.frame(
    row = rownames(all)[non_na_indices[, 1]],
    column = colnames(all)[non_na_indices[, 2]],
    value = all[non_na_indices])
  
  outliars[["all1"]] <- non_na_values
  
  #plotting all of the variables across all stations
  df<- env%>%
    select(1:8)
  data_long <- df %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")
  
  outliars[["plot"]]<- ggplot(data_long, aes(x = Variable, y = Value)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
  return(outliars)
}

# the problem is that it will remove NA rows
Site<- out("Site")
Habitat<- out("Habitat")
Season<- out("Season")
(all <- out("factors", 2))



################################# boxplots ####################################
env2<- env
env2[row.names(env2)== "C523_w", colnames(env2) == "TOM"] <- NA

ggplot(env2, aes(factors, Sal))+
  geom_boxplot(aes(color = Season, fill = Site))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

env3<- env1%>%
  filter(facts$Site == "Tis")

# Function to visualize data using boxplots
visualize_data <- function(data){
  data_long <- data %>%
    pivot_longer(cols = everything(), names_to = "Variable", values_to = "Value")
  
  ggplot(data_long, aes(x = Variable, y = Value)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Boxplot of Environmental Data")
}


visualize_data(env3)
# detect_outliars_z_score(env1)


####################### bar plots ###################
meanSAH<-list()
for(i in seq(ncol(env1))){ #mean and std of each varable
  dum<-data.frame(var = env1[,i],Season = env$Season,
                  Site = env$Site, Habitat = env$Habitat)

  SAH<- dum %>%  #grouped by Season, Site, and Habitat
    group_by(Season, Site, Habitat)%>%
    summarise(mean = mean(var, na.rm = T),
              SD = sd(var,na.rm = T),
              n = sum(!is.na(var))
              )
  meanSAH[[length(meanSAH)+1]]<- SAH
  names(meanSAH)[length(meanSAH)]<-names(env1[i])
}

plotSAH<- function(var = "Temp"){
  colind <- which(colnames(env) == var)
  ############### barplot #######################.
  v<-as.data.frame(meanSAH[[as.character(var)]])
  
  factor<- interaction(v$Site, v$Season)
  fill <- v$Habitat
  v<- data.frame(fill,factor, v[4:5])
  
  bar<- ggplot(v, aes(factor, mean,
                        fill = fill))+
    geom_col(position = "dodge",col = "black")+ #by dodge the fill factor will be side to side instead of stacked
    labs(x = "Site.seasn", y = as.character(var))+
    geom_errorbar(aes(ymin = mean-SD, ymax = mean+SD),
                  position = position_dodge2(width = 0.5, padding = 0.5)) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1))+
    scale_fill_manual(values = c("#4682B4", "#228B22", "#F4A460"))+
    theme(text = element_text(family= "Cambria", face="bold", size=12))
  return(bar)
}

(Temp<- plotSAH("Temp"))
Sal<- plotSAH("Sal")
pH<- plotSAH("pH")
# meanSAH[[3]]<- meanSAH[[3]]%>% 
  # filter(!(Season == "w"& Site == "Tis" & Habitat == "bch"))
DO<- plotSAH("DO")
sand<- plotSAH("sand")
clay<- plotSAH("clay")
silt<- plotSAH("silt")
TOM<- plotSAH("TOM")

bars<- ggarrange(Temp, Sal, pH,TOM, common.legend = T )

############ sediment relative composition bars ####.
# Aggregate the data to get mean percentages per factor


sed_agg <- env %>%
  filter(Season == "W")%>%
  group_by(Site, Habitat) %>%
  summarise(
    sand = mean(sand),
    silt = mean(silt),
    clay = mean(clay)
  )

# Melt the data to long format for ggplot2
sed_long <- sed_agg %>%
  pivot_longer(cols = c(sand, silt, clay), names_to = "particle",
               values_to = "percentage")
sed_long$particle<- factor(sed_long$particle, levels= c("clay", "silt", "sand"))

# Create the stacked bar plot
# fonts()
sed<- ggplot(sed_long, aes(x = interaction(sed_long$Habitat, sed_long$Site), 
                           y = percentage, fill = particle)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Habitat.Site",y = "fracrion(%)",fill = "sediment type")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = c("#1F77B4", "#5E3C99", "#E66100"))+
  theme(text = element_text(family= "Cambria", face="bold", size=12))
  # theme_cleveland()


  ##########################################################################.
################################ PCA ########################################

#function for plotting with manual theme 
plot<- function(){
  ggplot(pc1, aes(PC1,PC2))+
    geom_point(aes(shape = Site , color = Habitat)
               ,stroke = 1, size = 2.5)+
    theme_classic()+
    labs(x= "PC1", y= "PC2")+
    scale_color_manual( values=c( "#1F77B4", "#1b9e77", "#d95f02"))+
                        # labels = c("Mudflat", "Vegetated", "beach"))+
    theme( legend.title = element_text(size=7) ,legend.position = "right",
           legend.text = element_text(size = 6),
           axis.title = element_text(size= 14),
           axis.text = element_text(size = 12))+
    scale_shape_manual(values = c( 17, 16))+
                       # labels = c("Planted", "Natural"))+
    labs(colour="Habitat (color)")+
  labs(shape="Site (shape)")
}

plot1<- function(){
  ggplot(pc1, aes(PC1,PC2))+
    geom_point(aes(shape = Season.Site , color = Habitat)
               ,stroke = 1, size = 2.5)+
    theme_classic()+
    labs(x= "PC1", y= "PC2")+
    scale_color_manual( values=c( "#101820", "#1a7a4c", "#d95f02"))+
    theme( legend.title = element_text(size=10) ,#legend.position = "none",
           legend.text = element_text(size = 8),
           axis.title = element_text(size= 12),
           axis.text = element_text(size = 10))+
    scale_shape_manual(values = c(2, 17, 1, 16))
    }

################### two Seasons, all data, sediment copied for summer####
# 
# "#e41a1c"
# "#377eb8"
# "#22892B" '#07553B', '#CED46A'
# "#984ea3" "#101820", "#1a7a4c"
# "#ff7f00" "#343148", "#D7C49E"
# "#7070db" "#4B878B", "#921416"
# "#a65628" '#331B3F', '#ACC7B4'
# "#00cc66"
# "#E69F00", "#4995B8", "#009E73"
#  "#7570b3", "#1b9e77", "#d95f02"

env_all <- na.omit(env%>%
                     select(-DO, -sand))

pca<- prcomp(env_all%>%
               select(-names(facts)), scale = T)

psum<-summary(pca) #extracting percentages
perc1<- as.character(round(psum$importance[2,1],3)*100)
perc2<- as.character(round(psum$importance[2,2],3)*100)

pc<- as.data.frame(pca[["x"]])[1:2]
pc1<- cbind(pc, env_all[names(env_all) %in% names(facts)])
load<- data.frame(var = rownames(pca$rotation), pca$rotation) 

(p4<- ggplot(pc1, aes(PC1,PC2))+
  geom_point(aes(shape = Season.Site , color = Habitat),
             stroke = 1, size = 2.2, alpha = 0.8)+
  theme_classic()+
  labs(x= paste("PC1","(",perc1,"%",")"),
       y= paste("PC2","(",perc2,"%",")"))+
  scale_color_manual( values=c( "#101820", "#4B878B", "#d95f02"))+
  theme( legend.title = element_text(size=10) ,#legend.position = "bottom",
         legend.text = element_text(size = 8),
         axis.title = element_text(size= 12),
         axis.text = element_text(size = 10))+
  scale_shape_manual(values = c(2, 17, 1, 16))+ 
  geom_segment(data = load, aes(x = 0, y = 0, xend = (PC1*3),yend = (PC2*3)),
               arrow = arrow(length = unit(1/2, "picas")),
               alpha = 1, color = "black")+
  annotate("text", x = (load$PC1*4),alpha = 1, y = (load$PC2*3), label = load$var))

# (p1<- autoplot(pca, data = env_all , color = "Season", loadings = T,
#                loadings.colour = 'blue',
#                loadings.label = TRUE, loadings.label.size = 4))
# (p2<- autoplot(pca, data = env_all , color = "Site", loadings = T,
#                loadings.colour = 'blue',
#                loadings.label = TRUE, loadings.label.size = 4))
# (p3<- autoplot(pca, data = env_all , color = "Habitat", loadings = T,
#                loadings.colour = 'blue',
#                loadings.label = TRUE, loadings.label.size = 4)) 

  # 3 factors in one plot
# plot1() +stat_ellipse(aes(shape = Habitat))
# stat_ellipse(aes(shape =Site),type = "t", level = 0.8, linetype = 1, )
# (p4<- plot1() + 
#     geom_segment(data = load, aes(x = 0, y = 0, xend = (PC1*3),yend = (PC2*3)),
#                  arrow = arrow(length = unit(1/2, "picas")), color = "black")+
#     annotate("text", x = (load$PC1*4), y = (load$PC2*3), label = load$var))

  # two factors 

  
    
pca_all<- p4
# pca_all4<- ggarrange(p1, p2, p3, p4)
########################## for data of winter with sediment size ####

env1 <- select(env, -DO, -sand ) #dropping DO since it has not been properly measured
env_w<- na.omit(env1[env1$Season == "W",])

pca<- prcomp(env_w%>%
               select(-names(facts)), scale = T)
pca
summary(pca)
psum<-summary(pca) #extracting percentages
perc1<- as.character(round(psum$importance[2,1],3)*100)
perc2<- as.character(round(psum$importance[2,2],3)*100)

pc<- as.data.frame(pca[["x"]])[1:2]
pc1<- cbind(pc, env_w[names(env_w) %in% names(facts)])

# (p2<- autoplot(pca, data = env_w , color = "Site", loadings = T,
#                loadings.colour = 'blue',
#                loadings.label = TRUE, loadings.label.size = 4))
# (p3<- autoplot(pca, data = env_w , color = "Habitat", loadings = T,
#                loadings.colour = 'blue',
#                loadings.label = TRUE, loadings.label.size = 4))

load<- data.frame(var = rownames(pca$rotation), pca$rotation) 

(p4<- plot() +
    geom_segment(data = load, aes(x = 0, y = 0, xend = (PC1*2.5),yend = (PC2*2.5)),
                 arrow = arrow(length = unit(1/2, "picas")), color = "gray30")+
    stat_ellipse(aes(shape =Site), alpha = 0.3 )+
    labs(x= paste("PC1","(",perc1,"%",")"),
         y= paste("PC2","(",perc2,"%",")"))+
    scale_color_manual( values=c( "#264653", "#2a9d8f", "#e9c46a"))+
                        # labels = c("Mudflat", "Vegetated", "beach"))+
    theme( legend.title = element_text(size=12) ,#legend.position = 'bottom',
           legend.text = element_text(size = 10),
           axis.title = element_text(size= 12),
           axis.text = element_text(size = 10))+
    annotate("text", x = (load$PC1*c(2.3,3,3,3,3,3.5)), y = (load$PC2*c(3.5,2.8,2.5,3,2.5,2.3)),
             label = load$var))

# pca_win_4<- ggarrange(p2, p3, p4)
pca_win<- p4

##################### for data of all Seasons without sediment size ####
env1 <- na.omit(env%>%
                  select(-clay, -sand, -silt, -DO))

pca<- prcomp(env1%>%
               select(-names(facts)), scale = T)
pc<- as.data.frame(pca[["x"]])[1:2]
pc1<- cbind(pc, env1[names(env1) %in% names(facts)])

(p1<- autoplot(pca, data = env1 , color = "Season", loadings = T,
              loadings.colour = 'blue',
              loadings.label = TRUE, loadings.label.size = 4))
(p2<- autoplot(pca, data = env1 , color = "Site", loadings = T,
               loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size = 4))
(p3<- autoplot(pca, data = env1 , color = "Habitat", loadings = T,
               loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size = 4))
(p23<- autoplot(pca, data = env1 , color = "Habitat",shape = "Site", loadings = T,
               loadings.colour = 'blue',
               loadings.label = TRUE, loadings.label.size = 4))

load<- data.frame(var = rownames(pca$rotation), pca$rotation) 

(p4<- plot1()+
    geom_segment(data = load, aes(x = 0, y = 0, xend = (PC1*3),yend = (PC2*3)),
                 arrow = arrow(length = unit(1/2, "picas")), color = "black")+
    annotate("text", x = (load$PC1*4), y = (load$PC2*3), label = load$var))


sed_drop4<- ggarrange(p1, p2, p3, p4)
sed_drop<- p4


############### output ##########################
#PCA all
# ggsave("figs/PCA.all.png", plot =  pca_all , width = 8, height = 5)
# ggsave("figs/PCA.all(4).png", plot = pca_all4 , width = 8, height = 5)

#PCA sediment droped
# ggsave("figs/PCA.sediment droped.png", plot = sed_drop , width = 8, height = 5)
# ggsave("figs/PCA.sediment droped(4).png", plot = sed_drop4 , width = 8, height = 5)

#PCA winer
 ggsave("figs/PCA.all.winter.png", plot = pca_win , width = 7, height = 5)
# ggsave("figs/PCA.all.winter(4).png", plot = pca_win_4 , width = 8, height = 5)

#bar plots
# ggsave("figs/env_bars.png", plot = bars, width = 9, height = 10)
# ggsave("figs/sed_bar.png", plot = sed, width = 6, height = 5)

#table of means
write.csv(av_sd,"output/env.ave.csv")


