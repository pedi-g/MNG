#Pedram Ghahramani
#5/12/2025
#ABC curve (abundance - biomass)

# Load required libraries
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(forams)

###########################  Prepare data  ###########################################
#input
source("R-scripts/factors.R")
ab<-  read.csv("raw data/ab.csv", header = T, row.names = 1)%>% # abundance of all species
  filter(habitat != "bch")
ab<- ab%>% select(-station)

mass<-  read.csv("raw data/mass.csv", header = T, row.names = 1)# biomass of all species
mass<- cbind(mass,factors[2:4])%>%
  filter(factors$habitat!= "bch")
# Reshape the abundance data to long format
abundance <- ab %>%
  pivot_longer(cols = -c(area, season, habitat), names_to = "species", values_to = "abundance")

# Reshape the biomass data to long format
biomass<- mass %>%
  pivot_longer(cols = -c(area, season, habitat), names_to = "species", values_to = "biomass")

############## all in one function #########################

#a function that selects data by factors and computes cumulative sum of abundance or biomass
filt <- function(s= NULL, a = NULL,hab = NULL){
  ab<-abundance
  mass<- biomass
  #filtering accordinf to the input
  if (!is.null(s)) {
    ab <- ab %>% filter(season == s)
    mass <- mass %>% filter(season == s)
  }
  if (!is.null(a)) {
    ab <- ab %>% filter(area == a)
    mass <- mass %>% filter(area == a)
  }
  
  if (!is.null(hab)) {
    ab <- ab %>% filter(habitat == hab)
    mass <- mass %>% filter(habitat == hab)
  }
  
  # sum of abundance and biomass of each species
  ab<- ab%>%
    group_by(species) %>%
    summarize(absum = sum(abundance)) %>%
    arrange(desc(absum))%>%
    filter(absum > 0)
  mass<- mass%>%
    group_by( species) %>%
    summarize(msum = sum(biomass))%>%
    arrange(desc(msum))%>%
    filter(msum > 0)
  
  #saving for another method of ploting ( forams::abs() function)
  y<- as.data.frame(cbind(N = ab$absum, Biomass = mass$msum))
  
  #cumulative relative sum of each value
  ab<- ab%>%
    reframe( cum = cumsum(absum/sum(absum)))
  mass<- mass%>%
    reframe( cum = cumsum(msum/sum(msum)))
  
  #saving it as an input for the ggplot bellow 
  x<- as.data.frame(cbind(ab.cum = ab$cum, mass.cum = mass$cum))
  
  return(list(y,ggplot(x, aes(x = 1:nrow(x), y = ab.cum, color = "Abundance")) +
                geom_line() +
                geom_line(aes(y = mass.cum, color = "Biomass")) +
                labs(x = "Species rank", y = "Cumulative percentage", color = "") +
                scale_x_log10()+
                theme_classic()))
}


######## filtering by factors  ################

# all in one
all<- filt()

# separating by area
t<-filt(a = "Tis")
g<- filt(a = "Gw")

#separated by season
s<-  filt(s = "s")
w<-  filt(s = "w")

#separated by habitat
b<- filt(hab = "bch")
m<- filt(hab = "mng")
c<- filt(hab = "crk")

#separated by season and area
st<- filt("s", "Tis")
wt<- filt("w", "Tis")
sg<- filt("s", "Gw")
wg<- filt("w", "Gw")

#separated by season and habitat
sc<- filt(s = "s", hab = 'crk')
sm<- filt(s = "s", hab = 'mng')
sb<- filt(s = "s", hab = 'bch')
wc<- filt(s = "w", hab = 'crk')
wm<- filt(s = "w", hab = 'mng')
wb<- filt(s = "w", hab = 'bch')

#separated by area and habitat
tb<- filt(a = "Tis", hab = "bch")
tm<- filt(a = "Tis", hab = "mng")
tc<- filt(a = "Tis", hab = "crk")
gb<- filt(a = "Gw", hab = "bch")
gm<- filt(a = "Gw", hab = "mng")
gc<- filt(a = "Gw", hab = "crk")

#separated by season, area, and habitat
stb<- filt("s", "Tis", "bch")
stm<- filt("s", "Tis", "mng")
stc<- filt("s", "Tis", "crk")
sgb<- filt("s", "Gw", "bch")
sgm<- filt("s", "Gw", "mng")
sgc<- filt("s", "Gw", "crk")
wtb<- filt("w", "Tis", "bch")
wtm<- filt("w", "Tis", "mng")
wtc<- filt("w", "Tis", "crk")
wgb<- filt("w", "Gw", "bch")
wgm<- filt("w", "Gw", "mng")
wgc<- filt("w", "Gw", "crk")

# usind forams::abc() to plot function

#no factorization
png("figs/ABC_all.png",  width=2300, height=2000, res=300)
par(mfrow = c(1,1))
plot(abc(all[[1]]))
dev.off()

# area
png("figs/ABC_area.png",  width=3000, height=2000, res=300)
par(mfrow = c(1,2))  
plot(abc(t[[1]]))
mtext("Tis", side = 3, adj = 0)
plot(abc(g[[1]]))
mtext("Gwatr", side = 3, adj = 0)
dev.off()

#season
png("figs/ABC_season.png",  width=3000, height=2000, res=300)
par(mfrow = c(1,2))
plot(abc(s[[1]]))
mtext("summer", side = 3, adj = 0)
plot(abc(w[[1]]))
mtext("winter", side = 3, adj = 0)
dev.off()

#habitat
png("figs/ABC_habitat.png",  width=2000, height=1000, res=300)
par(mfrow = c(1,2))
# plot(abc(b[[1]]))
# mtext("beach", side = 3, adj = 0)
plot(abc(m[[1]]))
mtext("mangrove", side = 3, adj = 0)
plot(abc(c[[1]]))
mtext("creek", side = 3, adj = 0)
dev.off()

#season and area
png("figs/ABC_season-area.png",  width=2300, height=2000, res=300)
par(mfrow = c(2,2))
plot(abc(wt[[1]]))
mtext("wt", side = 3, adj = 0)
plot(abc(wg[[1]]))
mtext("wg", side = 3, adj = 0)
plot(abc(st[[1]]))
mtext("st", side = 3, adj = 0)
plot(abc(sg[[1]]))
mtext("sg", side = 3, adj = 0)
dev.off()

#season and habitat
png("figs/ABC_season-habtiat.png",  width=3000, height=2000, res=300)
par(mfrow = c(2,2))
plot(abc(wc[[1]]))
mtext("wc", side = 3, adj = 0)
plot(abc(wm[[1]]))
mtext("wm", side = 3, adj = 0)
# plot(abc(wb[[1]]))
# mtext("wb", side = 3, adj = 0)
plot(abc(sc[[1]]))
mtext("sc", side = 3, adj = 0)
plot(abc(sm[[1]]))
mtext("sm", side = 3, adj = 0)
# plot(abc(sb[[1]]))
# mtext("sb", side = 3, adj = 0)
dev.off()

# area and habtiat
png("figs/ABC_area-habtiat.png",  width=3000, height=2000, res=300)
par(mfrow = c(2,2))
plot(abc(tc[[1]]))
mtext("tc", side = 3, adj = 0)
plot(abc(tm[[1]]))
mtext("tm", side = 3, adj = 0)
# plot(abc(tb[[1]]))
# mtext("tb", side = 3, adj = 0)
plot(abc(gc[[1]]))
mtext("gc", side = 3, adj = 0)
plot(abc(gm[[1]]))
mtext("gm", side = 3, adj = 0)
# plot(abc(gb[[1]]))
# mtext("gb", side = 3, adj = 0)
dev.off()
# season. area. habitat
png("figs/ABC_season-area-hab.png",  width=3000, height=2000, res=200)
par(mfrow = c(2,4))
plot(abc(wtc[[1]])); mtext("wtc", side = 3, adj = 0)
plot(abc(wgc[[1]])); mtext("wgc", side = 3, adj = 0)
plot(abc(stc[[1]])); mtext("stc", side = 3, adj = 0)
plot(abc(sgc[[1]])); mtext("sgc", side = 3, adj = 0)
plot(abc(wtm[[1]])); mtext("wtm", side = 3, adj = 0)
plot(abc(wgm[[1]])); mtext("wgm", side = 3, adj = 0)
plot(abc(stm[[1]])); mtext("stm", side = 3, adj = 0)
plot(abc(sgm[[1]])); mtext("sgm", side = 3, adj = 0)
# plot(abc(wtb[[1]])); mtext("wtb", side = 3, adj = 0)
# plot(abc(wgb[[1]])); mtext("wgb", side = 3, adj = 0)
# plot(abc(stb[[1]])); mtext("stb", side = 3, adj = 0)
# plot(abc(sgb[[1]])); mtext("sgb", side = 3, adj = 0)
dev.off()

#arrangment##########
# area<- ggarrange(t[[2]], g[[2]],  common.legend = T)
# season<- ggarrange(s[[2]], w[[2]],  common.legend = T)
# habitat<- ggarrange(b[[2]], m[[2]], c[[2]], ncol = 3, common.legend = T)
# seas.area<- ggarrange(wt[[2]], wg[[2]], st[[2]], sg[[2]],
#                       common.legend = T)
# seas.hab<- ggarrange(wc[[2]], wm[[2]], wb[[2]], sc[[2]],
#                      sm[[2]], sb[[2]],  common.legend = T)
# area.habitat<- ggarrange(tc[[2]], tm[[2]], tb[[2]], gc[[2]],
#                          gm[[2]], gb[[2]],  common.legend = T)
# facts<- ggarrange(wtc[[2]], wtm[[2]], wtb[[2]], stc[[2]],
#                 stm[[2]], stb[[2]], wgc[[2]], wgm[[2]],
#                 wgb[[2]], sgc[[2]], sgm[[2]],sgb[[2]],
#                 common.legend = T, legend = "bottom", ncol = 3, nrow = 4)

# ggsave("figs/ABC(3factors).png", plot = facts , width = 10, height = 10)
# ggsave("figs/ABC(area).png", plot = area , width = 12, height = 8)
# ggsave("figs/ABC(season).png", plot = season , width = 12, height = 8)
# ggsave("figs/ABC(habitat).png", plot = habitat , width = 15, height = 15)
# ggsave("figs/ABC(area.season).png", plot = seas.area , width = 10, height = 9)
# ggsave("figs/ABC(season.habitat).png", plot = seas.hab , width = 12, height = 8)
# ggsave("figs/ABC(area.habitat).png", plot = area.habitat , width = 12, height = 10)
# ggsave("figs/ABC(all).png", plot = all[[2]] , width = 8, height = 8)
