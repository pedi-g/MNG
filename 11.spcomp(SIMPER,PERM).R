#species and trait composition analysis
#SIMPER and PERMANOVA
#Pedram.ghrm@gmail.com
#8-3-2024

#libraries
library(vegan)
library(tidyverse)
library(pairwiseAdonis)

# inputdata
source('R-scripts/factors.R') #factors
#gower distance for community-trait composition
dist.t <- read.csv("output/dist.comXtrait(gowdis).csv", row.names = 1)
#community-trait matrix
trait<- read.csv('output/comXtrait(LoGab).csv', row.names = 1)%>%
  select(startS_with(c('Wei','FD','Mv','BT','T','RT','LS')))
#community-species matrix
ab_raw<- read.csv("raw data/ab.csv", row.names = 1)
ab<-log( 1 + ab_raw[,-((ncol(ab_raw)-3):ncol(ab_raw))] )
ab_raw<- cbind(ab, factors)

####### station selection #########

######### seasonal filter
# factors<- factors%>% filter(transect != "1", transect != "2", transect!= "7", transect!= "8")

######## selecting habitats
factors<- factors%>% filter(habitat != "bch")

 ab_raw<- ab_raw[row.names(ab_raw) %in% row.names(factors),]
ab<- ab_raw%>%
  select(-names(factors))
dist.t<- as.dist(dist.t[rownames(dist.t) %in% rownames(factors),
                    colnames(dist.t) %in% rownames(factors)])
trait<- trait[rownames(trait) %in% rownames(factors),]



simp_sp<- simper(ab, factors$season.area.habitat)
perm_sp<- adonis2(ab ~ season * area * habitat, data = ab_raw,
           method = "bray", permutations = 9999)
pairwise_sp <- pairwise.adonis(ab, ab_raw$season.area.habitat, p.adjust.m = "fdr", perm = 9999)

simp_trait<- simper(trait, factors$season.area.habitat)
perm_trait<- adonis2(dist.t ~ season * area * habitat, data = factors, permutations = 9999)
pairwise_t <- pairwise.adonis(dist.t, ab_raw$season.area.habitat, p.adjust.m = "fdr", perm = 9999)


simp_sp$Mud.Tis.W_Mud.Tis.S$cusum
simp_sp$Mng.Tis.W_Mng.Tis.S$cusum
# simp_sp$bch.Tis.W_bch.Tis.S$cusum
simp_sp$Mud.Gw.W_Mud.Gw.S$cusum
simp_sp$Mng.Gw.W_Mng.Gw.S$cusum
# simp_sp$bch.Gw.W_bch.Gw.S$cusum
simp_sp$Mud.Tis.W_Mng.Tis.W$cusum
# simp_sp$Mud.Tis.W_bch.Tis.W$cusum
# simp_sp$Mng.Tis.W_bch.Tis.W$cusum
simp_sp$Mud.Tis.S_Mng.Tis.S$cusum
# simp_sp$Mud.Tis.S_bch.Tis.S$cusum
# simp_sp$Mng.Tis.S_bch.Tis.S$cusum
simp_sp$Mud.Gw.W_Mng.Gw.W$cusum
# simp_sp$Mud.Gw.W_bch.Gw.W$cusum
# simp_sp$Mng.Gw.W_bch.Gw.W$cusum
simp_sp$Mud.Gw.S_Mng.Gw.S$cusum
# simp_sp$Mud.Gw.S_bch.Gw.S$cusum
# simp_sp$Mng.Gw.S_bch.Gw.S$cusum
simp_sp$Mud.Tis.W_Mud.Gw.W$cusum
simp_sp$Mng.Tis.W_Mng.Gw.W$cusum
# simp_sp$bch.Tis.W_bch.Gw.W$cusum
simp_sp$Mud.Tis.S_Mud.Gw.S$cusum
simp_sp$Mng.Tis.S_Mng.Gw.S$cusum
# simp_sp$bch.Tis.S_bch.Gw.S$cusum

simp_trait$Mud.Tis.W_Mud.Tis.S$cusum
simp_trait$Mng.Tis.W_Mng.Tis.S$cusum
# simp_trait$bch.Tis.W_bch.Tis.S$cusum
simp_trait$Mud.Gw.W_Mud.Gw.S$cusum
simp_trait$Mng.Gw.W_Mng.Gw.S$cusum
# simp_trait$bch.Gw.W_bch.Gw.S$cusum
simp_trait$Mud.Tis.W_Mng.Tis.W$cusum
# simp_trait$Mud.Tis.W_bch.Tis.W$cusum
# simp_trait$Mng.Tis.W_bch.Tis.W$cusum
simp_trait$Mud.Tis.S_Mng.Tis.S$cusum
# simp_trait$Mud.Tis.S_bch.Tis.S$cusum
# simp_trait$Mng.Tis.S_bch.Tis.S$cusum
simp_trait$Mud.Gw.W_Mng.Gw.W$cusum
# simp_trait$Mud.Gw.W_bch.Gw.W$cusum
# simp_trait$Mng.Gw.W_bch.Gw.W$cusum
simp_trait$Mud.Gw.S_Mng.Gw.S$cusum
# simp_trait$Mud.Gw.S_bch.Gw.S$cusum
# simp_trait$Mng.Gw.S_bch.Gw.S$cusum
simp_trait$Mud.Tis.W_Mud.Gw.W$cusum
simp_trait$Mng.Tis.W_Mng.Gw.W$cusum
# simp_trait$bch.Tis.W_bch.Gw.W$cusum
simp_trait$Mud.Tis.S_Mud.Gw.S$cusum
simp_trait$Mng.Tis.S_Mng.Gw.S$cusum
# simp_trait$bch.Tis.S_bch.Gw.S$cusum


testsL1<- list()
testsL1$perm_sp<- perm_sp
testsL1$perma_t<- perm_trait

pair_ord<- c(
  'Mud.Tis.W vs Mud.Tis.S','Mng.Tis.W vs Mng.Tis.S','bch.Tis.W vs bch.Tis.S', #seasonal Tis
  'Mud.Gw.W vs Mud.Gw.S', 'Mng.Gw.W vs Mng.Gw.S','bch.Gw.W vs bch.Gw.S',   #seasonal G
  'Mud.Tis.W vs Mng.Tis.W', 'Mud.Tis.W vs bch.Tis.W','Mng.Tis.W vs bch.Tis.W', #habitat T win
  'Mud.Tis.S vs Mng.Tis.S', 'Mud.Tis.S vs bch.Tis.S','Mng.Tis.S vs bch.Tis.S', #habitat T sum
  'Mud.Gw.W vs Mng.Gw.W','Mud.Gw.W vs bch.Gw.W', 'Mng.Gw.W vs bch.Gw.W', # hab G win
  'Mud.Gw.S vs Mng.Gw.S', 'Mud.Gw.S vs bch.Gw.S', 'Mng.Gw.S vs bch.Gw.S', #hab G sum
  'Mud.Tis.W vs Mud.Gw.W', 'Mng.Tis.W vs Mng.Gw.W', 'bch.Tis.W vs bch.Gw.W', # area win
  'Mud.Tis.S vs Mud.Gw.S', 'Mng.Tis.S vs Mng.Gw.S','bch.Tis.S vs bch.Gw.S' #area sum
)

  filtered_rows <- pairwise_sp[
    pairwise_sp[,1] %in% pair_ord,]
  testsL1$pairwise_sp <- filtered_rows[match(pair_ord, filtered_rows[, 1]), ]

filtered_rows <- pairwise_t[
  pairwise_t[,1] %in% pair_ord,]
testsL1$pairwise_t <- filtered_rows[match(pair_ord, filtered_rows[, 1]), ]

testsL1$pairwise_sp <- na.omit(testsL1$pairwise_sp)
testsL1$pairwise_t <- na.omit(testsL1$pairwise_t)
testsL1
