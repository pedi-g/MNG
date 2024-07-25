#Ped_gh
#pedram.ghrm@gmail.com
#env data normality and test of significancy
#3-6-2024

#libraries----------
library(tidyverse)
library(vegan)
library(psych)
library(pairwiseAdonis)
#input-----------
env <- read.csv("raw data/env1.csv", row.names = 1)
#defining factor variables--------
env$station <- factor(env$station)
env$season <- factor(env$season)
env$habitat<- factor(env$habitat, levels= c("crk", "mng", "bch"))
env$area <- factor (env$area, levels = c("Tis", "Gw"))
env$season.area <- factor( paste(env$season, env$area, sep = "_"), 
                           levels = c("w_Tis", "w_Gw", "s_Tis", "s_Gw"))
env$factors <- factor( paste(env$season, env$area, env$habitat, sep = "."), 
                       levels = c("w.Tis.crk", "w.Tis.mng", "w.Tis.bch",
                                  "w.Gw.crk", "w.Gw.mng", "w.Gw.bch",
                                  "s.Tis.crk", "s.Tis.mng", "s.Tis.bch",
                                  "s.Gw.crk", "s.Gw.mng", "s.Gw.bch"))

env<- env%>%
  filter(habitat != "bch") #removing beach habitat
# removing outliers
 #env["C423_w",	"Sal" ]<- NA
# env["C221_w",	"DO"	]<- NA
# env['C222_w',	'DO'  ]<- NA
# env['C711_w',	'DO'  ]<- NA
# env['C511_w',	'clay']	<- NA
# env['C511_s',	'clay']	<- NA
 env['C523_w',	'TOM'	]<- NA


env1<- env%>%
  select(Temp, Sal, DO, pH, clay, silt, sand, TOM)


# if normal -> use ANOVA
  # if significant -> use post hoc test -> Tukey's HSD
# if not normal -> non parametric tests ->
  # if multiple groups -> Kruskal-wallis
    # if significant -> post hoc -> Dunn's test
  # if two indipendent grups -> Mann-Whitney U (Wilcoxon Rank-Sum Test)
  # if paired samples (e.g., before and after measurements within the same habitat)
   # -> use Wilcoxon Signed-Rank Test
  
#PERMANOVA is ANOVA for multivariate, non-normal data. MANOVA is the same for normal data


############ Normality ################

all<- list()
for(var in 1:ncol(env1)){ # the outer loop for chosing each env variable
  l<- data.frame()
  
  for(group in unique(env$factors)){ # the inner loop chooses each of the groups 

    data<- env[env$factors == group,var]
    a <- shapiro.test(data)
    l[nrow(l)+1,1]<- a[[1]]
    l[nrow(l),2]<- a[[2]]
    rownames(l)[nrow(l)] <- as.character(group) 
    
    #     l[[length(l)+1]]<- as.data.frame(w = a[[1]], p = a[[2]])
    # names(l)[length(l)]<- names(group)
  }
  
  names(l) <- c("w", "p-value")
  all[[length(all)+1]]<- l
  names(all)[length(all)]<- names(env1[var])
}

    
############## ANOVA and Tukey's ###################

anov<- list()
tuk<- list()
for (i in names(env1)){
  var<- env[c(i, 'season', 'area', 'habitat')]
  names(var)[1]<- "a"
  
  an<- aov(a ~ season * area * habitat, data = var)
  t<- TukeyHSD(an)
  
  anov[[length(anov)+1]]<- an
  tuk[[length(tuk)+1]]<- t
  
  names(anov)[length(anov)]<- names(env[i])
  names(tuk)[length(tuk)]<- names(env[i])
}

summary(anov$Sal)
tuk$Sal

############ post hoc for normal distribution##############
# 1. tukey's Honest Significant Difference:
  # Ideal for comparing all possible pairs of groups.

# 2.Bonferroni Correlation
  # A more conservative approach that adjusts p-values to control the family-wise error rate.
pairwise.t.test(env$Temp, env$factors, p.adjust.method = "bonferroni")
# print(pairwise_result)

# Scheffé's Method:
  # Useful when you have unequal sample sizes or are making many comparisons.
# library(DescTools)
# scheffe_result <- ScheffeTest(aov_result)
# print(scheffe_result)

####################### permanova ##############################

perm<- list()
for (i in names(env1)){

  var<- env[c(i, 'season', 'area', 'habitat')]
  names(var)[1]<- "a"
  var<- na.omit(var)
  a<- as.matrix(var[1])
  # a <- as.matrix(var[1])
  per <- adonis2(a ~ season * area * habitat, data = var,
                 method = "euclidean", permutations = 9999)
  
  perm[[length(perm)+1]]<- per
  names(perm)[length(perm)]<- names(env[i])
  
}

# PERMANOVA for sediments
sed<-na.omit(
  env[c('sand', 'silt', 'clay', 'season', 'area', 'habitat', 'factors')])
perm[[length(perm)+1]]<- adonis2(sed[1:3] ~ season * area * habitat, data = sed,
                 method = "euclidean", permutations = 9999, na.rm = T)
names(perm)[length(perm)]<- "sed_size"

  ############# pair wise test ###################

## solution two (more data in the function)
pairwise<- list()
for (i in names(env1)){
var<- na.omit(env[c(i,"factors")])
names(var)[1]<- "a"
dis<- dist(var$a)
pair <- pairwise.adonis(dis, var$factors, p.adjust.m='holm', sim.method = "euclidean")

pairwise[[length(pairwise)+1]]<- pair
names(pairwise)[length(pairwise)]<- names(env[i])
}

# pairwise ofr sediment 
dis_sed<- dist(sed[1:3])
pairwise[[length(pairwise)+1]]<- pairwise.adonis(dis_sed, sed$factors,
                                                 sim.method = "euclidean",
                                                 p.adjust.m = 'fdr')
names(pairwise)[length(pairwise)]<- "sed_size"

perm$Sal
pairwise$Sal

#deleting excesive variables from the environment
variables <- ls()
variables<- variables[!variables %in% c('all', 'tuk', 'anov', 'pairwise', 'perm', 'env', 'env1')]
rm(list = variables)
rm(list = "variables")
# ...

## alternative solution 
# library(RVAideMemoire)
# pairwise_result_season1 <- pairwise.perm.manova(dis, temp$factors, env$factors, nperm = 999)

########### correlation #################
 ### spearman correlation for env data without trasnforamtion
 #correlations####
 # png("figs/corr.env.pairs.png", width = 900, height = 900, pointsize = 20)
 # pairs.panels(env1, gap=0,pch=20,method="spearman", 
 #              ellipses = F, lm= T, density = F, hist.col ="#236573", 
 #              jiggle = F,rug=F, smoother=F,stars=T,ci=F,alpha=.05)
 # # dev.off()
 # ## with transformation
 # 
 # env_log<- log(env1+1)
 # pairs.panels(env_log, gap=0,pch=20,method="spearman", 
 #              ellipses = F, lm= T, density = F, hist.col ="#236573", 
 #              jiggle = F,rug=F, smoother=F,stars=T,ci=F,alpha=.05)
 # 
 