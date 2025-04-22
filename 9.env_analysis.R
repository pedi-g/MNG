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
env0 <- read.csv("raw data/env1.csv", row.names = 1)

#defining factor variables--------
env0$station <- factor(env0$station)
env0$transect <- factor(env0$transect)
env0$season <- factor(env0$season, levels = c("w", "s"))
env0$habitat<- factor(env0$habitat, levels= c("crk", "mng", "bch"))
env0$area <- factor (env0$area, levels = c("Tis", "Gw"))
env0$season.area <- interaction(env0$area, env0$season)
env0$season.habitat <- interaction(env0$season, env0$habitat)
env0$area.habitat <- interaction(env0$area, env0$habitat)
env0$season.area.habitat <- interaction(env0$habitat, env0$area, env0$season)


#################### station filters #######################

# ########## seasonal filter 
# env0<- env0%>%
#   filter(transect != "1", transect != "2", transect!= "7", transect!= "8")

######## filter bch out
# env0<- env0%>%
#    filter(habitat != "bch") #removing beach habitat


### removing outliers
 #env["C423_w",	"Sal" ]<- NA
# env["C221_w",	"DO"	]<- NA
# env['C222_w',	'DO'  ]<- NA
# env['C711_w',	'DO'  ]<- NA
# env['C511_w',	'clay']	<- NA
# env['C511_s',	'clay']	<- NA
 # env['C523_w',	'TOM'	]<- NA


sed<-env0[#env0$season.area.habitat!= 's.Tis.bch',
    c('sand', 'silt', 'clay', 'season', 'area', 'habitat', 'season.area.habitat')]
sed<- sed[sed$season == 'w',]
TOM<- env0[c("TOM", 'season', 'area', 'habitat', 'season.area.habitat')]

env<- na.omit(
  env0%>%
  select(-clay, -silt, -sand, -TOM)
  )

# env<- env0 

env1<- env[which(names(env) %in% c("Temp", "Sal", "DO", "pH"))]

### if normal -> use ANOVA
  ## if significant -> use post hoc test -> Tukey's HSD
###if not normal -> non parametric tests ->
  ## if multiple groups -> Kruskal-wallis
    # if significant -> post hoc -> Dunn's test
  ## if two indipendent grups -> Mann-Whitney U (Wilcoxon Rank-Sum Test)
  ## if paired samples (e.g., before and after measurements within the same habitat)
   # -> use Wilcoxon Signed-Rank Test
  
#PERMANOVA is ANOVA for multivariate, non-normal data. MANOVA is the same for normal data


############ Normality ################
env2<- env0[which(sapply(env0, is.numeric))]

# env2<-
#   sqrt(env2)
  # log(1+ env2)

env2<- cbind(env2, fact = env0$season.area.habitat)
 ### delete the NAs or the 0 NAs replaces by zero for shapiro to work
all<- list()
for(var in 1:(ncol(env2)-1)){ # the outer loop for choosing each env variable
  l<- data.frame()
 df<- na.omit(
cbind(env0[var], fact = env0$season.area.habitat))
  for(group in unique(df$fact)){ # the inner loop chooses each of the groups

    data<- na.omit(env2[env2$fact == group,var])
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
all
    
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

summary(anov$Temp)
tuk$Temp

############ post hoc for normal distribution##############
# 1. tukey's Honest Significant Difference:
  # Ideal for comparing all possible pairs of groups.

# 2.Bonferroni Correlation
  # A more conservative approach that adjusts p-values to control the family-wise error rate.
# pairwise.t.test(env$Temp, env$season.area.habitat, p.adjust.method = "bonferroni")
# print(pairwise_result)

# Scheffé's Method:
  # Useful when you have unequal sample sizes or are making many comparisons.
# library(DescTools)
# scheffe_result <- ScheffeTest(aov_result)
# print(scheffe_result)

####################### permanova ##############################

perm<- list()

perm$Temp<- adonis2(env$Temp ~ season * area * habitat, data = env,
                   method = "euclidean", permutations = 9999)
perm$Sal<- adonis2(env$Sal ~ season * area * habitat, data = env,
                   method = "euclidean", permutations = 9999)
perm$pH<- adonis2(env$pH ~ season * area * habitat, data = env,
                   method = "euclidean", permutations = 9999)
# PERMANOVA for sediments

perm$Sed_size<- adonis2(sed[1:3] ~ season * area * habitat,
                        data = sed,
                 method = "euclidean", permutations = 9999, na.rm = T)

perm$TOM<- adonis2(TOM[1] ~ season * area * habitat,data = TOM,
                   method = "euclidean", permutations = 9999, na.rm = T)
  
perm
############# pair wise test ###################

## solution two (more data in the function)
# use one of these as adjustment method for p value: "fdr", or "BY"
pairwise<- list()
pairwise$Temp <- pairwise.adonis(dist(env$Temp), env$season.area.habitat, perm = 9999,
                                p.adjust.m='fdr', sim.method = "euclidean")
pairwise$Sal <- pairwise.adonis(dist(env$Sal), env$season.area.habitat, perm = 9999,
                                p.adjust.m='fdr', sim.method = "euclidean")
pairwise$pH <-  pairwise.adonis(dist(env$pH), env$season.area.habitat, perm = 9999,
                                p.adjust.m='fdr', sim.method = "euclidean")


# pairwise for sediment 

pairwise$sed_size<- pairwise.adonis(dist(sed[1:3]), sed$season.area.habitat, perm = 9999,
                                                 sim.method = "euclidean",
                                                 p.adjust.m = 'fdr')
pairwise$TOM <-  pairwise.adonis(dist(TOM$TOM), TOM$season.area.habitat, perm = 9999,
                                 p.adjust.m='fdr', sim.method = "euclidean")

pairwise$Temp
pairwise$Sal
pairwise$pH  
pairwise$sed_size
pairwise$TOM
#deleting excesive variables from the environment ----
# variables <- ls()
# variables<- variables[!variables %in% c('all', 'tuk', 'anov', 'pairwise', 'perm', 'env', 'env0')]
# rm(list = variables)
# rm(list = "variables")


## alternative solution 
# library(RVAideMemoire)
# pairwise_result_season1 <- pairwise.perm.manova(dis, temp$season.area.habitat, env$season.area.habitat, nperm = 999)








############## kruskal-Wallis##############

# krus<- kruskal.test(env$Sal ~ season * area * habitat, data = env) #Kruskal-Wallis
# wilcox.test() # Mann-Whitney U

########### correlation #################
 ### spearman correlation for env data without trasnforamtion
 #correlations

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



############ selecting desired pairs and in order#################
tests<- list()
tests$anova<- anov
tests$permanova<- perm

pair_ord<- c(
  'crk.Tis.w vs crk.Tis.s','mng.Tis.w vs mng.Tis.s','bch.Tis.w vs bch.Tis.s', #seasonal Tis
  'crk.Gw.w vs crk.Gw.s', 'mng.Gw.w vs mng.Gw.s','bch.Gw.w vs bch.Gw.s',   #seasonal G
  'crk.Tis.w vs mng.Tis.w', 'crk.Tis.w vs bch.Tis.w','mng.Tis.w vs bch.Tis.w', #habitat T win
  'crk.Tis.s vs mng.Tis.s', 'crk.Tis.s vs bch.Tis.s','mng.Tis.s vs bch.Tis.s', #habitat T sum
  'crk.Gw.w vs mng.Gw.w','crk.Gw.w vs bch.Gw.w', 'mng.Gw.w vs bch.Gw.w', # hab G win
  'crk.Gw.s vs mng.Gw.s', 'crk.Gw.s vs bch.Gw.s', 'mng.Gw.s vs bch.Gw.s', #hab G sum
  'crk.Tis.w vs crk.Gw.w', 'mng.Tis.w vs mng.Gw.w', 'bch.Tis.w vs bch.Gw.w', # area win
  'crk.Tis.s vs crk.Gw.s', 'mng.Tis.s vs mng.Gw.s','bch.Tis.s vs bch.Gw.s' #area sum
)
for (ind in 1:length(pairwise)){
  filtered_rows <- pairwise[[ind]][
    pairwise[[ind]][,1] %in% pair_ord,]
  
  tests$pairwise[[ind]] <- filtered_rows[match(pair_ord, filtered_rows[, 1]), ]
}
names(tests$pairwise)<- names(pairwise)


tuk_ord <- c(
  's:Tis:crk-w:Tis:crk', 's:Tis:mng-w:Tis:mng', 's:Tis:bch-w:Tis:bch', #seasonal Tis
  's:Gw:crk-w:Gw:crk', 's:Gw:mng-w:Gw:mng', 's:Gw:bch-w:Gw:bch',#seasonal G
  'w:Tis:mng-w:Tis:crk', 'w:Tis:bch-w:Tis:crk', 'w:Tis:bch-w:Tis:mng',#habitat T win
  's:Tis:mng-s:Tis:crk', 's:Tis:bch-s:Tis:crk', 's:Tis:bch-s:Tis:mng',#habitat T sum
  'w:Gw:mng-w:Gw:crk', 'w:Gw:bch-w:Gw:crk', 'w:Gw:bch-w:Gw:mng',# hab G win
  's:Gw:mng-s:Gw:crk', 's:Gw:bch-s:Gw:crk', 's:Gw:bch-s:Gw:mng', #hab G sum
  'w:Gw:crk-w:Tis:crk', 'w:Gw:mng-w:Tis:mng', 'w:Gw:bch-w:Tis:bch', # area win
  's:Gw:crk-s:Tis:crk', 's:Gw:mng-s:Tis:mng', 's:Gw:bch-s:Tis:bch'#area sum
)

for (ind in 1:length(tuk)){
  
  filtered_rows <- tuk[[ind]][[7]][
    row.names(tuk[[ind]][[7]]) %in% tuk_ord,]
  tests$tuk[[ind]]<- filtered_rows[match(tuk_ord, row.names(filtered_rows)), ]
}

names(tests$tuk)<- names(tuk)
  

############### excel output ############################################
library(openxlsx)

output<- tests

# Create a new workbook
wb <- createWorkbook()

# Add ANOVA results to the workbook
addWorksheet(wb, "ANOVA")
anova_data <- do.call(rbind, lapply(names(output$anova), function(name) {
  anova_summary <- summary(output$anova[[name]])
  df <- as.data.frame(anova_summary[[1]])
  df$Variable <- name
  df$RowName <- rownames(df)
  return(df)
}))
writeData(wb, sheet = "ANOVA", anova_data)

# Add Tukey test results to the workbook
addWorksheet(wb, "Tukey")
tukey_data <- do.call(rbind, lapply(names(output$tuk), function(name) {
  tukey_summary <- output$tuk[[name]]
  df <- as.data.frame(tukey_summary)
  df$Variable <- name
  df$RowName <- rownames(df)
  return(df)
}))
writeData(wb, sheet = "Tukey", tukey_data)

# Add PERMANOVA results to the workbook
addWorksheet(wb, "PERMANOVA")
peranova_data <- do.call(rbind, lapply(names(output$permanova), function(name) {
  df <- as.data.frame(output$permanova[[name]])
  df$Variable <- name
  df$RowName <- rownames(df)
  return(df)
}))
writeData(wb, sheet = "PERMANOVA", peranova_data)

# Add pairwise test results to the workbook
addWorksheet(wb, "Pairwise")
pairwise_data <- do.call(rbind, lapply(names(output$pairwise), function(name) {
  df <- as.data.frame(output$pairwise[[name]])
  df$Variable <- name
  return(df)
}))
writeData(wb, sheet = "Pairwise", pairwise_data)

# Save the workbook to a file (Save only one according to the initial setting)
saveWorkbook(wb, "output/env_anal.xlsx", overwrite = TRUE) # all
# saveWorkbook(wb, "output/env_anal_s.xlsx", overwrite = TRUE) # only seasonal


# saveWorkbook(wb, "output/env_anal_s_-bch.xlsx", overwrite = TRUE) # only seasonal, without beach
# saveWorkbook(wb, "output/env_anal_-bch.xlsx", overwrite = TRUE) # all, without beach

 