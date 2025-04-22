# hypothesis test for CWM on trait modalities
# pedram.ghrm@gmail.com
# last mod 9/5/2024

# libraries
library(vegan)
library(tidyverse)
library(pairwiseAdonis)

# input data
source('R-scripts/factors.R') #factors
cwm <- read.csv("output/CWM.csv", row.names = 1)
cwm1 <- cbind(cwm,factors)

############ data filtering ######
# ######### seasonal filter
cwm1<- cwm1%>%
  filter(transect != "1", transect != "2", transect!= "7", transect!= "8")
# 
# # ######## filter bch out
cwm1<- cwm1%>%
  filter(habitat != "bch") #removing beach habitat


cwm<- cwm[row.names(cwm) %in% row.names(cwm1),]
factors<- factors[row.names(factors) %in% row.names(cwm1),]

############### outliers ##################################
## Function to detect outliars using Z-scores
# 
# out <- function(df, factor = "season.area.habitat", threshold = 3) {
#   groups<- unique(df[[factor]])
#   
#   outliars<- list()
#   for (i in groups){
#     outliar <- df[df[[factor]] == i,]%>%
#       select(1:13)%>%
#       mutate(across(everything(), ~ ifelse(abs(scale(.)) > threshold, ., NA)))
#     
#     outliars[[as.character(i)]] <- outliar
#     
#     # Get non-NA values and their indices
#     non_na_indices <- which(!is.na(outliar), arr.ind = TRUE)
#     non_na_values <- data.frame(
#       row = rownames(outliar)[non_na_indices[, 1]],
#       column = colnames(outliar)[non_na_indices[, 2]],
#       value = outliar[non_na_indices])
#     
#     outliars[[paste(as.character(i),"1")]] <- non_na_values
#     
#   }
#   
#   # all of the variables across all stations
#   all<- df%>%
#     select(1:13)%>%
#     mutate(across(everything(), ~ ifelse(abs(scale(.)) > threshold, ., NA)))
#   
#   non_na_indices <- which(!is.na(all), arr.ind = TRUE)
#   non_na_values <- data.frame(
#     row = rownames(all)[non_na_indices[, 1]],
#     column = colnames(all)[non_na_indices[, 2]],
#     value = all[non_na_indices])
#   
#   outliars[["all1"]] <- non_na_values
#   
#   
#   return(outliars)
# }

# outs<- out(cwm1, threshold = 2)[["all1"]]

# cwm1$Weight[row.names(cwm1) %in% outs[outs$column == "Weight",]$row]<- NA
# 
# 
# ggplot(cwm1, aes(season.area.habitat, Weight))+
#   geom_boxplot()

############################ PERMANOVA & Pairwise #########################

tests <- function(df){
  # ######################### ANOVA and Tukey's 
  # anov<- list()
  # tuk<- list()
  # ttest<- list()
  # for (i in names(df[1:13])){
  #   var<- df[c(i, 'season', 'area', 'habitat', 'season.area.habitat')]
  #   names(var)[1]<- "a"
  #   
  #   an<- aov(a ~ season * area * habitat, data = var)
  #   tu<- TukeyHSD(an)
  #   t<- pairwise.t.test(var$a, var$season.area.habitat, p.adjust.method = "none")
  #   
  #   anov[[length(anov)+1]]<- an
  #   tuk[[length(tuk)+1]]<- tu
  #   ttest[[length(ttest)+1]]<- t
  #   
  #   names(anov)[length(anov)]<- names(df[i])
  #   names(tuk)[length(tuk)]<- names(df[i])
  #   names(ttest)[length(ttest)]<- names(df[i])
  # }
  # 
  # 
  
  ############ post hoc for normal distribution
  # 1. tukey's Honest Significant Difference:
  # Ideal for comparing all possible pairs of groups.
  
  # 2.Bonferroni Correlation
  # A more conservative approach that adjusts p-values to control the family-wise error rate.
  # pairwise.t.test(df$H, df$season.area.habitat, p.adjust.method = "bonferroni")
  # print(pairwise_result)
  
  # Scheffé's Method:
  # Useful when you have unequal sample sizes or are making many comparisons.
  # library(DescTools)
  # scheffe_result <- ScheffeTest(aov_result)
  # print(scheffe_result)
  
  ##################################### permanova
  perm<- list()
  for (i in names(select_if(df, is.numeric))){
    
    var<- df[c(i, 'season', 'area', 'habitat')]
    names(var)[1]<- "a"
    var<- na.omit(var)
    a<- as.matrix(var[1])
    # a <- as.matrix(var[1])
    per <- adonis2(a ~ season * area * habitat, data = var,
                   method = "euclidean", permutations = 9999)
    
    perm[[length(perm)+1]]<- per
    names(perm)[length(perm)]<- names(df[i])
    
  }
  
  ################################### pair wise test 
  ## solution two (more data in the function)
  pairwise<- list()
  for (i in names(select_if(df, is.numeric))){
    var<- na.omit(df[c(i,"season.area.habitat")])
    names(var)[1]<- "a"
    dis<- dist(var$a)
    pair <- pairwise.adonis(dis, var$season.area.habitat, p.adjust.m='fdr',
                            sim.method = "euclidean", perm = 9999)
    
    pairwise[[length(pairwise)+1]]<- pair
    names(pairwise)[length(pairwise)]<- names(df[i])
  }
  
  return(list(#anova = anov, Tukey = tuk, t_test = ttest,
    permanova = perm, pairwise = pairwise))
}

test0<- tests(cwm1)

#################### for modalities combined (each tarit)#################
names<- c("w","FD","MV","HB","BT","TN","FL","RT","LD","FR","S.","LS" )
  
whole<- function(){
  perm1<- list()
  pairwise1<- list()
  for (i in names) {
    t<- cwm1%>%
      select(a = starts_with(i))
    per<- adonis2(t ~ season *area * habitat, data = cwm1,
        method = "euclidean", permutations = 9999)
    
    dis<- dist(t)
    pair <- pairwise.adonis(dis, cwm1$season.area.habitat, p.adjust.m='fdr',
                            sim.method = "euclidean", perm = 9999)
    
    pairwise1[[length(pairwise1)+1]]<- pair
    names(pairwise1)[length(pairwise1)]<- i
    
    perm1[[length(perm1)+1]]<- per
    names(perm1)[length(perm1)]<- i
  }
  
  return(list(#anova = anov, Tukey = tuk, t_test = ttest,
    permanova = perm1, pairwise = pairwise1))
}

test1<- whole()



############ selecting desired pairs and in order#################

tests0<- list()
tests1<- list()
tests0$permanova<- test0$permanova
tests1$permanova<- test1$permanova

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
for (ind in 1:length(test0$pairwise)){
  filtered_rows <- test0$pairwise[[ind]][
    test0$pairwise[[ind]][,1] %in% pair_ord,]
  
  tests0$pairwise[[ind]] <- filtered_rows[match(pair_ord, filtered_rows[, 1]), ]
}
for (ind in 1:length(test1$pairwise)){
  filtered_rows <- test1$pairwise[[ind]][
    test1$pairwise[[ind]][,1] %in% pair_ord,]
  
  tests1$pairwise[[ind]] <- filtered_rows[match(pair_ord, filtered_rows[, 1]), ]
}

names(tests0$pairwise)<- names(test0$pairwise)
names(tests1$pairwise)<- names(test1$pairwise)


############### excel output ############################################
library(openxlsx)

# Assuming test0 is the output list from the function tests
output <- tests0
output1 <- tests1
# Create a new workbook
wb <- createWorkbook()

# Add PERMANOVA results to the workbook
addWorksheet(wb, "PERMANOVA")
permanova_data <- do.call(rbind, lapply(names(output$permanova), function(name) {
  df <- as.data.frame(output$permanova[[name]])
  df$Variable <- name
  return(df)
}))
writeData(wb, sheet = "PERMANOVA", permanova_data)

addWorksheet(wb, "PERMANOVA1")
permanova_data1 <- do.call(rbind, lapply(names(output1$permanova), function(name) {
  df <- as.data.frame(output1$permanova[[name]])
  df$Variable <- name
  return(df)
}))
writeData(wb, sheet = "PERMANOVA1", permanova_data1)

# Add pairwise test results to the workbook
addWorksheet(wb, "Pairwise")

pairwise_data <- do.call(rbind, lapply(names(output$pairwise), function(name) {
  df <- as.data.frame(output$pairwise[[name]])
  df$Variable <- name
  return(df)
}))
writeData(wb, sheet = "Pairwise", pairwise_data)

addWorksheet(wb, "Pairwise1")
pairwise_data1 <- do.call(rbind, lapply(names(output1$pairwise), function(name) {
  df <- as.data.frame(output1$pairwise[[name]])
  df$Variable <- name
  return(df)
}))
writeData(wb, sheet = "Pairwise1", pairwise_data1)

# Save the workbook to a file
  ## all
# saveWorkbook(wb, "output/CWM_PERMANOVA.xlsx", overwrite = TRUE)
#   ## seasonal 
# saveWorkbook(wb, "output/CWM_PERMANOVA_s.xlsx", overwrite = TRUE)

#  ## all -bch
# saveWorkbook(wb, "C:/Users/pedi/docs/[01]Academia/[01]research papers/[1]Mangrove papers/past attempts/Without bch/output/CWM_PERMANOVA_-bch.xlsx", overwrite = TRUE)

  ## seasonal -bch
saveWorkbook(wb, "C:/Users/pedi/docs/[01]Academia/[01]research papers/[1]Mangrove papers/past attempts/Without bch/output/CWM_PERMANOVA_-bch_s.xlsx", overwrite = TRUE)

