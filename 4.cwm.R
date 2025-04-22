#ped
#1-12-2024
#cwm of my traits (each modality)

library(FD)
source('R-scripts/factors.R') #factors
#data input
traits<- read.csv("raw data/spxtraits.csv", row.names = 1)
ab <- read.csv("raw data/ab(for traits).csv",row.names = 1)
# com.trait<- read.csv('output/comxtrait(LOGab).csv', row.names = 1)

#transform for highly skewed numerical traits
# hist(traits$Weight)
traits$Weight <- log(1+traits$Weight)
# hist(traits$Weight)
#cumputation of cwm
cwm <- functcomp(traits, log(1+as.matrix(ab)), CWM.type = "all")

#scaling all traits to be between 0 and 1
# cwm[18] <- cwm[18]/max(cwm[18])
# cwm[19] <- cwm[19]/max(cwm[19])
# cwm[2:17] <- cwm[2:17]/5
# cwm[27:30] <- cwm[27:30]/4
# 
# cwm<- round(cwm,4)
# 
# cwm$Weight<- com.trait$Weight
# 
# cwm[1] <- cwm[1]/max(cwm[1])
cwm1<- cbind(cwm, factors)
boxplot(Weight ~ Season.Site.Habitat, data = cwm1)

higher_outliers <- rownames(cwm[cwm[1]>0.24,])
lower_outliers <- rownames(cwm[cwm[1]<0.004,])

ncwm <- cwm
ncwm[rownames(cwm) %in% higher_outliers | rownames(cwm) %in% lower_outliers,1]<- NA
w<- cbind( wAll = cwm$Weight, wNew = ncwm$Weight, factors)

w<- w %>%
  filter(Habitat != "bch")

cwm1<- cbind(cwm, factors)
par(mfrow = c(1,2))
boxplot(wAll ~ Season.Site.Habitat, data = w)
boxplot(wNew ~ Season.Site.Habitat, data = w)


cwm$Weight<- ncwm$Weight
#output
write.csv(cwm, "output/CWM.csv")


