#ped
#1-12-2024
#cwm of my traits (each modality)

library(FD)

#data input
traits<- read.csv("raw data/spxtraits(new).csv", row.names = 1)
ab <- read.csv("raw data/ab(for traits).csv",row.names = 1)

#transform for highly skewed numerical traits
traits$Weight <- log(1+traits$Weight)

#cumputation of cwm
cwm <- functcomp(traits, log(1+as.matrix(ab)), CWM.type = "all")

#scaling all traits to be between 0 and 1
cwm[18] <- cwm[18]/max(cwm[18])
cwm[19] <- cwm[19]/max(cwm[19])
cwm[2:17] <- cwm[2:17]/5

cwm<- round(cwm,4)

#output
write.csv(cwm, "output/CWM.csv")


