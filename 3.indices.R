# calculating the species diversity, functional diversity, and ecological quality indices

#Pedram Ghahramani pedram.ghrm@gmail.com
#fall 2023

#libraries
library(dplyr) #for filtering at the end
library(vegan)  # species diversity
# library(gawdis)   # gower dissimilarity
library(picante)  # functional diversity 
#library(ade4)    # functional diverity "divc" function for Rao
#library(data.table) # functional diverity
library(FD)         # functional diverity
library(adiv)   # for Rao and Redundancy
  ##importing raw data ----
ab_raw<- read.csv("raw data/ab.csv", header = T, row.names = 1)# abundance of all species (for taxonomic diversity)
# traits<- read.csv("raw data/spXtraits.csv", header = T, row.names = 1)#sp-trait matrix
ab_t<- read.csv("raw data/ab(for traits).csv", header = T, row.names = 1)# abundance with exclusion of rare species
ecoQ<- read.csv("raw data/EcoQab.csv", header = T, row.names = 1) # abundance of each ecoQ group
M_AMBI<- read.csv("raw data/M-AMBI.csv", header = T, row.names = 1) # M-AMBI index and grouping (computed in "Azti.exe" software)
dist<-as.dist(read.csv("output/dist.spxtrait(gowdis).csv", row.names = 1))
#dist<- read.csv("output/dist.spxtrait(gawer).csv", header = T, row.names = 1)
#dist <- read.csv("output/dist.spXtrait(bray).csv", header = T, row.names = 1) #bray cusrtis species trait distance (computed in another script)

#deleting factor columns and transformation
ab<- ab_raw[-((ncol(ab_raw)-3):ncol(ab_raw))]
ab.L<-log(1 + ab)
ab_t.L <- log(1+ ab_t)

  ##Species diversity ----

#without transformation of the abundance
S<- specnumber(ab) # Species richness (S) 
N<- rowSums(ab) #total abundance
Smg=(S-1)/log(N) #Margalef’s diversity index
H=diversity(ab) #Shannon-Wiener’s diversity index
D=diversity(ab,"simpson") #Simpson’s diversity index
J=H/log(S) #Pielou’s evenness index
H<- replace(H, H == 0,  NA)
#with log transformation of the abundance
S.L<- specnumber(ab.L) # Species richness (S) 
N.L<- rowSums(ab.L) #total abundance
Smg.L<- (S.L-1)/log(N.L) #Margalef’s diversity index
H.L<- diversity(ab.L) #Shannon-Wiener’s diversity index
D.L<- diversity(ab.L,"simpson") #Simpson’s diversity index
J.L<- H.L/log(S.L) #Pielou’s evenness index
H.L<- replace(H.L, H.L == 0,  NA)
  ##functional diversity ----

    #IF WE DON'T USE BRAY CURTIS WE CAN USE GOWER, COMPUTED AS BELLOW
#trait distance between species
#we use gower distance for traits with different number of modalities
# dist<- (gowdis(traits[1]) + gowdis(traits[2:6])/max(gowdis(traits[2:6]))+ 
#           gowdis(traits[7:10])/max(gowdis(traits[7:10])) + 
#           gowdis(traits[11:13])/max(gowdis(traits[11:13])) + 
#           gowdis(traits[14:17])/max(gowdis(traits[14:17])) + 
#           gowdis(traits[18]) + gowdis(traits[19]) + 
#           gowdis(traits[20:22])/max(gowdis(traits[20:22])) + 
#           gowdis(traits[23:25])/max(gowdis(traits[23:25])) + 
#           gowdis(traits[26:29])/max(gowdis(traits[26:29])) )/9


#without transformation of the abundance
# mpd<- mpd (ab_t, as.matrix(dist), abundance.weighted = F)#Mean pairwise distance
# mntd<- mntd (ab_t, as.matrix(dist), abundance.weighted = F)#Mean nearest taxon distance
FD<- dbFD(dist,ab_t, message = F, 
          calc.CWM = FALSE, stand.FRic = TRUE) #sp.num, FRich, FEve, FDiv, FDis, Rao

#with transformation 
# mpd.L<- mpd (ab_t.L, as.matrix(dist), abundance.weighted = F)#Mean pairwise distance
# mntd.L<- mntd (ab_t.L, as.matrix(dist), abundance.weighted = F)#Mean nearest taxon distance
FD.L<- dbFD(dist, ab_t.L, message = F, 
            calc.CWM = FALSE, stand.FRic = TRUE)

  ## ecological quality  ----

EG=ecoQ[1:4] #AMBI groups
G=ecoQ[c(5,6)] #BENTIX groups
AMBI<- (0*EG[,1])+(1.5*EG[,2])+(3*EG[,3])+(4.5*EG[,4])
BENTIX<- 6*G[,1]+2*G[,2]

#M-AMBI
#it is calculated by combining the AMBI, Shannon, and species richness
#first we must standardize the each variable
##std.data1 <- data.frame(AMBI=scale(AMBI),shannon= scale(H),richness= scale(S))
#pc=prcomp(std.data)
#library(psych)
#pairs.panels(indices,gap=0,bg = c("red", "yellow", "blue","purple","orange","violet","brown"),pch=21)

####### functional redundancy #############
##indices with adive package
adiv<- uniqueness(ab_t, dist)[[3]]
#log transformed
adiv.L<- uniqueness(ab_t.L, dist)[[3]]


##### station with less that 2 species should be eliminated 
adiv$Q[specnumber(ab_t) == 1]<- NA
adiv.L$Q[specnumber(ab_t) == 1]<- NA
AMBI[specnumber(ab) == 1]<- NA
#export  ----

indices<- data.frame( S, S2 = adiv$N,  H , D, J ,Dm= Smg,
                      # mpd, mntd, 
                      FRic=FD$FRic, FEve= FD$FEve,
                      Rao = adiv$Q, redun = adiv$R,
                      AMBI, BENTIX, M_AMBI)

indices.L <- data.frame(S= S.L, S2 = adiv.L$N, H = H.L ,D= D.L,
                        J= J.L ,Dm= Smg.L,
                        # mpd.L,mntd.L,
                        FRic = FD.L$FRic, FEve = FD.L$FEve,
                        Rao = adiv.L$Q, redun = adiv.L$R,  
                        AMBI, BENTIX, M_AMBI)

# indices[indices$S <3,]$AMBI <- 7
# indices.L[indices.L$S <3,]$AMBI <- 7
#adding factor columns
indices<- cbind(indices,ab_raw[(ncol(ab_raw)-3):ncol(ab_raw)])
indices.L<- cbind(indices.L, ab_raw[(ncol(ab_raw)-3):ncol(ab_raw)])

#excliding plots with richness of 1
#indices<- indices[indices$S!=1,]
#indices.L<- indices.L[indices.L$S!=1,]

#replacing NA with 0
# indices[is.na(indices)]<- 0
# indices.L[is.na(indices.L)]<- 0


write.csv(indices, file = "output/indices.CSV")
write.csv(indices.L, file = "output/indices(Log).CSV")


