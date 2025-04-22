# stations functional dissimilarity (Bray)
#PED
# 1/21/2024
library(tidyverse)
library(vegan)
library(gawdis)

#input the trait-species and species-community matrices
sp.trait<- read.csv("raw data/spXtraits.csv", row.names = 1) 
  #first trait is nominal, 6th and 7th are ordinal, the rest are fuzzy
sp.trait1<- read.csv("raw data/spXtraits(new).csv", row.names = 1)
  #fist ,numinal, 2nd to 5th area fuzzy, 6th and 7th are ordinal, 8th, 9th and 10th are categorical

rownames(sp.trait1) <- rownames(sp.trait)
ab <- read.csv("raw data/ab(for traits).csv",row.names = 1)


# calculating community-trait matrix
com.trait <- as.data.frame( log(1+(as.matrix(ab) %*% as.matrix(sp.trait))))

################################### gower distance ############################

    #using gowdis() function from FD package
feed<- sp.trait%>%select(starts_with("FD"))
move<- sp.trait%>%select(starts_with("MV"))
hab<- sp.trait%>%select(starts_with("HB"))
  size<- sp.trait%>%select(starts_with("S."))
  bioturb<- sp.trait%>%select(starts_with("BT"))
  repT<- sp.trait%>%select(starts_with("RS"))
  Tough<- sp.trait%>%select(starts_with("T"))
  Flex<- sp.trait%>%select(starts_with("FL"))
  Form<- sp.trait%>%select(starts_with("FR"))
  LS<- sp.trait%>%select(starts_with("LS"))
  w<- sp.trait%>%select(starts_with("w"))
  w1<- sp.trait1%>%select(starts_with("w"))
  
dis.sp.trait<- ( 
                  gowdis(w1) + # weight
                  gowdis(feed)/max(gowdis(feed))+ # Feeding 
                  gowdis(move)/max(gowdis(move)) + # Movement
                 # gowdis(sp.trait[11:13])/max(gowdis(sp.trait[11:13])) + # Habitat
                  gowdis(bioturb)/max(gowdis(bioturb)) + # Bioturbation
                  gowdis(Tough) + #Toughness
                  # gowdis(Flex) + # flexibility
                  gowdis(repT)/max(gowdis(repT)) +# Reporductive technique
                  # gowdis(sp.trait[26:29])/max(gowdis(sp.trait[26:29]))+ # Form
                   # gowdis(size)/max(gowdis(size))+ #size
                   gowdis(LS)/max(gowdis(LS)) # lifespan
                )/7



feed<- com.trait%>%select(starts_with("FD"))
move<- com.trait%>%select(starts_with("MV"))
hab<- com.trait%>%select(starts_with("HB"))
size<- com.trait%>%select(starts_with("S."))
bioturb<- com.trait%>%select(starts_with("BT"))
repT<- com.trait%>%select(starts_with("RS"))
Tough<- com.trait%>%select(starts_with("T"))
Flex<- com.trait%>%select(starts_with("FL"))
Form<- com.trait%>%select(starts_with("FR"))
LS<- com.trait%>%select(starts_with("LS"))
w<- com.trait%>%select(starts_with("w"))

dis.com.trait<- (
  gowdis(log(1+w)) + # weight
  gowdis(feed)/max(gowdis(feed))+ # Feeding 
    gowdis(move)/max(gowdis(move)) + # Movement
    # gowdis(sp.trait[11:13])/max(gowdis(sp.trait[11:13])) + # Habitat
    gowdis(bioturb)/max(gowdis(bioturb)) + # Bioturbation
    gowdis(Tough) + #Toughness
    # gowdis(Flex) + # flexibility
    gowdis(repT)/max(gowdis(repT)) + # Reporductive technique
    # gowdis(sp.trait[26:29])/max(gowdis(sp.trait[26:29]))+ # Form
    # gowdis(size)/max(gowdis(size)) #size
    gowdis(LS)/max(gowdis(LS)) # lifespan
)/7

dis.sp.trait <- as.matrix(dis.sp.trait)
dis.com.trait <- as.matrix(dis.com.trait)

    #using gawdis() function from gawdis package

#species-trait distance
#dis.sp.trait0 <-as.matrix(round( gawdis(sp.trait, w.type = "equal",
  #                  groups = rep(1:10, c(1,5,4,3,4,1,1,3,3,4) ) , fuzzy = c(2,3,4,5,8,9,10)), 3))

#community-trait distance
#dis.com.trait0 <-as.matrix(round( gawdis(com.trait, w.type = "equal",
  #                 groups = rep(1:10, c(1,5,4,3,4,1,1,3,3,4)), fuzzy = c(2,3,4,5,8,9,10)), 3))

###################################### Bray curtis ############################

# #species-trait distance
# #compute each trait separately
# S <- vegdist(sp.trait[1]) #wight
# FD <- vegdist(sp.trait[2:6]) #feeding
# M <- vegdist(sp.trait[7:10]) #movement
# # H <- vegdist(sp.trait[11:13]) #habitat
# B <- vegdist(sp.trait[14:17]) #bioturbation
# t <- vegdist(sp.trait[18]) #shell toughness
# FL <- vegdist(sp.trait[19]) #flexibility
# # LD <- vegdist(sp.trait[20:22]) #Larval developement
# RT <- vegdist(sp.trait[23:25]) #Reproductive technique
# #SD <- vegdist(sp.trait[35:36]) #sexual differensiation
# FR <- vegdist(sp.trait[26:29]) #Form
# #RS <- vegdist(sp.trait[42:43]) #semelparity
# 
# dis.spXt.list<- list(as.matrix(S), as.matrix(FD), as.matrix(M), as.matrix(H), as.matrix(B),
#             as.matrix(t), as.matrix(FL), as.matrix(LD), as.matrix(RT), as.matrix(FR))
# mean.dist.spXt <- as.data.frame(apply(simplify2array(dis.spXt.list), c(1, 2),
#                                            mean, na.rm = T), 2)
# rownames(mean.dist.spXt)<- colnames(mean.dist.spXt)
# 
# #mean.dall <- (S+FD+M+H+B+t+FL+LD+RT+SD+FR+RS)/12 #another way (if no NA)
# 
# #community-trait distance
# #the same process as species-trait
# S1 <- vegdist(com.trait[1])
# FD1 <- vegdist(com.trait[2:6])
# M1 <- vegdist(com.trait[7:10])
# H1 <- vegdist(com.trait[11:13])
# B1 <- vegdist(com.trait[14:17])
# t1 <- vegdist(com.trait[18])
# FL1 <- vegdist(com.trait[19])
# # LD1 <- vegdist(com.trait[20:22])
# RT1 <- vegdist(com.trait[23:25])
# #SD1 <- vegdist(com.trait[35:36])
# FR1 <- vegdist(com.trait[26:29])
# #RS1 <- vegdist(com.trait[42:43])
# 
# dis.comXt.list<- list(as.matrix(S1), as.matrix(FD1), as.matrix(M1), as.matrix(H1), as.matrix(B1),
#             as.matrix(t1), as.matrix(FL1), #as.matrix(LD1),
#             as.matrix(RT1), as.matrix(FR1))
# mean.dist.comXt <- as.data.frame(apply(simplify2array(dis.comXt.list), c(1, 2),
#                            mean, na.rm = T), 2)
# rownames(mean.dist.comXt)<- colnames(mean.dist.comXt)

  #export ################################################

write.csv(com.trait ,"output/comXtrait(LOGab).csv") # community-trait matrix
#write.csv(dis.sp.trait0 , "output/dist.spXtrait(gawdis).csv") #species-trait distance  with gower
write.csv(dis.sp.trait, "output/dist.spXtrait(gowdis).csv") #species-trait distance  with gower
#write.csv(dis.com.trait0 ,"output/dist.comXtrait(gawdis).csv") #community-trait distance with gower
write.csv(dis.com.trait ,"output/dist.comXtrait(gowdis).csv") #community-trait distance with gower
# write.csv(mean.dist.spXt , "output/dist.spXtrait(bray).csv") #species-trait distance with bray curtis
# write.csv(mean.dist.comXt ,"output/dist.comXtrait(bray).csv") #community-trait distance with bray curtis
