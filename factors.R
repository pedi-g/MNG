# Community abundance matrix
ab_raw <- read.csv("raw data/ab.csv", row.names = 1)

factors <- data.frame(station = factor(ab_raw$station),
                    season = factor(ab_raw$season, levels = c("w", "s")),
                    area = factor(ab_raw$area, levels = c("Tis", "Gw")),
                    habitat = factor(ab_raw$habitat, levels = c("crk", "mng", "bch")))
rownames(factors)<- rownames(ab_raw)
factors$season.area <- interaction(factors$season, factors$area)
factors$season.habitat <- interaction(factors$season, factors$habitat)
factors$area.habitat <- interaction(factors$area, factors$habitat)
factors$season.area.habitat <- interaction(factors$season, factors$area, factors$habitat)


s <- as.character(factors$station)
m <- c()
for (i in seq(length(s))) {
  substring(s[i], 1, 1)
  m <- append(m, substring(s[i], 1, 1))
}
factors$transect <- as.factor(m)

rm(ab_raw,i,m,s)

