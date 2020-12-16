library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("gghighlight")

world <- ne_countries(scale = "small", returnclass = "sf")

lakes = read.csv("coords2.csv", h=T, col.names=c("Lake", "longitude", "latitude"))

lakes_sf = st_as_sf(x = lakes, coords = c("longitude", "latitude"), crs=  "+proj=longlat")
lakes_sf = data.table(lakes_sf)
lakes_derep = lakes_sf[,.(nb_bases = as.numeric(sum(nb_bases)), taxon = taxon[1], types = length(levels(factor(taxon))), nb_samples = length(taxon)),  by=coord]


ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=lakes_sf, col="violetred1", size=5, shape=18)+
  coord_sf(xlim=c(5, 30), ylim=c(57,70))+
  geom_sf_text_repel(data=lakes_sf, aes(label=Lake), size=3, force = 20, nudge_x = -4, seed = 17)+xlab("")+ylab("")
ggsave("scandi_map.pdf")
