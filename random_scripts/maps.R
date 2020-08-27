library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("gghighlight")

world <- ne_countries(scale = "small", returnclass = "sf")
proj="+proj=laea +lat_0=45 +lon_0=-30 +ellps=WGS84 +units=m"

lakes = read.csv("coords2.csv", h=T, col.names=c("Lake", "longitude", "latitude"))
world_pretty = st_transform(world, crs = proj)

lakes_sf = st_as_sf(x = lakes, coords = c("longitude", "latitude"), crs=  "+proj=longlat")
lakes_pretty = st_transform(lakes_sf, crs = proj)

world_label_data =  lakes_sf[lakes_sf$Lake %in% c("La Plata","Loclat","Toolik","B4"),]
world_label_data$Lake = c("Kuujjuarapik-Whapmagoostui",as.vector(world_label_data$Lake)[2:4])

make_box = function(blx,bly,tlx,tly,brx,bry,trx,try){
  box = st_sfc(st_linestring(matrix(c(tlx, blx,brx, trx, tlx, tly,bly, bry,try, tly), ncol=2)), crs = "+proj=longlat")
  box = st_transform(box,proj)
  box
}

box = st_sfc(st_linestring(matrix(c(5, 5,30,30,5,70,57, 57,70,70), ncol=2)), crs = "+proj=longlat")
box = st_transform(box,proj)


zoom = make_box(-100,-10,-170, 20, 50,-15, 120,20)
xwind = c(min(st_coordinates(lakes_pretty)[,"X"]), max(st_coordinates(lakes_pretty)[,"X"]))*1.2
ywind = c(min(st_coordinates(lakes_pretty)[,"Y"]), max(st_coordinates(lakes_pretty)[,"Y"]))*1.2

ggplot(data = world_pretty) + geom_sf(fill="floralwhite")+geom_sf(data=lakes_pretty, col="violetred1", size=8, shape=18)+
  geom_sf_text_repel(data=world_label_data, aes(label=Lake),force =40, size=7, seed = 23,segment.size=1, nudge_x=-1000000, nudge_y=-500000)+
  xlab("")+ylab("")+geom_sf(data=box,size=1)+coord_sf(xlim=xwind, ylim=ywind)
ggsave("world_map.pdf")

canada_label_data =  lakes_sf[lakes_sf$Lake %in% c("B1_2", "B4", "C2_4", "C5", "F5", "G1", "H1", "I1_2", "I3", "SAS2A", "SAS2C", "SAS2D",'SAS2B'),]

ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=lakes_sf, col="violetred1", size=5, shape=18)+coord_sf(xlim=c(-77.698, -77.694), ylim=c(55.225,55.227))+
geom_sf_text_repel(data=canada_label_data, aes(label=Lake),force =30, size=4, seed = 42)+xlab("")+ylab("")
ggsave("canada_map.pdf")

scandi_label_data =  lakes_sf[lakes_sf$Lake %in% c("Malstasjön", "Bengtgölen", "Långsjön", "Odrolstjärnen", "Alinen Mustajärvi", "Haukijärvi",
"Keskinen Rajajärvi", "Mekkojärvi", "Ki1", "Ki2", "Fyrsån", "Valkea Kotinen", "Lillsjön", "Erken", "754378-169136", "Nästjärnen", "Stortöveln", "Parsen", "Björntjärnen West",
"Plåten", "Glimmingen", "Lumpen", "Lovojärvi","Björntjärnen East", "Lotsjön", "Lomtjärnan", "Ylinen Rajajärvi"),]

ggplot(data = world) + geom_sf(fill="floralwhite", alpha=0.5, col="lightgray")+geom_sf(data=lakes_sf, col="violetred1", size=5, shape=18)+
  coord_sf(xlim=c(5, 30), ylim=c(57,70))+
  geom_sf_text_repel(data=scandi_label_data, aes(label=Lake), size=3, force = 20, nudge_x = -4, seed = 17)+xlab("")+ylab("")
ggsave("scandi_map.pdf")
