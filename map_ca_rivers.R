
## adapted from R Peek's code
## https://ryanpeek.org/2016-09-28-static_maps_in_r/

library(viridis) # nice color palette
library(ggplot2) # plotting
library(ggmap) # ggplot functionality for maps
library(ggsn) # for scale bars/north arrows in ggplots
library(maps)
library(mapdata)
library(cowplot)
library(sf)
library(data.table)
library(scatterpie)


## Make coordinate data spatial
# read in data
sample_coords_admix<-read.csv("refs/CCGP-trout_sample-coordinate-match.csv", stringsAsFactors = FALSE)
sample_coords_admix_nonNA <- subset(sample_coords_admix,!is.na(longitude))

# make this a spatial data frame 
# (i.e., spatial data needs to be in numeric form (use UTM, lat/long, or X/Y columns, doesn't matter you can convert later) 
# using lat/long
# need to figure out what the geodetic datum (or reference for representing the positions of locations on Earth)
# these are the CRS datum/projections "strings"
# in this case, will use the default for lat/long
# using the sf package instead of sp
# https://stackoverflow.com/questions/68157684/set-crs-for-latitude-longitude-point-data
sites.sf <- st_as_sf(sample_coords_admix_nonNA, coords=c("longitude", "latitude"), crs="EPSG:4326")
# geodetic CRS: WGS 84 (EPSG:4326)
st_crs(sites.sf)

# read in shapefile
ca_rivers <- st_read("NHD_Major_Rivers_and_Creeks/NHD_Major_Rivers_and_Creeks.shp")
# geodetic CRS: NAD83 (NAD83 + NAVD88 height) 
st_crs(ca_rivers)
# transform to add projection so that the map matches the coordinates
ca_rivers <- st_transform(ca_rivers, 4326)
st_crs(ca_rivers)
ca_rivers.sf_Rivers <- ca_rivers[ca_rivers$gnis_name %like% "River", ]

# convert sf objects to spatial object for downstream mapmaking
sites.sp <- as(sites.sf,"Spatial")
ca_rivers.sp <- st_zm(ca_rivers)
ca_rivers.sp <- as(ca_rivers.sp,"Spatial")
# subset rivers
library(data.table)
ca_rivers.sp_Rivers <- ca_rivers.sp[ca_rivers.sp$gnis_name %like% "River", ]


## Now map

# pull HEX color codes from iridis palette
pal <- viridis_pal()(6)
# pal  "#440154FF" "#414487FF" "#2A788EFF" "#22A884FF" "#7AD151FF" "#FDE725FF"

# this is using the `maps` and `mapsdata` package
# start by plotting the state/county

map("state",region=c('CA'), xlim = c(-122.5,-119.5), ylim=c(38,40.5))
map.axes()
map("county",region=c('CA'),boundary=FALSE,lty=3,
    xlim = c(-122.5,-119.5), ylim=c(38,40.5), add=TRUE)
lines(ca_rivers.sp_Rivers, col="#2A788EFF", lwd=1.1, xlim = c(-122.5,-119.5), ylim=c(38,40.5))

# add point data from earlier
points(sites.sp, cex=1.4, pch=21, bg="#FDE725FF")

# add labels (using trial and error for placement)
#text(sites.sp, labels=as.character(sites.sp@data$coord_name), col="gray20", cex=0.6, font=2, offset=0.5, adj=c(0,2))

# add a plot legend
legend("topright", border = "gray80",box.lwd = 1.5, box.col = "gray80",
       legend=c("Rivers", "Study Sites"),
       title="Study Map", bty="o", inset=0.05,
       lty=c( 1,-1,-1), pch=c(-1,21, 1), cex=0.8,
       col=c("#2A788EFF", "black"), pt.bg=c(NA, "#FDE725FF"))

# add a white box for the scale bar
rect(-122.47, ybottom = 38.02, xright = -121.52, ytop = 38.06, col = "white", border = "transparent")

# add north arrows
library(prettymapr)
addnortharrow(
  pos = "topleft",
  padin = c(0.15, 0.15),
  scale = 1,
  lwd = 1,
  border = "black",
  cols = c("white", "black"),
  text.col = "black"
)
addscalebar(
  plotunit = NULL,
  plotepsg = NULL,
  widthhint = 0.25,
  unitcategory = "metric",
  htin = 0.1,
  padin = c(0.15, 0.15),
  style = "bar",
  bar.cols = c("black", "white"),
  lwd = 1,
  linecol = "black",
  tick.cex = 0.7,
  labelpadin = 0.08,
  label.cex = 0.8,
  label.col = "black",
  pos = "bottomleft"
)
#map.scale(xc = -122.05, yc=38.3, len=0.7, units = "km",subdiv = 5, sfcol = "white", ndivs = 3)



## with stadia map base

# register stadia API
register_stadiamaps(key = "0ecd68c8-88ed-4c05-ae11-e18e5613cd58")

# set the bounding box of the map
# c(lowerleftlon, lowerleftlat, upperrightlon, upperrightlat)
bbox <- c(-122.5,38.0,-119.5,40.5) 

map_base <- get_stadiamap(bbox=bbox,crop = F,
                color="bw",
                maptype="stamen_terrain",
                zoom=8)

sitemap <- ggmap(map_base, extent = 'device') 
sitemap

# make river data spatial for ggplot
#rivers_df <- fortify(ca_rivers.sf_Rivers)
ca_rivers.sf_Rivers_coords.df <- as.data.frame(sf::st_coordinates(ca_rivers.sf_Rivers))


# we can just use the dataframe "sample_coords_admix_nonNA" with the lon/lat columns 
# (since we know they are in WGS84)
# If we wanted to use the sites.SP we would need to fortify first

nicemap<-
  sitemap + 
  labs(x="Longitude (WGS84)", y="Latitude") + 
  geom_path(data=ca_rivers.sf_Rivers_coords.df, aes(X, Y, group=L1, fill=NULL), 
            color="#2A788EFF", alpha=0.7) + 
  geom_label(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude, label=coord_name), 
             nudge_x=0.17, nudge_y=0.05, label.size=0.1,size=3, 
             fontface = "bold.italic", label.r=unit(0.20, "lines"))+
  geom_point(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude), pch=21, size=4, fill="#FDE725FF")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

nicemap


### overlay the ancestry piecharts

## california regional subset
nicemap_ancestries <-
  sitemap + 
  labs(x="Longitude (WGS84)", y="Latitude") + 
  geom_path(data=ca_rivers.sf_Rivers_coords.df, aes(X, Y, group=L1, fill=NULL), 
            color="#2A788EFF", alpha=0.7) + 
#  geom_label(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude, label=coord_name), 
#             nudge_x=0.17, nudge_y=0.05, label.size=0.1,size=3, 
#             fontface = "bold.italic", label.r=unit(0.20, "lines"))+
  geom_scatterpie(aes(x=longitude, y=latitude, group=region), data=sample_coords_admix_nonNA, cols=paste0("q",1:chosen_K), color=NA) +
#  geom_point(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude), pch=21, size=4, fill="#FDE725FF")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
nicemap_ancestries



## all of california

# set the bounding box of the map
# c(lowerleftlon, lowerleftlat, upperrightlon, upperrightlat)
ca_bbox <- c(-125,32.0,-114.5,42.5) 

ca_map_base <- get_stadiamap(bbox=ca_bbox,crop = F,
                          color="bw",
                          maptype="stamen_terrain",
                          zoom=8)

ca_sitemap <- ggmap(ca_map_base, extent = 'device') 
ca_sitemap

# grab river coordinates, put into dataframe
ca_rivers.sf_Rivers_coords.df <- as.data.frame(sf::st_coordinates(ca_rivers.sf_Rivers))

# prep the admixture piecharts
# read the sample_list.txt file; the qopt file is ordered by the sample order in this file
ngsadmix_dir <- "ngsadmix/filt_snps05-thin100_0-maf05" # set NGSadmix outputs directory
sample_list <- read.table(file.path(ngsadmix_dir,"sample_list.txt"))
colnames(sample_list) <- c("sample_ID")

# read pop and coordinate metadata
metadata <- read.csv("refs/CCGP-trout_sample-coordinate-match.csv",header=T)

# merge sample list with metadata (in order of sample list)
sample_info <- merge(sample_list,metadata,by="sample_ID")

# read in the qopt ancestry proportions 
chosen_K   <- 8
chosen_rep <- 1
admix_props <- read.table(file.path(ngsadmix_dir,paste0("K_",chosen_K,"_rep_",chosen_rep),"output.qopt"))
names(admix_props) = paste0("q",1:chosen_K)

# merge sample-coordinate metadata with ancestry proportions data
sample_coords_admix <- as.data.frame(cbind(sample_info, admix_props))

# temporary prep steps
sample_coords_admix_nonNA <- subset(sample_coords_admix,!is.na(longitude))
sample_coords_admix_nonNA$region <- factor(1:nrow(sample_coords_admix_nonNA))

# remove the hatchery fish
sample_coords_admix_nonNA_nohatchery <- subset(sample_coords_admix_nonNA, subspecies_name != "hatchery")

ca_nicemap_ancestries <-
  ca_sitemap + 
  labs(x="Longitude (WGS84)", y="Latitude") + 
  geom_path(data=ca_rivers.sf_Rivers_coords.df, aes(X, Y, group=L1, fill=NULL), 
            color="#2A788EFF", alpha=0.7) + 
  #  geom_label(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude, label=coord_name), 
  #             nudge_x=0.17, nudge_y=0.05, label.size=0.1,size=3, 
  #             fontface = "bold.italic", label.r=unit(0.20, "lines"))+
  geom_scatterpie(aes(x=longitude, y=latitude, group=region), data=sample_coords_admix_nonNA_nohatchery, cols=paste0("q",1:chosen_K), pie_scale=0.7, color=NA) +
  #  geom_point(data=sample_coords_admix_nonNA, aes(x=longitude, y=latitude), pch=21, size=4, fill="#FDE725FF")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
ca_nicemap_ancestries
ggsave(ca_nicemap_ancestries, filename='Figures/ccgp-trout_admixture_map_stadia-bw.pdf', device = cairo_pdf, width=6, height=8, bg="transparent")




## extra
## CA rivers and creeks
library(sf)
ca_rivers <- st_read("NHD_Major_Rivers_and_Creeks/NHD_Major_Rivers_and_Creeks.shp")
ca_rivers_nozm <- st_zm(ca_rivers) # 324886 features and 16 fields
# can show specific features by subsetting them like this
print(ca_rivers_nozm[9:15], n = 3) 
ca_rivers_base <- ggplot(ca_rivers_nozm) + geom_sf() # takes a while

## Historic trout ranges
historic_trout_sfdf <- st_read("Historic_Trout_Watersheds_[ds440]/Historic_Trout_Watersheds_[ds440].shp")
plot(st_geometry(historic_trout_sfdf))
plot(historic_trout_sfdf["Species"])
ggplot() + geom_sf(historic_trout_sfdf, aes(fill = "Species"))  


# https://wifire-data.sdsc.edu/dataset/california-river-and-stream-hydrography/resource/a16a75a2-c5f3-4cb5-a15b-8fa414566a59
ca_hydrography <- st_read("California_River_and_Stream_Hydrography/California_River_and_Stream_Hydrography.shp")
ggplot(ca_hydrography) + geom_sf() # takes a while


# if you want at some point, you can use ArcMap to make your own shapefiles from
# google earth images
# https://gis.stackexchange.com/questions/313662/how-can-i-generate-shapefile-from-google-earth

