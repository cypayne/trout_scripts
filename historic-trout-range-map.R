# Read this shape file with the rgdal library. 
library(sf)
library(ggplot2)
library(ggmap) # ggplot functionality for maps
library(data.table)



## prepare historic trout ranges and site coordinates

# read in historic trout watershed shapefile
hist_trout_ranges <- st_read("Historic_Trout_Watersheds_[ds440]/Historic_Trout_Watersheds_[ds440].shp")
# geodetic CRS: WGS 84 
st_crs(hist_trout_ranges)
# transform to add projection so that the map matches the sample coordinates
hist_trout_ranges <- st_transform(hist_trout_ranges, 4326)
st_crs(hist_trout_ranges)

# read in site data
sample_coords<-read.csv("refs/CCGP-trout_sample-coordinate-match.csv", stringsAsFactors = FALSE)
sample_coords <- subset(sample_coords,!is.na(longitude))
# make this a spatial data frame 
# (i.e., spatial data needs to be in numeric form (use UTM, lat/long, or X/Y columns, doesn't matter you can convert later) 
# using lat/long
# need to figure out what the geodetic datum (or reference for representing the positions of locations on Earth)
# these are the CRS datum/projections "strings"
# in this case, will use the default for lat/long
# using the sf package instead of sp
# https://stackoverflow.com/questions/68157684/set-crs-for-latitude-longitude-point-data
sites.sf <- st_as_sf(sample_coords, coords=c("longitude", "latitude"), crs="EPSG:4326")
# geodetic CRS: WGS 84 (EPSG:4326)
st_crs(sites.sf)

# convert sf objects to spatial object for downstream mapmaking
hist_trout_ranges.sp <- st_zm(hist_trout_ranges)
hist_trout_ranges.sp <- as(hist_trout_ranges.sp,"Spatial")
sites.sp <- as(sites.sf,"Spatial")

# pull out the coastal rainbow/steelhead range so that it can be laid down first
hist_trout_ranges.steelhead <- subset(hist_trout_ranges, Species=="Coastal Rainbow Trout")
# pull out all but steelhead
hist_trout_ranges.nonsteelhead <- subset(hist_trout_ranges, Species!="Coastal Rainbow Trout")


## make map
# rainbow/cutthroat colorscheme
# rainbow more blue/purple, cutthroat pink/brown colors, clarkii purple
# save salmon for outgroup
# "Warner Lake Redband Trout" 
# "Bull Trout (extirpated)"   
# "Coastal Cutthroat Trout"   
# "Kern River Rainbow Trout"  
# "Goose Lake Redband Trout" 
# "California Golden Trout"   
# "McCloud Redband Trout"     
# "Little Kern Golden Trout"  
# "Lahontan Cutthroat Trout"  
# "Eagle Lake Rainbow Trout" 
# "Paiute Cutthroat Trout"  
# "Coastal Rainbow Trout"
# extra colors: 
"#959ABF"
# slim palette: steelhead, redband, mccloud, kern, coastal cutts, lahontan, paiute, salmon
trout_palette <- c("#7794A1","#546A7C", "#BBD3D1", "#6D7046", "#BFBDD1", "#D3928C", "#8E7068","#F4AD8B")
hist_range_map_palette <- c("#0082BA","#FBFFBA","#BFBDD1","#6D7046","#ACD6E8","#467047","#BBD3D1","#95BF96","#D3928C","#253333","#8E7068")
hist_range_map_palette<-c("Coastal Rainbow Trout"="#7794A1","Warner Lake Redband Trout"="#0082BA","Bull Trout (extirpated)"="#FBFFBA",
                "Coastal Cutthroat Trout"="#BFBDD1","Kern River Rainbow Trout"="#6D7046","Goose Lake Redband Trout"="#ACD6E8",
                "California Golden Trout"="#467047","McCloud Redband Trout"="#BBD3D1","Little Kern Golden Trout"="#95BF96",
                "Lahontan Cutthroat Trout"="#D3928C","Eagle Lake Rainbow Trout"="#253333","Paiute Cutthroat Trout"="#8E7068")
counties <- map_data("county")
ca_county <- subset(counties, region == "california")
ca_base <- ggplot(data = ca_county, mapping = aes(x = long, y = lat, group = group)) + theme_bw() +
  geom_polygon(color = "white", fill = "#ECF1ED") + 
  geom_sf(data = hist_trout_ranges.steelhead, aes(color=Species, fill=Species),
        inherit.aes = FALSE, show.legend="polygon") +
#        color = "#7794A1", fill=adjustcolor("#7794A1", alpha.f=1)) +
  geom_sf(data = hist_trout_ranges.nonsteelhead, aes(color=Species, fill=Species),
          inherit.aes = FALSE, show.legend="polygon") +
#          color = hist_range_map_palette, fill=adjustcolor( hist_range_map_palette, alpha.f=1)) +
#  geom_path(data=ca_rivers.sf_Rivers_coords.df, aes(X, Y, group=L1, fill=NULL), 
#            color="white", alpha=0.7) + 
  scale_fill_manual(values=hist_range_map_palette) +
  scale_color_manual(values=hist_range_map_palette) +
  theme(legend.position = c(.85, .81)) 

# could also color points by species
hist_trout_range_map <- 
  ca_base +
  geom_point(data=sample_coords, aes(x=longitude, y=latitude, group=geo_group)) +
  geom_label(data=sample_coords, aes(x=longitude, y=latitude, group=geo_group, label=sample_ID), 
             nudge_x=0.3, nudge_y=0, label.size=0.05, size=2,  
             fontface = "bold.italic", label.r=unit(0.20, "lines")) 

ggsave(ca_base, filename='Figures/historical-trout-range-map.pdf', device = cairo_pdf, width=10, height=14, bg="transparent")
ggsave(hist_trout_range_map, filename='Figures/historical-trout-range-map_ccgp-trout-sites_labelled.pdf', device = cairo_pdf, width=10, height=14, bg="transparent")



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



## CA rivers and creeks
ca_rivers <- st_read("NHD_Major_Rivers_and_Creeks/NHD_Major_Rivers_and_Creeks.shp")
ca_rivers_nozm <- st_zm(ca_rivers) # 324886 features and 16 fields
# can show specific features by subsetting them like this
print(ca_rivers_nozm[9:15], n = 3)
ggplot(ca_rivers_nozm) + geom_sf()

## Historic trout ranges
historic_trout_sfdf <- st_read("Historic_Trout_Watersheds_[ds440]/Historic_Trout_Watersheds_[ds440].shp")
plot(st_geometry(historic_trout_sfdf))
plot(historic_trout_sfdf["Species"])
ggplot() + geom_sf(historic_trout_sfdf, aes(fill = "Species")) 


## Map with mapview
## to choose map.types: https://leaflet-extras.github.io/leaflet-providers/preview/

library(mapview)

clarkii_coords <- read.csv("metadatas/Oclarkii_sample_coordinates.csv", header=T)
mykiss_coords <- read.csv("metadatas/Omykiss_sample_coordinates.csv", header=T)
trout_coordinates <- rbind(clarkii_coords,mykiss_coords)


# the base is nice
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", crs = 4269, grid = FALSE)
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="Stadia.OSMBright", crs = 4269, grid = FALSE)
# I think this is my favorite
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="Stadia.Outdoors", crs = 4269, grid = FALSE)
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="Stadia.StamenTerrain", crs = 4269, grid = FALSE)
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="Esri.NatGeoWorldMap", crs = 4269, grid = FALSE)
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="Esri.WorldTopoMap", crs = 4269, grid = FALSE)
# this one's good
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="CartoDB.VoyagerLabelsUnder", crs = 4269, grid = FALSE)
mapview(trout_coordinates, zcol="subspecies", xcol = "latitude", ycol = "longitude", map.types="USGS.USImageryTopo", crs = 4269, grid = FALSE)


## this is a fabulous tutorial for making cool california river maps!
## https://ryanpeek.org/2016-09-28-static_maps_in_r/

