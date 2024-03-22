# Read this shape file with the rgdal library. 
library(sf)
library(ggplot2)

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