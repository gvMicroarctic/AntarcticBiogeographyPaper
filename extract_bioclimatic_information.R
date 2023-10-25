################################################################
###Retrieve data from CHELSA for sub-Antractic/Antarctic soil###
################################################################

#set working directory
setwd("/home/gilda/all")

library(terra)

#for reference, here is the link to the CHELSA files: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2F

#import raster and create RasterStack
rasStack_bio1 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio1_1981-2010_V.2.1.tif")
rasStack_bio2 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio2_1981-2010_V.2.1.tif")
rasStack_bio4 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio4_1981-2010_V.2.1.tif")
rasStack_bio5 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio5_1981-2010_V.2.1.tif")
rasStack_bio10 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio10_1981-2010_V.2.1.tif")
rasStack_bio12 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio12_1981-2010_V.2.1.tif")
rasStack_bio14 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio14_1981-2010_V.2.1.tif")
rasStack_bio15 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio15_1981-2010_V.2.1.tif")
rasStack_bio17 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio17_1981-2010_V.2.1.tif")
rasStack_bio18 <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_bio18_1981-2010_V.2.1.tif")
rasStack_scd <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_scd_1981-2010_V.2.1.tif")
rasStack_swe <- rast("/home/gilda/all/extract_geochemical_data/db5/CHELSA_swe_1981-2010_V.2.1.tif")

#import geographical points
pointCoordinates <- data.frame(read.table("/home/gilda/all/all_coordinates.txt"))
dim(pointCoordinates)
#[1] 990   2
#coordinates(pointCoordinates)= ~ longitude+ latitude

#extract vales from raster for the geographical points
rasValue_bio1 = extract(rasStack_bio1, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio2 = extract(rasStack_bio2, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio4 = extract(rasStack_bio4, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio5 = extract(rasStack_bio5, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio10 = extract(rasStack_bio10, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio12 = extract(rasStack_bio12, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio14 = extract(rasStack_bio14, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio15 = extract(rasStack_bio15, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio17 = extract(rasStack_bio17, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_bio18 = extract(rasStack_bio18, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_scd = extract(rasStack_scd, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude
rasValue_swe = extract(rasStack_swe, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude

#save data in csv file
combinePointValue=cbind(sample=row.names(pointCoordinates),pointCoordinates, rasValue_bio1,rasValue_bio2,rasValue_bio4,rasValue_bio5,
                        rasValue_bio10,rasValue_bio12,rasValue_bio14,rasValue_bio15,
                        rasValue_bio17,rasValue_bio18,rasValue_scd,rasValue_swe)
dim(combinePointValue)
#[1] 990  15

###################################################################
###Retrieve data from bioregion for sub-Antarctic/Antarctic soil###
###################################################################

#for reference, here is the link to the CHELSA files: https://envicloud.wsl.ch/#/?prefix=chelsa%2Fchelsa_V2%2FGLOBAL%2Fclimatologies%2F

#import raster and create RasterStack
rasStack <- rast("/home/gilda/all/extract_geochemical_data/db3/ACBR_WGS84_30s.tif")

#extract vales from raster for the geographical points
bioregion = extract(rasStack, pointCoordinates[,c(2,1)])[,2] #column 2 reports longitude and column 1 reports latitude

#save data in csv file
combinePointValue_new = cbind(combinePointValue, bioregion)
dim(combinePointValue_new)
#[1] 990  16

#add dataset information
load("all_dataset_ps_srs_bacteria.RData")
library(phyloseq)

row.names(sample_data(ps_srs)) == combinePointValue$sample #all TRUE

#save data in csv file
combinePointValue_new = cbind(combinePointValue, bioregion, dataset = sample_data(ps_srs)$dataset)
dim(combinePointValue_new)
#[1] 990  17

write.csv(combinePointValue_new,"bioclim_bioregion_antarctic_study.csv", row.names = FALSE)

#bioclim_bioregion_antarctic_study_curated.csv : NA substituted with bioregions and curated in ArcMap

#############################
###Retrieve elevation data###
#############################

#download elevation data from WorldClim
#https://www.worldclim.org/data/worldclim21.html

#import geographical points
coordinates <- read.csv("bioclim_bioregion_antarctic_study_curated.csv")

#import raster and create RasterStack
rasStack_ele <- rast("/home/gilda/all/extract_geochemical_data/db1/wc2.1_30s_elev/wc2.1_30s_elev.tif")

#extract vales from raster for the geographical points
rasValue_ele = extract(rasStack_ele, coordinates[,c(3,2)])[,2] #column 3 reports longitude and column 2 reports latitude

#N.B. This method does not really work. So I am getting data from https://www.gpsvisualizer.com/elevation; just need to upload latitude and longitude

######################################
###Retrieve distance from coastline###
######################################

#https://github.com/mdsumner/distancetocoast
#https://stackoverflow.com/questions/35555709/global-raster-of-geographic-distances

#devtools::install_github("mdsumner/distancetocoast")
library(distancetocoast)
library(terra)

#prepare longitude, latitude table
long_lat <- data.frame(coordinates[,3], coordinates[,2]) #all longitudes + latitudes
colnames(long_lat) <- c("long", "lat")

#extract values
rasValue_coast <- terra::extract(distance_to_coastline_lowres, cbind(long_lat))
#many NA; I will turn NA to 0 because those are the samples close to coastline

#Integrate new data
pointCoordinates_curated <- cbind(coordinates, "elevation" = rasValue_ele, "dist_coast" = rasValue_coast)
#not perfect especially for islands so I am only use distance to coast information for samples on continent.

#Write new file
write.csv(pointCoordinates_curated,"bioclim_bioregion_antarctic_study_curated_new.csv", row.names = FALSE)
