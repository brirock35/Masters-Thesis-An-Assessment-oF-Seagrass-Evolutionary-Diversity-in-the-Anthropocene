
# Seagrass Biogeography and Global Change

The purpose of this study is to use machine learning approaches to model the distributions of seagrass species under various scenarios of climate change and combine these predictions with a dated molecular phylogeny to predict the impact of climate change on seagrass evolutionary diversity.

### Authors

[Brianna M. Rock](https://github.com/brirock35)

[Barnabas H. Daru](https://barnabasdaru.com/) 

#### Texas A&M University - Corpus Christi, MARB Program

Daru Lab

![alt text](https://barnabasdaru.files.wordpress.com/2018/05/title-picture.png?w=300&h=300g)

## Contents

1. [Introduction](#introduction)
1. [Methodology](#methodology)
1. [Results](#results)
1. [Discussion](#discussion)
1. [Conclusion](#conclusion)
1. [References](#references)
1. [Acknowledgments](#acknowledgments)

## Introduction

Seagrass meadows provide enormous ecosystem goods and services, ultimately establishing complex coastal habitats that teem with marine diversity. However, the future of seagrass communities under rapid shifts in climate in the Anthropocene—a period of most pronounced human impact on biotic communities—has recently been questioned. Seagrasses face elevated extinction pressures from environmental change and to a greater extent anthropogenic threats such as eutrophication, nutrient sedimentation, marine dredging, overfishing, and irresponsible recreational boating practices. Non-random extinction patterns coupled with human-induced climate change has the potential to negatively affect entire clades of seagrass evolutionary history, which we suggest is a more defining feature of biotic change in the Anthropocene. Over time, such rapid environmental changes may lead to shifts in the distributions of these marine plant taxa. Concurrent to potential range shifts in seagrasses is the threat of complete extinctions of specific species.

A large body of existing research supports the idea that current ecological processes, e.g. competition, grazing, light and nutrient availability, temperature, and physical exposure, shape seagrass distribution [@greve2004factors]. However, the significance of seagrass biogeography and phylogenetic structure as important components for biodiversity remain largely unexplored. This study will utilize carefully curated occurrence data to generate species distribution models under various scenarios of climate change. We will then generate new DNA sequences from fresh and herbarium material to construct a dated molecular phylogeny of extant seagrass species to explore the impact of global change on the evolutionary diversity of seagrasses.

Ultimately, this study will contribute to our understanding of how seagrass diversity will be impacted under global environmental change. The dated molecular phylogeny generated in this study will contribute to resolving taxonomic relationships among taxa while at the same time providing insights into climate change impacts on shared evolutionary history. As ecosystem engineers, the loss of seagrass species would be detrimental to the overall diversity of marine coastal systems. Therefore, the findings of this research can inform conservation management how best to maintain seagrass diversity in an uncertain future.

![image](https://www.earthisland.org/images/uploads/seagrass_crop.jpg) Photograph by: Christian Gloor (https://www.flickr.com/photos/christian_gloor/30774220094/).

## Methodology
### I. Species Distribution Modeling 
#### Seagrass Data Compiling

For this study, we compiled spatial data records for all seagrass species from four separate open-source databases: the Global Biodiversity Information Facility (https://www.gbif.org/), iDigBio (https://www.idigbio.org/), SEINet (http://swbiodiversity.org/seinet/), and the Seagrass Watch (www.seagrasswatch.org/). 

From these raw datasets, we selected specific variables that would later be useful for our analyses. Below is an example of how this was accomplished using one of our downloaded datasets:

```{r}
#Below is example of editing the SEINet raw data

#Seinet
library(data.table)
d <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Seagrasses/Seinet/SymbOutput_2020-06-08_095857_DwC-A/occurrences.csv")`

d1 <- d[,c("basisOfRecord", "family", "scientificName", "year", "month", "day", "habitat", "verbatimAttributes", "decimalLongitude", "decimalLatitude", "locality", "country",  "county", "stateProvince")]
names(d1) <- c("basis", "family", "species", "year", "month", "day", "habitat", "description", "lon", "lat", "locality", "country", "county", "state")

d1$species <- gsub("_"," ", d1$species)
```

An additional column was added into the dataset labeled “source” to depict the source of the spatial data records.

```
d1$source <- "SEINET"
```

Once all columns were ordered, the spatial data from all four sources were combined into one comma delimited file (.csv) together using the “rbind” function. Each dataset contributed the same categories of information.

```
# Combine SEINet data with GBIF data
g <- fread("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Seagrasses/GBIF/seagrass occurrances jan_17 2020.csv", sep = "\t")
g1 <- g[,c("basisOfRecord", "family", "species", "year", "month", "day", "decimalLongitude", "decimalLatitude", "locality", "stateProvince")]
names(g1) <- c("basis", "family", "species", "year", "month", "day", "lon", "lat", "locality", "state")

g1$source <- "GBIF"
z <- rbind(d1, g1)
```
Finally, the newly combined dataset were written to the desired folder. This process was repeated for each dataset until we had one master dataset possessing all spatial information for each seagrass species.

``` 
write.csv(z, "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/raw_occurences/master_occrrences.csv", row.names= F)
``` 

#### Climate Data 
The necessary climate layers were acquired from the Bio- ORACLE (http://www.bio-oracle.org/index.php) database. We downloaded benthic and surface layers for the current year (2020), and two future decades (2040-2050 and 2090-2100) for the following variables (mean and range): temperature, salinity, and currents velocity. These datasets were also gathered under four representative concentration pathway scenarios (RCP26, RCP45, RCP60, RCP85) that depict various green house gas emission levels and concentrations.


#### Georeferencing the Master Data
Before georeferencing the data, we took the master data and converted the spatial points into shapefiles for each of the 68 included species. We then utilized QGIS (https://qgis.org/en/site/) to manually prune, or remove, all obviously incorrect points through comparing the location of the spatial points with a shapefile layer of polygons where the species were likely to be found based off of expert opinion and prior research. Next, we assessed which species had very little to no points, for which the associated raw data should be further consulted to preserve more spatial data. We also determined which species possessed large numbers of spatial points that should be kept, but fell outside of known polygon ranges. These points would also need to be georeferenced. In total, this gave us 20 total species to isolate for georeferencing. 

To begin the geolocation process, we first utilized the GEOLocate Web Client (https://www.geo-locate.org/web/default.html) to batch georeference the locality information for the 20 species in need of georeferencing to provide lat/lon coordinates. Once this was accomplished, we exported the results as comma delimited files for each species. We later converted these into shapefiles containing the spatial points (coordinates) obtained from the geolocation process.

```
# The below example shows how we converted a .csv file with geolocated points for a seagrass (Zostera chilensis) unto a shapefile for that species

library(phyloregion)
fx <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/GeoLocated/GeoLocated_points_added/Zostera chilensis.csv")
fx <- fx[,c("Lat", "Lon")]
fx <- na.omit(fx)
fx <- as.data.frame(fx)
fx$species <- "Zostera chilensis"
fx <- fx[,c(3,1,2)]
coordinates(fx) <- ~Lon + Lat

rgdal::writeOGR(fx, dsn = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/GeoLocated/GeoLocated_points_added/shapefiles_GeoLocate_final", layer = "Zostera_chilensis", driver = "ESRI Shapefile", overwrite_layer = TRUE)
```
Next, we obtained a Google API and created a workflow to geocode these 20 species again using the same original files with the locality information on each species. The outputs of this process were saved as shapefiles consisting of the spatial points (coordinates) for each species that had been geocoded.

```
library(ggmap)
library(tmaptools)
library(RCurl)
library(jsonlite)
library(tidyverse)
library(leaflet)

register_google(key = "yourAPIhere")

fx <- list.files(path = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points", pattern = ".csv")
fx <- gsub(".csv", "", fx)
f <- list.files(path = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_to_more_points", pattern = ".csv", full.names = TRUE)
names(f) <- fx

for (i in seq_along(f)) {
  r <- read.csv(f[[i]])
  m <- geocode(location = as.character(r$LocalityString), source = "google") 
  res <- cbind(species = labels(f[i]), m, locality = r[, c("LocalityString", "Country",  "State", "County")])
  res <- res[!is.na(res$lon), ]
  coordinates(res) <- ~lon + lat
  rgdal::writeOGR(res, dsn = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/geocode_R", layer = labels(f[i]), driver = "ESRI Shapefile", overwrite_layer = TRUE)
}
```

We then went through another round of QGIS point pruning with both the outputs of the GEOLocate site and the outputs of the google geocode R function. These points were validated by checking the locality information listed in the associated attribute tables and were layered with the originally pruned shapefiles for additional comparison. Once the shapefiles were pruned and the points validated, the final shapefiles were once again converted to comma delimited files and combined with the rest of the originally pruned species (48 species) and saved to the master data .csv file. 

```
library(raster)
# Converting geocoded species to a .csv file

#geocode
f1 <- list.files(path = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/GeoLocated/shapefiles_GeoLocate_final", pattern = ".shp", full.names = TRUE)

out <- NULL
for (i in seq_along(f1)) {
  tryCatch({
    m <- as.data.frame(shapefile(f1[[i]]))
    m <- m[, c("species", "coords.x1", "coords.x2")]
    names(m)[c(2,3)] <- c("lon", "lat")
    m$source <- "GeoLocate"
    out <- rbind(out, m)
    print(i)
  }, error= function(e){})
  
}

write.csv(out, "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/GeoLocated/shapefiles_GeoLocate_final/cleaned_GeoLocate.csv", row.names = FALSE)
```


```
# Converting the pruned master dataset with newly georeferenced species

d1 <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Master_occurrances/pointpruned_cleaned_masterlist_7_14_2020.csv")
d1$source <- "originally_ok"

d2 <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/geocode_R/cleaned_geocode.csv")
d3 <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_more_points/GeoLocated/GeoLocated_points_added/shapefiles_GeoLocate_final/cleaned_GeoLocate.csv")
d4 <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_to_relocate_points/geocode_R/cleaned_geocode.csv")
d5 <- read.csv("/Users/darulab/Desktop/BriannaR/Research/SDMs/data/georeference/need_to_relocate_points/GeoLocated/shapefiles_GeoLocate_final/cleaned_GeoLocate.csv")

m <- rbind(d1, d2, d3, d4, d5)
mx <- unique(m)

write.csv(mx, "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Master_occurrances/complete_cleaned_masterdata.csv", row.names = FALSE)
```


#### Preprocessing the Master Data
In order to run the data within the modeling code, we first need to convert the spatial data (geocoordinates) into community data. This was accomplished using a combination of the “data.table”, “phyloregion”, and “raster” packages.  We selected a climate raster that we had previously rescaled to create raw grids (created only from the rescaled climate raster layer) and also created pruned grids (using the poly shapefiles created from the points2comm function with the species data) to use for producing our spatial models. The community data was produced from our master dataset using the "points2comm" function. The resulting community data was saved as a final comma delimited file that will be used in the SDM analyses.

```
library(biomod2)
library(data.table)
library(phyloregion)
library(raster)

p <- shapefile("/Users/darulab/Desktop/BriannaR/Research/SDMs/rasters/clim_shp/resampled_clim_shp.shp")
p <- p[!is.na(p@data$bnthc__), ]
p$grids <- paste0("v", 1:nrow(p))
pp <- p
pp <- coordinates(pp)
p1 <- cbind(p, pp)
p1 <- p1[, c("grids", "X1", "X2")]
names(p1) <- c("grids", "lon", "lat")

# Get community data
m <- points2comm(d, lon = "lon", lat = "lat", shp.grids = p1)
plot_swatch(m$poly_shp, m$poly_shp$richness, border = NA)
com <- sparse2long(m$comm_dat)
write.csv(com,"/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Master_occurrances/master_final/resampled_comm_final_data.csv",  row.names= F)

rgdal::writeOGR(m$poly_shp, dsn = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Master_occurrances/master_final/shapefiles", layer = "resampled_pruned_grids", driver = "ESRI Shapefile")

rgdal::writeOGR(p1, dsn = "/Users/darulab/Desktop/BriannaR/Research/SDMs/data/Master_occurrances/master_final/shapefiles", layer = "resampled_raw_grids", driver = "ESRI Shapefile")
```

### Modeling 
To build our species distribution models for each seagrass species, we first read in the necessary packages. The primary package that we will be using for the species modeling is "phyloregion" (https://cran.r-project.org/web/packages/phyloregion/index.html).  Afterwards, we also read in our preprocessed datasets for both the climate data as well as the seagrass spatial data and grid cells. We also ensure that our working directory is set to the location where we want the results written to. 

#### Modeling the Master Data Using "phyloregion"



### Analyses 

#### Alpha Diversity
Changes in α-diversity was determined by separately computing species richness (SR), weighted endemism (WE), phylogenetic diversity (PD), and phylogenetic endemism (PE) between current and future climate scenarios. Thus, shifts in -diversity is expressed as:

 =(_j  – _i)/_i 

where i is species composition under current climate and j is species composition under future scenarios. 

```
rm(list = ls()) 
library(phyloregion)
library(scico)
library(raster)

r <- raster("/Users/BriRock/Desktop/Seagrass Research/Data/worldRaster_SDM_50km.tif")
d <- read.csv("/Users/BriRock/Desktop/fixed_alpha_csvs/2100_RCP85.csv")

r[1:ncell(r)] <- paste0("v", seq_len(ncell(r)))
index <- match(values(r), d$grids)
z1 <- setValues(r, d$cur_sr[index])
z2 <- setValues(r, d$fut_sr[index])
z3 <- setValues(r, d$dx_sr[index])
z4 <- setValues(r, d$cur_pd[index])
z5 <- setValues(r, d$fut_pd[index])
z6 <- setValues(r, d$dx_pd[index])
z7 <- setValues(r, d$sesPDcur[index])
z8 <- setValues(r, d$sesPDfut[index])
z9 <- setValues(r, d$dx_sespd[index])

CLR <- scico(255, palette = 'batlow')

range(d$cur_sr)
range(d$fut_sr)
range(d$dx_sr)
range(d$cur_pd)
range(d$fut_pd)
range(d$dx_pd)
range(d$sesPDcur)
range(d$sesPDfut)
range(d$dx_sespd)

postscript("~/Desktop/2100_RCP85_PD_SR.ps", height = 8, width = 12)
par(mfrow=c(3,3))
par(mar=rep(1,4))
plot(z1,  col=CLR, axes=FALSE, legend=FALSE)
plot(z1, legend.only=TRUE, col=CLR, legend.width=1, legend.shrink=0.5,
     smallplot=c(0.05,0.07, 0.4,0.5))

plot(z2,  col=CLR, axes=FALSE, legend=FALSE)
plot(z2, legend.only=TRUE,  col=CLR, legend.width=1, legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z3, zlim=c(-0.875, 7.0000000), col=CLR, axes=FALSE, legend=FALSE)
plot(z3, legend.only=TRUE, zlim=c(-0.875, 7.0000000), col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z4, col=CLR, axes=FALSE, legend=FALSE)
plot(z4, legend.only=TRUE,  col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z5,  col=CLR, axes=FALSE, legend=FALSE)
plot(z5, legend.only=TRUE,  col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z6, col=CLR, axes=FALSE, legend=FALSE)
plot(z6, legend.only=TRUE, col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z7, zlim=c(-4.009199, 1.952685), col=CLR, axes=FALSE, legend=FALSE)
plot(z7, legend.only=TRUE, zlim=c(-4.009199, 1.952685), col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z8, zlim=c(-4.009199, 1.952685), col=CLR, axes=FALSE, legend=FALSE)
plot(z8, legend.only=TRUE, zlim=c(-4.009199, 1.952685), col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

plot(z9, zlim=c(-19506.44, 3442957.67), col=CLR, axes=FALSE, legend=FALSE)
plot(z9, legend.only=TRUE, zlim=c(-19506.44, 3442957.67), col=CLR, legend.width=1, 
     legend.shrink=0.75, smallplot=c(0.05,0.07, 0.4,0.5))

dev.off()


```

#### Beta Diversity



