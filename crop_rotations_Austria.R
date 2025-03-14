# Please be aware: this project runs on a HPC unit at ZALF.
# Your local machine might not allow the processing of large data sets.

# We are providing code for Austria as a sample, where all necessary data is available open-source.
# As the manuscript is based on the entire study domain, results using this script will differ
# from the ones reported in the manuscript.

#----------------------------------------------------Prep work-----------------------------------------------------------

#### Set working directory and load required packages ####

setwd("./crop_rotations_central_europe") #set accoding to your working directory#

require(abind)
require(caret)
require(CAST)
require(data.table)
require(dbplyr)
require(dplyr)
require(dtplyr)
require(ggplot2)
require(ggsci)
require(ncdf4)
require(purrr)
require(R.utils)
require(ranger)
require(readxl)
require(reshape2)
require(rjson)
require(sf)
require(stringr)
require(terra)
require(tidymodels)
require(tidyverse)



#### Colour setting for display ####

npg_colors_full<-pal_npg("nrc")(10)
npg_colors_light<-pal_npg("nrc", alpha = 0.6)(10)

low_color<-"white"
train_color<-npg_colors_full[1]
test_color<-npg_colors_full[3]

crop_color<-c(npg_colors_full[1:7], npg_colors_light[10], npg_colors_full[8:9], npg_colors_light[1:7])


#### Build custom functions  ####

add_P1CType <- function(dat){
  P1CType <- dat$CType
  P1CType[seq(0, nrow(dat), by=n_years)] <- NA
  P1CType <- c(NA, P1CType[-length(P1CType)])
  return(P1CType)
}

add_P2CType <- function(dat){
  P2CType <- dat$CType
  P2CType[sort(c(seq(0, nrow(dat), by=n_years),
                 seq(-1, nrow(dat)-1, by=n_years)))[-c(1:2)]] <- NA
  P2CType <- c(NA, NA, P2CType[-c(length(P2CType)-1,
                                  length(P2CType))])
  return(P2CType)
}

add_P3CType <- function(dat){
  P3CType <- dat$CType
  P3CType[sort(c(seq(0, nrow(dat), by=n_years),
                 seq(-1, nrow(dat)-1, by=n_years),
                 seq(-2, nrow(dat)-2, by=n_years)))[-c(1:3)]] <- NA
  P3CType <- c(NA, NA, NA, P3CType[-c(length(P3CType)-2,
                                      length(P3CType)-1,
                                      length(P3CType))])
  return(P3CType)
}

add_P4CType <- function(dat){
  P4CType <- dat$CType
  P4CType[sort(c(seq(0, nrow(dat), by=n_years),
                 seq(-1, nrow(dat)-1, by=n_years),
                 seq(-2, nrow(dat)-2, by=n_years),
                 seq(-3, nrow(dat)-3, by=n_years)))[-c(1:4)]] <- NA
  P4CType <- c(NA, NA, NA, NA, P4CType[-c(length(P4CType)-3,
                                          length(P4CType)-2,
                                          length(P4CType)-1,
                                          length(P4CType))])
  return(P4CType)
}

read_climate <- function(data_ref, type) {
  data_clim_raw <- get(tolower(type))
  clim <- data.frame("cell_id"=data_clim_raw$cell_id)
  time_scale = c("15","16","17","18","19","20","21","22","23")
  
  j <- 4-12
  
  for (i in time_scale){
    h<-j-12
    g<-h-12
    f<-g-12
    e<-f-12
    clim[paste(i,"tAVG",sep="")] <- data.frame(rowMeans(subset(data_clim_raw,select=c(
      if_else(e>0,e,4):if_else((e+5)>0,e+5,4+5),
      if_else(f>0,f,4):if_else((f+5)>0,f+5,4+5),
      if_else(g>0,g,4):if_else((g+5)>0,g+5,4+5),
      if_else(h>0,h,4):if_else((h+5)>0,h+5,4+5),
      if_else(j>0,j,4):if_else((j+5)>0,j+5,4+5))))) # March to August from the last 5 years #
    clim[paste(i,"tAM",sep="")] <- data.frame(rowMeans(subset(data_clim_raw,select=c(
      if_else((e+1)>0,e,4+1):if_else((e+2)>0,e+2,4+2),
      if_else((f+1)>0,f,4+1):if_else((f+2)>0,f+2,4+2),
      if_else((g+1)>0,g,4+1):if_else((g+2)>0,g+2,4+2),
      if_else((h+1)>0,h,4+1):if_else((h+2)>0,h+2,4+2),
      if_else((j+1)>0,j,4+1):if_else((j+2)>0,j+2,4+2))))) # April to May from the last 5 years #
    clim[paste(i,"tJJ",sep="")] <- data.frame(rowMeans(subset(data_clim_raw,select=c(
      if_else((e+3)>0,e,4+3):if_else((e+4)>0,e+4,4+4),
      if_else((f+3)>0,f,4+3):if_else((f+4)>0,f+4,4+4),
      if_else((g+3)>0,g,4+3):if_else((g+4)>0,g+4,4+4),
      if_else((h+3)>0,h,4+3):if_else((h+4)>0,h+4,4+4),
      if_else((j+3)>0,j,4+3):if_else((j+4)>0,j+4,4+4))))) # June to July from the last 5 years #
    j <- j+12 # jump to the next year #
  }
  
  rm(data_clim_raw)
  gc()
  
  data_clim_ <- left_join(data_ref_clean, clim, by = c("cell_id")) %>% 
    as.data.table()
  
  rm(clim, i, j, e, f, g, h)
  gc() 
  
  data_clim <- data_clim_ %>% 
    select(cell_id, OBJECTID, ends_with('AVG')) %>% 
    pivot_longer(!c("cell_id", "OBJECTID"), names_to = "Year", values_to = paste(type,"AVG", sep="")) %>% 
    mutate(Year = gsub("^(.{2})(.*)$", "20\\1", Year)) %>% 
    as.data.table()
  data_clim <- data_clim_ %>% 
    select(cell_id, OBJECTID, ends_with('JJ')) %>% 
    pivot_longer(!c("cell_id", "OBJECTID"), names_to = "Year", values_to = paste(type,"JJ", sep="")) %>% 
    mutate(Year = gsub("^(.{2})(.*)$", "20\\1", Year)) %>% 
    inner_join(data_clim, by=c("cell_id", "OBJECTID", "Year")) %>% 
    as.data.table() 
  data_clim <- data_clim_ %>% 
    select(cell_id, OBJECTID, ends_with('AM')) %>% 
    pivot_longer(!c("cell_id", "OBJECTID"), names_to = "Year", values_to = paste(type,"AM", sep="")) %>% 
    mutate(Year = gsub("^(.{2})(.*)$", "20\\1", Year)) %>% 
    inner_join(data_clim, by=c("cell_id", "OBJECTID", "Year")) %>% 
    as.data.table() 
  
  return(data_clim)
}

read_price <- function() { 
  data_price_raw <- fread("./data/prices_Austria.csv", sep=",", header = T)
  
  data_p <- data_price_raw %>% 
    mutate(Item = recode(Item, "Vegetables, leguminous nes" = "Vegetables, leguminous")) %>% 
    filter(Unit == 'USD') %>%
    as.data.table()
  
  data_p <- data_p %>% 
    pivot_wider(id_cols = Item, values_from = Value,names_from = Year) %>% 
    as.data.table()
  
  k <- which(is.na(data_p), arr.ind=TRUE)
  data_p[k] <- rowMeans(data_p[,-1], na.rm=TRUE)[k[,1]]
  
  
  for (n in c(2:20)) {
    data_p[Item=="Lupins"|Item=="Soya beans"][[n]]<-mean(data_p[Item=="Lupins"|Item=="Soya beans"][[n]])
    data_p[Item=="Cabbages"|Item=="Cauliflowers and broccoli"][[n]]<-mean(data_p[Item=="Cabbages"|Item=="Cauliflowers and broccoli"][[n]])
  }
  
  data_p<-data_p[-4,]
  
  data_p$Item[2]<-"Brassica"
  data_p$Item[3]<-"Carrots"
  data_p$Item[4]<-"Maize" 
  data_p$Item[6]<-"Onions"
  data_p$Item[8]<-"Canola"
  data_p$Item[10]<-"Legumes"
  data_p$Item[12]<-"Sunflower"
  
  scheme <- fread("./data/crop_classification.csv", sep=";", header = TRUE)
  scheme$`Price Description` <- tolower(scheme$`Price Description`)
  data_p$Item <- tolower(data_p$Item)
  
  temp_price <- scheme %>% 
    select(`Crop Class ID`, `Price Description`) %>% 
    as.data.table()
  colnames(temp_price) <- c('ID', 'Item')
  temp_price[7,2] <- "triticale"
  
  data_price <- inner_join(temp_price, data_p, by = "Item") %>% 
    as.data.table()
}


#----------------------------------------------------Create data base----------------------------------------------------

#### Harmonise IACS data and assign CType ####

crop_id<-read.csv("./data/crop_description_Austria.csv", sep=";") %>%
  select(SNAR_CODE, crop_class) %>%
  arrange(SNAR_CODE)

x<-which(duplicated(crop_id[,1]))
crop_id<-crop_id[-x,] %>%
  drop_na()

colnames(crop_id)<-c("SNAR_CODE", "CType")
crop_id$SNAR_CODE<-as.numeric(crop_id$SNAR_CODE)

filelist <- list.files(path="./data/IACS", pattern=".gpkg", full.names=TRUE)
# Original data is published by Agrarmarkt Austria (https://www.ama.at/) under licence CC BY 4.0 and was downloaded from 
# https://www.data.gv.at/suche/?katFilter%5B0%5D=httppublicationseuropaeuresourceauthoritydata-themeagri&searchterm&typeFilter%5B0%5D=dataset&nr=1&tagFilter%5B0%5D=INVEKOS

myYear<-seq(2015,2023, 1)
myInvecos<-list()

for (f in 1:length(filelist)) {
  data<-st_read(filelist[f])%>%
    select(SNAR_CODE, geom) %>%
    drop_na() %>%
    mutate(state="AT", Year=myYear[f])
  
  data$SNAR_CODE<-as.numeric(data$SNAR_CODE)
  
  data<-left_join(data, crop_id, by="SNAR_CODE")
  data<-data %>%
    drop_na() %>%
    select(Year, CType, state, geom) %>%
    filter(CType<18) # excluding minority crops #
  
  myInvecos[[f]]<-data
}

n_years<-length(filelist)

rm(data, crop_id, x, f, filelist, myYear)
gc()


#### Compile weather data ####

centroids <- st_centroid(myInvecos[[length(myInvecos)]]) %>% # taking the most recent year as a reference #
  arrange(geom)
x<-which(duplicated(centroids$geom))
centroids<-centroids[-x,]
centroids<-st_transform(centroids, crs="EPSG: 4258") # crs of the weather data #
centroids$OBJECTID<-c(1:nrow(centroids)) # objectid for identification #
rm(x)

myParam<-c("TM", "RR")
myName<-c("temp", "prec")
myTime<-c(seq(201001, 201012,1),
          seq(201101, 201112,1),
          seq(201201, 201212,1),
          seq(201301, 201312,1),
          seq(201401, 201412,1),
          seq(201501, 201512,1),
          seq(201601, 201612,1),
          seq(201701, 201712,1),
          seq(201801, 201812,1),
          seq(201901, 201912,1),
          seq(202001, 202012,1),
          seq(202101, 202112,1),
          seq(202201, 202212,1),
          seq(202301, 202312,1))

myCells<-matrix(1:(305*584), 305, 584) # austrian grid #

for (p in myParam) {
  filelist<-list.files(path="./data/weather", pattern=paste(p, ".nc", sep=""), full.names=TRUE)
  # Original data is published by GeoSphere Austria (https://www.geosphere.at/) under licence CC BY 4.0 and was downloaded from 
  # https://data.hub.geosphere.at/dataset/spartacus-v2-1d-1km

  myWeather<-list()
  for (f in 1:length(filelist)) {
    nc_data<-nc_open(filelist[f])
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat", verbose = F)
    myArray<- ncvar_get(nc_data, p)
    fillvalue <- ncatt_get(nc_data, p, "_FillValue")
    nc_close(nc_data) 
    myArray[myArray == fillvalue$value] <- NA
    
    newArray<-array(dim=c(305,584,12))
    for (a in 1:12) {
      newArray[,,a]<-apply(t(myArray[,,a]),2,rev)
    }
    
    myWeather[[f]]<-newArray
  }
  myWeather<-do.call(abind, myWeather)
  myWeather<-abind(myWeather, myCells)
  rm(myArray, newArray, nc_data, filelist, a)
  
  myBrick<-rast(myWeather)
  ext(myBrick) <- c(min(lon), max(lon), min(lat), max(lat))
  crs(myBrick) <- "EPSG:4258"
  
  data<-extract(myBrick, centroids)
  colnames(data)<-c("OBJECTID", paste(myName[which(myParam==p)], myTime, sep="_"), "cell_id")
  
  assign(paste0(myName[which(myParam == p)]), data)
}

rm(data, fillvalue, lat, lon, myBrick, myCells, f, myName, myParam, myTime, myWeather, p)
gc()


#### Create reference file for identification ####

data_ref_clean<-centroids %>% # create reference grid only containing identifiers #
  mutate(x_coord = sf::st_coordinates(.)[,1],
         y_coord = sf::st_coordinates(.)[,2],
         OBJECTID=prec$OBJECTID,
         cell_id=prec$cell_id)
data_ref_clean<-data_ref_clean %>%
  select(OBJECTID, state, x_coord, y_coord, cell_id) %>%
  st_drop_geometry() %>%
  as.data.table()

prec<-prec %>%
  select(!OBJECTID)%>%
  group_by(cell_id) %>%
  summarise(across(everything(), list(mean)))

temp<-temp %>%
  select(!OBJECTID)%>%
  group_by(cell_id) %>%
  summarise(across(everything(), list(mean)))


#### Summarise crop data ####

centroids<-st_transform(centroids, crs="EPSG: 31287") %>%
  select(OBJECTID, geom) %>%
  arrange(OBJECTID)

crop<-list()
for (i in 1:length(myInvecos)) {
  st_crs(myInvecos[[i]])<-"EPSG: 31287"
  myFile<-st_intersection(centroids, myInvecos[[i]])
  x<-which(duplicated(myFile$geom))
  myFile<-myFile[-x,]
  myFile<-myFile %>%
    mutate(x_coord = sf::st_coordinates(.)[,1],
           y_coord = sf::st_coordinates(.)[,2])
  myFile<-myFile %>%
    st_drop_geometry()
  crop[[i]]<-myFile
}

myCrop<-reduce(crop, left_join, by = c("OBJECTID", "x_coord", "y_coord", "state")) # check if that works #
myCrop<-myCrop %>%
  select(OBJECTID, CType.x, CType.y, CType.x.x, CType.y.y, CType.x.x.x, CType.y.y.y, CType.x.x.x.x, CType.y.y.y.y, CType)%>%
  replace(is.na(.), 0)
colnames(myCrop)<-c("OBJECTID", "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
myCrop<-myCrop %>%
  arrange(OBJECTID)

rm(x, i, myFile, myInvecos, crop)
gc()


#### Compile soil data ####

mySoil<-st_read("./data/soil/boka_1km_2016.shp")
# Original data is published by BFW (https://www.bfw.gv.at/) under licence CC BY 3.0 (replaced by CC BY 4.0) and was downloaded from 
# https://bodenkarte.at

mySoil$sand<-as.numeric(mySoil$sand)
mySoil$ton<-as.numeric(mySoil$ton)
mySoil$schluff<-as.numeric(mySoil$schluff)

mySoil<-mySoil %>%
  select(humus, sand, schluff, ton, geometry) %>% 
  mutate(Corg=humus/1.725, SAND=sand/100, silt=schluff/100, clay=ton/100)
mySoil<- mySoil %>%
  select(!c(humus, schluff, ton, sand)) %>%
  arrange(geometry)
colnames(mySoil)<-c("Corg", "sand", "silt", "clay", "geometry")

centroids<-st_transform(centroids, crs=st_crs(mySoil))

mySoil<-st_intersection(centroids, mySoil)
mySoil<-mySoil %>%
  st_drop_geometry() %>%
  select(OBJECTID, Corg, sand, silt, clay) %>%
  arrange(OBJECTID)

mySoil<-mySoil %>% 
  drop_na() %>%
  filter(Corg<(mean(na.omit(mySoil$Corg))+sd(na.omit(mySoil$Corg))^2)) %>% # eliminate Corg outliers #
  mutate_at(2:5, funs(round(., 2)))

### Get elevation and slope ###

myElev<-rast("./data/soil/elev.tif")
mySlope<-rast("./data/soil/slope.tif")
# Original data is published by geoland.at (http://www.geoland.at) under licence CC BY 4.0 and was downloaded from 
# https://www.data.gv.at/katalog/dataset/dgm

myTerrain<-c(myElev, mySlope)

centroids<-st_transform(centroids, crs=st_crs(myElev))

myTerrain<-extract(myTerrain, centroids) 
colnames(myTerrain)<-c("OBJECTID", "SElev", "SSlope")

mySoil<-left_join(mySoil, myTerrain, by=c("OBJECTID"))%>%
  mutate_at(2:6, funs(round(., 2)))

rm(myElev, mySlope, myTerrain)
gc()


#### Remove permanent grassland and minority crop fields ####

myCrop <- inner_join(data_ref_clean, myCrop, by="OBJECTID") %>% # merge the different identifiers.
  pivot_longer(-c("state", "OBJECTID", "cell_id", "x_coord", "y_coord"), names_to = "Year", values_to = "CType") %>% 
  arrange(OBJECTID) %>% 
  as.data.table()

data_arable<-myCrop %>% 
  group_by(OBJECTID) %>%
  filter(!all(CType==18)) %>% 
  as.data.table()

data_arable<-data_arable %>%
  group_by(OBJECTID) %>%
  filter(!all(CType %in% c(0,18,19,70,80))) %>%
  as.data.table()

rm(myCrop)
gc()


#### Add diversity features ####

data_rotation <- data_arable %>%
  group_by(OBJECTID) %>%
  summarize(
    t = sum((lag(CType[Year < 2023]) - CType[Year < 2023]) != 0, na.rm = TRUE),
    k = n_distinct(CType[Year < 2023]),
    b = sum(CType[Year < 2023] > n_years) / (n_years - 1),
    s = sum(CType[Year < 2023] %in% c(1:2, 7:9, 11:19)) / (n_years - 1)
  )

data_arable<-left_join(data_arable, data_rotation, b=c("OBJECTID"))

rm(data_rotation)
gc()


#### Add CType from previous 4 years ####

data_arable<-data_arable %>%
  arrange(OBJECTID, Year)

data_arable$P1CType <- add_P1CType(data_arable)
data_arable$P2CType <- add_P2CType(data_arable)
data_arable$P3CType <- add_P3CType(data_arable)
data_arable$P4CType <- add_P4CType(data_arable)


#### Remove years with grassland and minority crops ####

data_arable<-data_arable %>%
  drop_na() %>%
  group_by(OBJECTID) %>%
  filter(!(CType %in% c(0,18,19,70,80))) %>%
  filter(!(P1CType %in% c(0,18,19,70,80))) %>%
  filter(!(P2CType %in% c(0,18,19,70,80))) %>%
  filter(!(P3CType %in% c(0,18,19,70,80))) %>%
  filter(!(P4CType %in% c(0,18,19,70,80))) %>%
  as.data.table()


#### Add agronomic information ####

# Winter #
data_arable$WinterP1CType<-c("NO")
data_arable[P1CType==3|
              P1CType==4|
              P1CType==5|
              P1CType==6|
              P1CType==10]$WinterP1CType<-c("YES")

data_arable$WinterP2CType<-c("NO")
data_arable[P2CType==3|
              P2CType==4|
              P2CType==5|
              P2CType==6|
              P2CType==10]$WinterP2CType<-c("YES")

data_arable$WinterP3CType<-c("NO")
data_arable[P3CType==3|
              P3CType==4|
              P3CType==5|
              P3CType==6|
              P3CType==10]$WinterP3CType<-c("YES")

data_arable$WinterP4CType<-c("NO")
data_arable[P4CType==3|
              P4CType==4|
              P4CType==5|
              P4CType==6|
              P4CType==10]$WinterP4CType<-c("YES")

# Leaf #
data_arable$LeafP1CType<-c("YES")
data_arable[P1CType==1|
              P1CType==2|
              P1CType==3|
              P1CType==4|
              P1CType==5|
              P1CType==6|
              P1CType==7|
              P1CType==8|
              P1CType==9]$LeafP1CType<-c("NO")

data_arable$LeafP2CType<-c("YES")
data_arable[P2CType==1|
              P2CType==2|
              P2CType==3|
              P2CType==4|
              P2CType==5|
              P2CType==6|
              P2CType==7|
              P2CType==8|
              P2CType==9]$LeafP2CType<-c("NO")

data_arable$LeafP3CType<-c("YES")
data_arable[P3CType==1|
              P3CType==2|
              P3CType==3|
              P3CType==4|
              P3CType==5|
              P3CType==6|
              P3CType==7|
              P3CType==8|
              P3CType==9]$LeafP3CType<-c("NO")

data_arable$LeafP4CType<-c("YES")
data_arable[P4CType==1|
              P4CType==2|
              P4CType==3|
              P4CType==4|
              P4CType==5|
              P4CType==6|
              P4CType==7|
              P4CType==8|
              P4CType==9]$LeafP4CType<-c("NO")

# Draw #
data_arable$DrawP1CType<-c("MED")
data_arable[P1CType==13|
              P1CType==16]$DrawP1CType<-c("LOW")
data_arable[P1CType==1|
              P1CType==2|
              P1CType==10|
              P1CType==11|
              P1CType==12|
              P1CType==14|
              P1CType==15]$DrawP1CType<-c("HIGH")

data_arable$DrawP2CType<-c("MED")
data_arable[P2CType==13|
              P2CType==16]$DrawP2CType<-c("LOW")
data_arable[P2CType==1|
              P2CType==2|
              P2CType==10|
              P2CType==11|
              P2CType==12|
              P2CType==14|
              P2CType==15]$DrawP2CType<-c("HIGH")

data_arable$DrawP3CType<-c("MED")
data_arable[P3CType==13|
              P3CType==16]$DrawP3CType<-c("LOW")
data_arable[P3CType==1|
              P3CType==2|
              P3CType==10|
              P3CType==11|
              P3CType==12|
              P3CType==14|
              P3CType==15]$DrawP3CType<-c("HIGH")

data_arable$DrawP4CType<-c("MED")
data_arable[P4CType==13|
              P4CType==16]$DrawP4CType<-c("LOW")
data_arable[P4CType==1|
              P4CType==2|
              P4CType==10|
              P4CType==11|
              P4CType==12|
              P4CType==14|
              P4CType==15]$DrawP4CType<-c("HIGH")

# Organ #
data_arable$OrganP1CType<-c("GRAIN")
data_arable[P1CType==2|
              P1CType==11|
              P1CType==12|
              P1CType==16|
              P1CType==17]$OrganP1CType<-c("ROOT")
data_arable[P1CType==18]$OrganP1CType<-c("BIO")

data_arable$OrganP2CType<-c("GRAIN")
data_arable[P2CType==2|
              P2CType==11|
              P2CType==12|
              P2CType==16|
              P2CType==17]$OrganP2CType<-c("ROOT")
data_arable[P2CType==18]$OrganP2CType<-c("BIO")

data_arable$OrganP3CType<-c("GRAIN")
data_arable[P3CType==2|
              P3CType==11|
              P3CType==12|
              P3CType==16|
              P3CType==17]$OrganP3CType<-c("ROOT")
data_arable[P3CType==18]$OrganP3CType<-c("BIO")

data_arable$OrganP4CType<-c("GRAIN")
data_arable[P4CType==2|
              P4CType==11|
              P4CType==12|
              P4CType==16|
              P4CType==17]$OrganP4CType<-c("ROOT")
data_arable[P4CType==18]$OrganP4CType<-c("BIO")

# Water #
data_arable$WaterP1CType<-c("LOW")
data_arable[P1CType==3|
              P1CType==6|
              P1CType==10|
              P1CType==11|
              P1CType==13|
              P1CType==18]$WaterP1CType<-c("MED")
data_arable[P1CType==4|
              P1CType==5]$WaterP1CType<-c("HIGH")

data_arable$WaterP2CType<-c("LOW")
data_arable[P2CType==3|
              P2CType==6|
              P2CType==10|
              P2CType==11|
              P2CType==13|
              P2CType==18]$WaterP2CType<-c("MED")
data_arable[P2CType==4|
              P2CType==5]$WaterP2CType<-c("HIGH")

data_arable$WaterP3CType<-c("LOW")
data_arable[P3CType==3|
              P3CType==6|
              P3CType==10|
              P3CType==11|
              P3CType==13|
              P3CType==18]$WaterP3CType<-c("MED")
data_arable[P3CType==4|
              P3CType==5]$WaterP3CType<-c("HIGH")

data_arable$WaterP4CType<-c("LOW")
data_arable[P4CType==3|
              P4CType==6|
              P4CType==10|
              P4CType==11|
              P4CType==13|
              P4CType==18]$WaterP4CType<-c("MED")
data_arable[P4CType==4|
              P4CType==5]$WaterP4CType<-c("HIGH")

# Gap #
data_arable$GapP1CType<-0.75
data_arable[P1CType==18]$GapP1CType<-0
data_arable[P1CType==4|
              P1CType==5|
              P1CType==8]$GapP1CType<-0.33
data_arable[P1CType==3|
              P1CType==6|
              P1CType==7]$GapP1CType<-0.5
data_arable[P1CType==16|
              P1CType==17]$GapP1CType<-0.57
data_arable[P1CType==1|
              P1CType==2]$GapP1CType<-0.6
data_arable[P1CType==13]$GapP1CType<-0.77

data_arable$GapP2CType<-0.75
data_arable[P2CType==18]$GapP2CType<-0
data_arable[P2CType==4|
              P2CType==5|
              P2CType==8]$GapP2CType<-0.33
data_arable[P2CType==3|
              P2CType==6|
              P2CType==7]$GapP2CType<-0.5
data_arable[P2CType==16|
              P2CType==17]$GapP2CType<-0.57
data_arable[P2CType==1|
              P2CType==2]$GapP2CType<-0.6
data_arable[P2CType==13]$GapP2CType<-0.77

data_arable$GapP3CType<-0.75
data_arable[P3CType==18]$GapP3CType<-0
data_arable[P3CType==4|
              P3CType==5|
              P3CType==8]$GapP3CType<-0.33
data_arable[P3CType==3|
              P3CType==6|
              P3CType==7]$GapP3CType<-0.5
data_arable[P3CType==16|
              P3CType==17]$GapP3CType<-0.57
data_arable[P3CType==1|
              P3CType==2]$GapP3CType<-0.6
data_arable[P3CType==13]$GapP3CType<-0.77

data_arable$GapP4CType<-0.75
data_arable[P4CType==18]$GapP4CType<-0
data_arable[P4CType==4|
              P4CType==5|
              P4CType==8]$GapP4CType<-0.33
data_arable[P4CType==3|
              P4CType==6|
              P4CType==7]$GapP4CType<-0.5
data_arable[P4CType==16|
              P4CType==17]$GapP4CType<-0.57
data_arable[P4CType==1|
              P4CType==2]$GapP4CType<-0.6
data_arable[P4CType==13]$GapP4CType<-0.77


#### Add climate data ####

data_temp <- read_climate(data_ref_clean, "Temp")
data_temp$Year<-as.numeric(data_temp$Year)
data_arable$Year<-as.numeric(data_arable$Year)
data_arable<- left_join(data_arable, data_temp,by = c("cell_id", "Year", "OBJECTID")) %>% 
  drop_na() %>% 
  as.data.table()

data_prec <- read_climate(data_ref_clean, "Prec")
data_prec$Year<-as.numeric(data_prec$Year)
data_arable$Year<-as.numeric(data_arable$Year)
data_arable <- left_join(data_arable, data_prec,by = c("cell_id", "Year", "OBJECTID")) %>% 
  drop_na() %>% 
  as.data.table()

rm(data_temp, data_prec, temp, prec)
gc()


#### Add price data ####

data_price <- read_price()
# Original data is published by FAO (https://www.fao.org/) under licence CC BY 4.0 and was downloaded from 
# https://www.fao.org/faostat/en/#data/PP

data_price<-data_price[-2,]

av_price<-list()
my_price<-vector()

for (i in seq_along(data_price$ID)) { 
  for(y in 8:21){
    my_price<-c(my_price, sum(data_price[i,(y-5):(y-1)])/length(data_price[i,(y-5):(y-1)]))
  }
  av_price[[i]]<-my_price
  my_price<-vector()
}

rm(my_price, i, y)
gc()

av_price<-do.call(rbind, av_price)
av_price<-cbind(data_price[,1:2], av_price)
colnames(av_price)<-c(colnames(data_price)[1:2], colnames(data_price)[8:21])

delta_price<-as.matrix(data_price[,-c(1:7)])-as.matrix(av_price[,-c(1:2)])
delta_price<-cbind(data_price[,1:2], delta_price)

av_price<-as.data.table(t(av_price[,3:16]))
colnames(av_price)<-paste("avPrice_CType", data_price$ID, sep="")
av_price$Year<-seq(2010, 2023,1)

delta_price<-as.data.table(t(delta_price[,3:16]))
colnames(delta_price)<-paste("deltaPrice_CType", data_price$ID, sep="")
delta_price$Year<-seq(2010, 2023,1)

data_arable$Year<-as.numeric(data_arable$Year)
av_price$Year<-as.numeric(av_price$Year)
delta_price$Year<-as.numeric(delta_price$Year)

data_arable<-left_join(data_arable, av_price, by="Year")
data_arable<-left_join(data_arable, delta_price, by="Year")

rm(av_price, data_price, delta_price)
gc()


#### Add politics ####

politics<-data_arable %>%
  select("OBJECTID", "Year", "state")%>%
  as.data.table

politics<-politics %>% 
  group_by(Year, state) %>%
  mutate(REL="NO",
         BST=ifelse(Year>2013, "NEG", "NO"),
         PCP=ifelse(Year>2019, "POS", "NO"),
         Greening=ifelse(Year>2014, "POS", "NO"),
         EQS=ifelse(Year>2016, "NEG", "NO"),
         BNN=ifelse(Year>2019, "NEG", "NO"),
         DCP="NO") %>%
  as.data.table()

politics<-politics %>%
  select(-c("OBJECTID", "Year", "state")) %>%
  as.data.table()

data_arable<-data_arable %>%
  bind_cols(politics) %>%
  as.data.table()

rm(politics)
gc()


#### Add soil data ####

data_arable <- inner_join(data_arable, mySoil,by = c("OBJECTID")) %>% 
  as.data.table()

rm(mySoil)
gc()


#### Add rotational typology ####

data_arable$struc <- "test"
data_arable$func <- "test"

data_arable[data_arable$t == 0 & data_arable$k == 1, "struc"] <- "A"
data_arable[data_arable$t %in% 1:2, "struc"] <- "B"
data_arable[data_arable$t %in% 3:4 & data_arable$k == 2, "struc"] <- "C"
data_arable[data_arable$t > 4 & data_arable$k == 2, "struc"] <- "D"
data_arable[data_arable$t %in% 3:4 & data_arable$k == 3, "struc"] <- "E"
data_arable[data_arable$t > 4 & data_arable$k == 3, "struc"] <- "eF"
data_arable[data_arable$t %in% 3:4 & data_arable$k == 4, "struc"] <- "G"
data_arable[data_arable$t > 4 & data_arable$k == 4, "struc"] <- "H"
data_arable[data_arable$t > 2 & data_arable$k > 4, "struc"] <- "I"

data_arable[data_arable$b == 0 & data_arable$s == 0, "func"] <- 1
data_arable[data_arable$b == 0 & data_arable$s > 0 & data_arable$s <= 0.5, "func"] <- 2
data_arable[data_arable$b == 0 & data_arable$s > 0.5, "func"] <- 3
data_arable[data_arable$b > 0 & data_arable$b <= 0.5 & data_arable$s == 0, "func"] <- 4
data_arable[data_arable$b > 0 & data_arable$b <= 0.5 & data_arable$s > 0 & data_arable$s <= 0.5, "func"] <- 5
data_arable[data_arable$b > 0 & data_arable$b <= 0.5 & data_arable$s > 0.5, "func"] <- 6
data_arable[data_arable$b > 0.5 & data_arable$s == 0, "func"] <- 7
data_arable[data_arable$b > 0.5 & data_arable$s > 0 & data_arable$s <= 0.5, "func"] <- 8
data_arable[data_arable$b > 0.5 & data_arable$s > 0.5, "func"] <- 9

data_arable<-data_arable %>% 
  select(!c(t,k,s,b))


#### Change numbered crop types to character ####

crops<-c("MG", "MS", "WW", "WB", "WR", "WTS", "SWTR", "SB", "SO", "WO", "SU", "PO", "LEG", "VEG", "SUN")
data_arable$CType<-as.factor(data_arable$CType)
levels(data_arable$CType)<-crops

data_arable$P1CType<-as.factor(data_arable$P1CType)
levels(data_arable$P1CType)<-crops
data_arable$P2CType<-as.factor(data_arable$P2CType)
levels(data_arable$P2CType)<-crops
data_arable$P3CType<-as.factor(data_arable$P3CType)
levels(data_arable$P3CType)<-crops
data_arable$P4CType<-as.factor(data_arable$P4CType)
levels(data_arable$P4CType)<-crops

rm(crops)
gc()


#### Add neighboring information ####

neighbour <- data_arable %>%
  group_by(cell_id, Year) %>% 
  summarize(
    NP1CType = as.character(names(which.max(table(P1CType)))),
    NP2CType = as.character(names(which.max(table(P2CType)))),
    NP3CType = as.character(names(which.max(table(P3CType)))),
    NP4CType = as.character(names(which.max(table(P4CType))))
  )

data_arable <- left_join(data_arable, neighbour, by = c("cell_id", "Year")) %>% 
  mutate(across(is.numeric, round, digits=2)) %>%
  as.data.table()

data_arable<-data_arable %>%
  drop_na()


#----------------------------------------------------Intermediate storage and clean-up----------------------------------

write_rds(data_arable, "./data/data_arable.RDS")

rm(list=setdiff(ls(), "data_arable"))
gc()


#----------------------------------------------------Model training-----------------------------------------------------

#### Generate environmental clusters ####

# Perform PCA #
enviroPCA<-prcomp(data_arable[,c("TempAM", "TempJJ", "TempAVG",
                                 "PrecAM", "PrecJJ", "PrecAVG",
                                 "SElev", "SSlope",
                                 "Corg", "sand", "clay", "silt")],
                  center = TRUE, scale. = TRUE, retx = TRUE)
summary(enviroPCA)

# Clustering of PCs #
# The 25 clusters are indented to cover the entire study domain including Austria, Germany, and the Czech Republic. #
myCluster<-kmeans(enviroPCA$x[,1:5], centers = 25, nstart = 25)
data_arable$cluster<-myCluster$cluster

rm(enviroPCA, myCluster)
gc()


#### Build model_a: using random train-test split and all available features ####

# Data preparation #

# Disregard variables used for environmental cluster #
myVariables<-colnames(data_arable)[!(colnames(data_arable) %in% c("OBJECTID", "x_coord", "y_coord", "cell_id",
                                                                  "TempAM", "TempJJ", "TempAVG",
                                                                  "PrecAM", "PrecJJ", "PrecAVG",
                                                                  "SElev", "SSlope",
                                                                  "Corg", "sand", "clay", "silt"))]

# Select a random 10%-sample #
set.seed(123)
x<-sample(1:nrow(data_arable), nrow(data_arable)*0.1)
y<-which(colnames(data_arable) %in% myVariables)

mySample<-data_arable[x, ..y]
mySample$CType<-as.character(mySample$CType)
mySample$CType<-as.factor(mySample$CType)

rm(x, y)
gc()

# Split sample into 80% for training and 20% for testing #
set.seed(123)
x<-sample(1:nrow(mySample), nrow(mySample)*0.8)

train_data<-mySample[x,]
test_data<-mySample[-x,]

rm(x)
gc()

# Train model #
model_a<- ranger(CType ~ ., data = train_data,
                 importance="impurity",
                 mtry = floor(ncol(train_data)/3),
                 num.trees=60, 
                 oob.error = TRUE,
                 probability = FALSE,
                 classification = TRUE)

# Testing accuracy #
pred_test<-predict(model_a, data=test_data)
(cm_test<-confusionMatrix(pred_test$predictions, test_data$CType)) 

# Training accuracy #
pred_train<-predict(model_a, data=train_data)
(cm_train<-confusionMatrix(pred_train$predictions ,train_data$CType))

# Variable importance #
# Variable groups for display #
var_group<-c(rep("General", 2),
             rep("IACS", 5),
             rep("Agronomy", 24),
             rep("Prices", 32),
             rep("Politics", 7),
             rep("Agronomy", 2),
             rep("IACS", 4),
             "Environment") %>% 
  as.data.frame()
var_group$Var<-myVariables
colnames(var_group)<-c("Group", "Var")
var_group<-var_group[-3,] %>% 
  as.data.frame()

var_importance_a <- as.data.frame(model_a$variable.importance) 
var_importance_a$ct <- rownames(var_importance_a)
colnames(var_importance_a) <- c("VI", "Var")
var_importance_a<-left_join(var_importance_a, var_group, by="Var") %>% 
  mutate(VI=VI/max(VI)*100) %>% 
  arrange(desc(VI)) %>% 
  filter(VI>0)
var_importance_a$Group<-as.factor(var_importance_a$Group)

var_imp_a <- ggplot(var_importance_a) + 
  geom_col(aes(y = VI, x = Var, fill = Group)) +
  labs(x = "", y = "Relative variable importance", title = "Model a") +
  theme(
    plot.title = element_text(margin = margin(10, 0, 10, 0), size = 14),
    legend.title = element_blank()) +
  theme_bw() +
  scale_fill_manual(values = npg_colors_full[1:6]) +
  scale_x_discrete(limits = var_importance_a$Var, guide = guide_axis(angle = 90))
var_imp_a

rm(cm_test, cm_train, pred_test, pred_train, test_data, train_data, var_importance_a)
gc()



#### Build model_b: using forward feature selection and spatio-temporal cross validation ####

# Pre-model built on 10% sample #

# Create space-time folds #
indices <- CreateSpacetimeFolds(mySample, spacevar="cluster", timevar = "Year", k=5)
# For the entire study domain, we use k=10.

# Train model #
pre_model_b<-ffs(mySample[,-3],mySample$CType, metric="Accuracy",
                 method="rf", tuneGrid=data.frame("mtry"=2),
                 verbose=FALSE,ntree=25,
                 trControl=trainControl(method="cv",
                                        index = indices$index,
                                        savePredictions = "final",
                                        classProbs = T))

global_validation(pre_model_b)


# Variable importance #
var_importance_pre_b <- as.data.frame(varImp(pre_model_b)[[1]]) 
var_importance_pre_b$ct <- rownames(var_importance_pre_b)
colnames(var_importance_pre_b) <- c("VI", "Var")
var_importance_pre_b<-left_join(var_importance_pre_b, var_group, by="Var") %>% 
  mutate(VI=VI/max(VI)*100) %>% 
  arrange(desc(VI)) %>% 
  filter(VI>0)
var_importance_pre_b$Group<-as.factor(var_importance_pre_b$Group)

var_imp_pre_b <- ggplot(var_importance_pre_b) + 
  geom_col(aes(y = VI, x = Var, fill = Group)) +
  labs(x = "", y = "Relative variable importance", title = "Model b") +
  theme(
    plot.title = element_text(margin = margin(10, 0, 10, 0), size = 14),
    legend.title = element_blank()) +
  theme_bw() +
  scale_fill_manual(values = npg_colors_full[c(1,4,6)]) +
  scale_x_discrete(limits = var_importance_pre_b$Var, guide = guide_axis(angle = 90))
var_imp_pre_b

rm(indices, mySample, myVariables, var_importance_pre_b)
gc()

# Final model #

# Data preparation #

# Create space-time folds #
indices <- CreateSpacetimeFolds(data_arable, spacevar="cluster", timevar = "Year", k=5)
# For the entire study domain, we use k=10.

# Select variables from pre_model_b as predictors #
myVariables<-c(rownames(pre_model_b$finalModel$importance), "CType")

y<-which(colnames(data_arable) %in% myVariables)

mySample<-data_arable[,..y]
mySample$CType<-as.character(mySample$CType)
mySample$CType<-as.factor(mySample$CType)

rm(y)
gc()

model_b<-train(mySample[,-1],mySample$CType, metric="Accuracy",
               method="rf", tuneGrid=data.frame("mtry"=5),
               verbose=FALSE,ntree=60,
               trControl=trainControl(method="cv",
                                      index = indices$index,
                                      savePredictions = "final",
                                      classProbs = T))

global_validation(model_b)


# Variable importance #
var_importance_b <- as.data.frame(varImp(model_b)[[1]]) 
var_importance_b$ct <- rownames(var_importance_b)
colnames(var_importance_b) <- c("VI", "Var")
var_importance_b<-left_join(var_importance_b, var_group, by="Var") %>% 
  mutate(VI=VI/max(VI)*100) %>% 
  arrange(desc(VI)) %>% 
  filter(VI>0)
var_importance_b$Group<-as.factor(var_importance_b$Group)

var_imp_b <- ggplot(var_importance_b) + 
  geom_col(aes(y = VI, x = Var, fill = Group)) +
  labs(x = "", y = "Relative variable importance", title = "Model b") +
  theme(
    plot.title = element_text(margin = margin(10, 0, 10, 0), size = 14),
    legend.title = element_blank()) +
  theme_bw() +
  scale_fill_manual(values = npg_colors_full[c(1,4,6)]) +
  scale_x_discrete(limits = var_importance_b$Var, guide = guide_axis(angle = 90))
var_imp_b


# Confusion matrix #

conf_df <- as.data.frame(as.table(confusionMatrix(model_b)$table))
colnames(conf_df) <- c("Predicted", "Actual", "Count")

conf_df <- conf_df %>%
  group_by(Actual) %>%
  mutate(Percentage = Count / sum(Count) * 100)

confusion_b<-ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Percentage)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Percentage, 1)), color = "black", size = 5) +
  scale_fill_gradient(low = "white", high = test_color) + 
  theme_minimal() +
  labs(title = "Confusion Matrix Heatmap",
       x = "Actual Class",
       y = "Predicted Class",
       fill = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
confusion_b

rm(conf_df, indices, mySample, myVariables, var_importance_b)
gc()


#### Build model_c: disregarding cropping history ####

# Data preparation #

# Disregard variables related to cropping history and environmental cluster #
myVariables<-colnames(data_arable)[!(colnames(data_arable) %in% c("OBJECTID", "x_coord", "y_coord", "cell_id",
                                                                  "P1CType", "P2CType", "P3CType", "P4CType",
                                                                  "WinterP1CType", "WinterP2CType", "WinterP3CType", "WinterP4CType",
                                                                  "LeafP1CType", "LeafP2CType", "LeafP3CType", "LeafP4CType",
                                                                  "DrawP1CType", "DrawP2CType", "DrawP3CType", "DrawP4CType",
                                                                  "OrganP1CType", "OrganP2CType", "OrganP3CType", "OrganP4CType",
                                                                  "DroughtP1CType", "DroughtP2CType", "DroughtP3CType", "DroughtP4CType",
                                                                  "GapP1CType", "GapP2CType", "GapP3CType", "GapP4CType",
                                                                  "NP1CType", "NP2CType", "NP3CType", "NP4CType",
                                                                  "cluster"))]

# Select a random 10%-sample #
set.seed(123)
x<-sample(1:nrow(data_arable), nrow(data_arable)*0.1)
y<-which(colnames(data_arable) %in% myVariables)
mySample<-data_arable[x, ..y]

mySample$CType<-as.character(mySample$CType)
mySample$CType<-as.factor(mySample$CType)

rm(x, y)
gc()

# Create space-time folds #
indices <- CreateSpacetimeFolds(mySample, spacevar="cluster", timevar = "Year", k=5)
# Same as above:for the entire study domain, we use k=10.

# Train model #
model_c<- ffs(mySample[,-3],mySample$CType, metric="Accuracy",
              method="rf", tuneGrid=data.frame("mtry"=2),
              verbose=FALSE,ntree=25,
              trControl=trainControl(method="cv",
                                     index = indices$index,
                                     savePredictions = "final",
                                     classProbs = T))

global_validation(model_c)


# Variable importance #
var_importance_c <- as.data.frame(varImp(model_c)[[1]]) 
var_importance_c$ct <- rownames(var_importance_c)
colnames(var_importance_c) <- c("VI", "Var")
var_importance_c<-left_join(var_importance_c, var_group, by="Var") %>% 
  mutate(VI=VI/max(VI)*100) %>% 
  arrange(desc(VI)) %>% 
  filter(VI>0)
var_importance_c[1:2,3]<-"Environment"
var_importance_c$Group<-as.factor(var_importance_c$Group)


var_imp_c <- ggplot(var_importance_c) + 
  geom_col(aes(y = VI, x = Var, fill = Group)) +
  labs(x = "", y = "Relative variable importance", title = "Model c") +
  theme(
    plot.title = element_text(margin = margin(10, 0, 10, 0), size = 14),
    legend.title = element_blank()) +
  theme_bw() +
  scale_fill_manual(values = npg_colors_full[c(1:3,5)]) +
  scale_x_discrete(limits = var_importance_c$Var, guide = guide_axis(angle = 90))
var_imp_c

rm(indices, mySample, myVariables, var_importance_c)
gc()


#----------------------------------------------------Intermediate storage and clean-up----------------------------------

write_rds(model_a, "./data/model_a.RDS")
write_rds(pre_model_b, "./data/pre_model_b.RDS")
write_rds(model_b, "./data/model_b.RDS")
write_rds(model_c, "./data/model_c.RDS")

rm(list=setdiff(ls(), c("data_arable", "model_b", "crop_color")))
gc()


#----------------------------------------------------Rotation projections---------------------------------------------

AT<-data_arable%>%
  group_by(OBJECTID) %>%
  filter(Year == max(Year)) %>%
  slice(1) %>%
  ungroup() 
AT<-AT %>%
  filter(Year>2022)
AT$Year<-as.numeric(AT$Year)

#### Get price information ####
  
read_price <- function() { 
  data_price_raw <- fread("./data/prices_Austria.csv", sep=",", header = T)
  
  data_p <- data_price_raw %>% 
    mutate(Item = recode(Item, "Vegetables, leguminous nes" = "Vegetables, leguminous")) %>% 
    filter(Unit == 'USD') %>%
    as.data.table()
  
  data_p <- data_p %>% 
    pivot_wider(id_cols = Item, values_from = Value,names_from = Year) %>% 
    as.data.table()
  
  k <- which(is.na(data_p), arr.ind=TRUE)
  data_p[k] <- rowMeans(data_p[,-1], na.rm=TRUE)[k[,1]]
  
  
  for (n in c(2:20)) {
    data_p[Item=="Lupins"|Item=="Soya beans"][[n]]<-mean(data_p[Item=="Lupins"|Item=="Soya beans"][[n]])
    data_p[Item=="Cabbages"|Item=="Cauliflowers and broccoli"][[n]]<-mean(data_p[Item=="Cabbages"|Item=="Cauliflowers and broccoli"][[n]])
  }
  
  data_p<-data_p[-4,]
  
  data_p$Item[2]<-"Brassica"
  data_p$Item[3]<-"Carrots"
  data_p$Item[4]<-"Maize" 
  data_p$Item[6]<-"Onions"
  data_p$Item[8]<-"Canola"
  data_p$Item[10]<-"Legumes"
  data_p$Item[12]<-"Sunflower"
  
  scheme <- fread("./data/crop_classification.csv", sep=";", header = TRUE)
  scheme$`Price Description` <- tolower(scheme$`Price Description`)
  data_p$Item <- tolower(data_p$Item)
  
  temp <- scheme %>% 
    select(`Crop Class ID`, `Price Description`) %>% 
    as.data.table()
  colnames(temp) <- c('ID', 'Item')
  temp[7,2] <- "triticale"
  
  data_price <- inner_join(temp, data_p, by = "Item") %>% 
    as.data.table()
}

myPrice<-read_csv("./data/globiom_Austria_RCP8p5.csv") %>%
  select(YEAR, delta_Price_CType13)%>%
  filter(YEAR>2023)
colnames(myPrice)<-c("Year", "deltaPrice_CType13")

data_price <- read_price()
data_price<-data_price[-2,]

av_price<-list()
my_price<-vector()

for (i in seq_along(data_price$ID)) {
  for(y in 8:21){
    my_price<-c(my_price, sum(data_price[i,(y-5):(y-1)])/length(data_price[i,(y-5):(y-1)]))
  }
  av_price[[i]]<-my_price
  my_price<-vector()
}

rm(my_price, i, y)
gc()

av_price<-do.call(rbind, av_price)
av_price<-cbind(data_price[,1:2], av_price)
colnames(av_price)<-c(colnames(data_price)[1:2], colnames(data_price)[8:21])

delta_price<-as.matrix(data_price[,-c(1:7)])-as.matrix(av_price[,-c(1:2)])
delta_price<-cbind(data_price[,1:2], delta_price)

av_price<-as.data.table(t(av_price[,3:16]))
colnames(av_price)<-paste("avPrice_CType", data_price$ID, sep="")
av_price$Year<-seq(2010, 2023,1)

delta_price<-as.data.table(t(delta_price[,3:16]))
colnames(delta_price)<-paste("deltaPrice_CType", data_price$ID, sep="")
delta_price$Year<-seq(2010, 2023,1)

delta_price<-delta_price %>%
  filter(Year<2024) %>%
  select(Year, deltaPrice_CType13)

myPrice<-rbind(delta_price, myPrice)

rm(av_price, data_price, delta_price, read_price)
gc()
  
  
#### Create future rotation predictions ####
  
myYear<-sort(unique(AT$Year))
myPrediction<-list()
for (m in 1:length(myYear)) {
  myAT<-AT%>%
    filter(Year==myYear[m]) %>%
    select(!deltaPrice_CType13)
  newYear<-seq(min(myAT$Year)+1,2069,1)
  
  AT_predictions<-list()
  for (y in 1:length(newYear)) {
    newCrop<-myAT %>%
      mutate(OBJECTID=myAT$OBJECTID,
             Year=newYear[y],
             CType=NA,
             P3CType=myAT$P2CType,
             P4CType=myAT$P3CType,
             P2CType=myAT$P1CType,
             P1CType=myAT$CType,
             func=myAT$func,
             GapP3CType=0.75,
             DrawP1CType="MED",
             GapP1CType=0.75,
             DrawP4CType="MED",
             WaterP4CType="LOW",
             OrganP4CType="GRAIN",
             GapP4CType=0.75,
             WinterP4CType="NO",
             LeafP4CType="YES") %>%
      as.data.table()
    
    levels(newCrop$P1CType)<-c(1:17)
    levels(newCrop$P2CType)<-c(1:17)
    levels(newCrop$P3CType)<-c(1:17)
    levels(newCrop$P4CType)<-c(1:17)
    
    newCrop[P3CType==4|
              P3CType==5|
              P3CType==8]$GapP3CType<-0.33
    newCrop[P3CType==3|
              P3CType==6|
              P3CType==7]$GapP3CType<-0.5
    newCrop[P3CType==16|
              P3CType==17]$GapP3CType<-0.57
    newCrop[P3CType==1|
              P3CType==2]$GapP3CType<-0.6
    newCrop[P3CType==13]$GapP3CType<-0.77
    
    newCrop[P1CType==13|
              P1CType==16]$DrawP1CType<-c("LOW")
    newCrop[P1CType==1|
              P1CType==2|
              P1CType==10|
              P1CType==11|
              P1CType==12|
              P1CType==14|
              P1CType==15]$DrawP1CType<-c("HIGH")
    
    newCrop[P1CType==4|
              P1CType==5|
              P1CType==8]$GapP1CType<-0.33
    newCrop[P1CType==3|
              P1CType==6|
              P1CType==7]$GapP1CType<-0.5
    newCrop[P1CType==16|
              P1CType==17]$GapP1CType<-0.57
    newCrop[P1CType==1|
              P1CType==2]$GapP1CType<-0.6
    newCrop[P1CType==13]$GapP1CType<-0.77
    
    newCrop[P4CType==13|
              P4CType==16]$DrawP4CType<-c("LOW")
    newCrop[P4CType==1|
              P4CType==2|
              P4CType==10|
              P4CType==11|
              P4CType==12|
              P4CType==14|
              P4CType==15]$DrawP4CType<-c("HIGH")
    
    newCrop[P4CType==3|
              P4CType==6|
              P4CType==10|
              P4CType==11|
              P4CType==13]$WaterP4CType<-c("MED")
    newCrop[P4CType==4|
              P4CType==5]$WaterP4CType<-c("HIGH")
    
    newCrop[P4CType==2|
              P4CType==11|
              P4CType==12|
              P4CType==16|
              P4CType==17]$OrganP4CType<-c("ROOT")
    
    newCrop[P4CType==4|
              P4CType==5|
              P4CType==8]$GapP4CType<-0.33
    newCrop[P4CType==3|
              P4CType==6|
              P4CType==7]$GapP4CType<-0.5
    newCrop[P4CType==16|
              P4CType==17]$GapP4CType<-0.57
    newCrop[P4CType==1|
              P4CType==2]$GapP4CType<-0.6
    newCrop[P4CType==13]$GapP4CType<-0.77
    
    newCrop[P4CType==3|
              P4CType==4|
              P4CType==5|
              P4CType==6|
              P4CType==10]$WinterP4CType<-c("YES")
    
    newCrop[P4CType==1|
              P4CType==2|
              P4CType==3|
              P4CType==4|
              P4CType==5|
              P4CType==6|
              P4CType==7|
              P4CType==8|
              P4CType==9]$LeafP4CType<-c("NO")
    
    newCrop<-left_join(newCrop, myPrice, by="Year")
    
    crops<-c("MG", "MS", "WW", "WB", "WR", "WTS", "SWTR", "SB", "SO", "WO", "SU", "PO", "LEG", "VEG", "SUN", "ON", "CA")
    newCrop$P1CType<-as.factor(newCrop$P1CType)
    levels(newCrop$P1CType)<-crops
    newCrop$P2CType<-as.factor(newCrop$P2CType)
    levels(newCrop$P2CType)<-crops
    newCrop$P3CType<-as.factor(newCrop$P3CType)
    levels(newCrop$P3CType)<-crops
    newCrop$P4CType<-as.factor(newCrop$P4CType)
    levels(newCrop$P4CType)<-crops
    
    newCrop$state<-as.factor(newCrop$state)
    newCrop$Year<-as.numeric(newCrop$Year)
    newCrop$CType<-as.factor(newCrop$CType)
    newCrop$P3CType<-as.factor(newCrop$P3CType)
    newCrop$P4CType<-as.factor(newCrop$P4CType)
    newCrop$P2CType<-as.factor(newCrop$P2CType)
    newCrop$P1CType<-as.factor(newCrop$P1CType)
    newCrop$func<-as.factor(newCrop$func)
    newCrop$GapP3CType<-as.numeric(newCrop$GapP3CType)
    newCrop$DrawP1CType<-as.factor(newCrop$DrawP1CType)
    newCrop$GapP1CType<-as.numeric(newCrop$GapP1CType)
    newCrop$DrawP4CType<-as.factor(newCrop$DrawP4CType)
    newCrop$WaterP4CType<-as.factor(newCrop$WaterP4CType)
    newCrop$OrganP4CType<-as.factor(newCrop$OrganP4CType)
    newCrop$GapP4CType<-as.numeric(newCrop$GapP4CType)
    newCrop$WinterP4CType<-as.factor(newCrop$WinterP4CType)
    newCrop$deltaPrice_CType13<-as.numeric(newCrop$deltaPrice_CType13)
    newCrop$LeafP4CType<-as.factor(newCrop$LeafP4CType)
    
    tryCatch({
      CType <- predict(model_b, newCrop, type = "prob")
      myCType <- apply(CType, 1, function(row) names(CType)[which.max(row)])
      
      x <- which(apply(CType, 1, max) < summary(apply(CType, 1, max))[[2]])
      
      if (length(x) > 0) {
        myCType[x] <- sample(names(CType),length(x), replace = T)
      }
      
      newCrop$CType <- as.factor(myCType)
      
    }, error = function(e) {
      message("Error occurred, using fallback prediction.")
      
      newCrop$CType <- predict(model_b, newCrop)
    })
    
    myAT<-newCrop %>%
      select(!deltaPrice_CType13)
    AT_prediction[[y]]<-newCrop %>%
      select(OBJECTID, Year, CType)
    
    rm(CType, myCType, newCrop, x)
  }
  
  AT_prediction<-do.call(rbind, AT_prediction)
  
  AT_prediction<-AT_predictions %>%
    pivot_wider(names_from=Year, values_from=CType)
  
  myPrediction[[m]]<-AT_predictions
}


AT_prediction<-data_arable %>%
  select(OBJECTID, Year, CType)


for (m in 1:length(myPrediction)) {
  myData<-myPrediction[[m]] %>%
    pivot_longer(!OBJECTID,names_to="Year", values_to="CType")
  
  AT_prediction<-rbind(AT_prediction, myData) %>%
    arrange(OBJECTID, Year)
  
  rm(myData)
}

x<-which(AT_prediction$CType=="WW" & AT_prediction$Year>2023)
x <- sample(x, length(x)*0.23)
AT_prediction$CType[x]<-sample(unique(AT_prediction$CType)[unique(AT_prediction$CType)!="WW"],length(x), replace = T)

#### Plot results ####

statistics<-AT_prediction %>%
  group_by(Year, CType) %>%
  summarise(crops=n()) %>%
  as.data.table()

statistics_new<-statistics%>%
  group_by(Year)%>%
  summarise(fields=sum(crops))%>%
  as.data.table()

statistics<-left_join(statistics, statistics_new, by="Year")
statistics$rel_crops<-statistics$crops/statistics$fields*100
statistics$CType<-as.factor(statistics$CType)

stat_graph<-ggplot(statistics, aes(x=Year, y=rel_crops, group=CType))+
  geom_line(aes(color=CType), linewidth=1.5)+
  scale_color_manual(values = crop_color)+
  theme_bw()+
  labs(x="Year", y="Relative crop type distribution", title="(b) Austria")+
  scale_x_discrete(breaks = seq(2013, 2063, 10),labels=seq(2013, 2063, 10))+
  theme(legend.title=element_blank())
stat_graph

  
#----------------------------------------------------Intermediate storage and clean-up----------------------------------
  
  write_rds(AT_prediction, "./data/AT_prediction.RDS")
  
  rm(list=setdiff(ls(), c("AT_prediction","crop_color")))
  gc()


