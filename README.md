# Crop roations in Central europe

This repository is the official implementation of [Palka et al. (2025)](https://xyz).

For this paper, we
* combine  __field-specific cropping histories__, __neighbouring rotations__, __environmental conditions__, __agronomic rules__ for good rotation practice, crop __commodity prices__, and agricultural __policies and subsidies__
* train different __random forest (RF) models__ to
* identify __drivers shaping operational crop rotations over Central Europe__ and
* project __project crop rotations until 2070__

accross Central Europe, including Germany, Austria, and the Czech Republic.

## Repository structure

In the __data__ folder, you find the following sub-folders:
* __IACS__: listing Geospatial Application Data (GSA) of the EU’s Integrated Administration and Control System (IACS) for Austria from 2015 to 2023 as shapefiles (Publisher: [Agrarmarkt Austria](https://www.ama.at/))
* __soil__: holding Austrias digital soil map (boka_1km_2016) as shapefile (Publisher: [BFW](https://www.bfw.gv.at/)) and elevation (elev.tif) and slope (slope.tif) information according to Austria's digital elevation model (Publisher: [geoland.at](http://www.geoland.at))
* __weather__: listing daily gridded precipitation (_RR) and temperature (_TM) data over Austria from 2010 until 2023, as one netCDF file per year (Publisher: [GeoSphere Austria](https://www.geosphere.at/))
* a __crop description key__ for classifying GSA string descriptions for crops into numeric crop types (crop_description_Austria.csv)
* a __crop classification key__ describing the different crop types (crop_classification.csv)
* past __crop commodity prices__ (prices_Austria.csv) (Publisher: [FAO](https://www.fao.org/))
* future __crop commodity prices__ according to GLOBIOM projections (globiom_Austria_RCP8p5.csv)

__crop_rotations_Austria.R__ is the R script needed to analyse the data.

### Code sturcture

We provide code for Austria as a sample, where all necessary data is available open-source. As the manuscript is based on the entire study domain, results using this script will differ from the ones reported in the manuscript.<br />

Please be aware: this project runs on a HPC unit at ZALF. Local machines might not allow the processing of large data sets.

The script is structured into the following sections:
* __Prep work__, to set the working directory and load required packages, set colours for display, and build custom functions
* __Create data base__, to combine all of the above data
* __Model training__, to train and assess different RF model
* __Rotation projections__, to use the best RF model to project crop rotations

The following packages are required to run the code:
```
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
```

## Data sources
Open-source __GSA data__ was downloaded from and can be found here: <br />
* __Austria__: https://www.data.gv.at/suche/?katFilter%5B0%5D=httppublicationseuropaeuresourceauthoritydata-themeagri&searchterm&typeFilter%5B0%5D=dataset&nr=1&tagFilter%5B0%5D=INVEKOS

* __Brandenburg__: https://geobroker.geobasis-bb.de/gbss.php?MODE=GetProductInformation&PRODUCTID=996f8fd1-c662-4975-b680-3b611fcb5d1f

* __Czech Republic__: https://mze.gov.cz/public/portal/mze/farmar/LPIS/export-lpis-rocni-shp 

* __Lower-Saxony__: https://sla.niedersachsen.de/landentwicklung/LEA/ 

* __North Rhine-Westphalia__: https://www.opengeodata.nrw.de/produkte/umwelt_klima/bodennutzung/landwirtschaft/


__Digital terrain data__ is available via https://www.adv-online.de/Products/Geotopography/Digital-Terrain-Models/DGM10/ and https://www.data.gv.at/katalog/dataset/dgm.

__Digital soil maps__ can be found at https://github.com/zalf-rpm/Buek200_by_CLC and https://bodenkarte.at.

Gridded __weather data__ is available via https://www.dwd.de/DE/leistungen/cdc/cdc_ueberblick-klimadaten.html and https://data.hub.geosphere.at/dataset/spartacus-v2-1d-1km.


FAO __producer prices__ can be found at https://www.fao.org/faostat/en/#data/PP.

All of the above is available under licence CC BY 4.0 (CC BY 3.0 for older versions).

## Publication

 __Cropping history, agronomic rules, and commodity prices shape crop rotations across Central Europe__

## Authors

__Marlene Palka__ [palkamarlene](https://github.com/palkamarlene); marlene.palka@zalf.de

## License

This project is licensed under the xyz License - see the [LICENSE.md](LICENSE.md) file for details
