# Inference for extreme spatial temperature events in a changing climate with application to Ireland
Here is the code for the paper "Inference for extreme spatial temperature events in a changing climate with application to Ireland", which you can read on arXiv. Feel free to use/adapt this code :) I cannot publish the data analysed for this project in this repository but details on how to download it are [here](#downloading-data). To get an overview of the project structure you can jump down here to the [file tree](#file-tree).


## Downloading data
--------------------- 
### Observational data
- Data for sites in the Republic of Ireland can be downloaded from [Met Éireann](https://www.met.ie/climate/available-data/historical-data).
- Data from the sites in Northern Irealnd from [Met Office Integrated Data Archive System (MIDAS)](https://data.ceda.ac.uk/badc/ukmo-midas-open/data/uk-daily-temperature-obs).

### Climate model data
- Climate model data used can be accessed through the [Copernicus Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/projections-cordex-domains-single-levels). Selecting `Europe` > `Historical` > `Daily mean` > ` Maximum 2m temperature in the last 24 hours` > `ICHEC-EC-EARTH (Ireland)` > `CLMcom-CLM-CCLM4-8-17 (EU)` will give the climate model data used in our analysis.


### Covariates
- The global mean temperature used in our analyis can be downloaded from [HadCRUT](https://crudata.uea.ac.uk/cru/data/temperature/).
- Data from Irish C0$_2$ emissions can be downloaded from [Our World in Data](https://github.com/owid/co2-data).

## Using the code
--------------------- 
