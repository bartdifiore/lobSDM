library(tidyverse)
library(stringr)
library(sf)
library(raster)
library(lubridate)
library(gmRi)
# devtools::install_github("rstudio/reticulate")
library(reticulate)
library(ncdf4)

source("Code/enhance_r_funcs.R")

#####
## Preliminary set up
#####
glorys_path<- cs_path(box_group = "RES_Data", subfolder = "GLORYs/NE_Shelf_TempSal/")

# Some quick exploration to get familiar with things
list.files(glorys_path)

t_nc <- nc_open(list.files(glorys_path, full.names = TRUE)[[1]])
print(t_nc)
t_nc$dim

#####
## Sourcing Adam's python script for extracting bottom temperature
#####
# Adam has a nice script for extracting bottom temperature here (https://github.com/adamkemberling/glorys_northeast/blob/main/py/bottom_layer_extraction.py). It was a bit nutty to try to get this all working, I am wayyy out of practice with this stuff and Adam is leap years ahead of me. For what it's worth, here's what I had to do (after making sure python, jupetyr and such were all set)

# Open terminal and type `python -m venv myenv`
# pip install ipykernel
# python -m ipykernel install --user --name=aallyn --display-name "Python (myenv)" 

# Next, I opened the .ipynb and in the top right I had to basically get the python thing to be working with the virtual environment I just set up.

# After doing that, I was then able to run .ipynb cells and debug as needed (for example, I had to run `pip install xarray` and do the same call for any other...modules? missing?)
glorys_new_path<- cs_path(box_group = "RES_Data", subfolder = "GLORYs/NE_Shelf_MonthlyTemps/")
t<- raster(paste0(glorys_new_path, "BT.nc"))
t
plot(t[[1]])

t2<- raster(paste0(glorys_new_path, "SST.nc"))
t2
plot(t2[[1]])

# Looks good! 
