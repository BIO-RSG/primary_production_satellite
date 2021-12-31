# Primary production

*00_background* subfolder contains notes / images / research papers in pdf format for reference.  

This repo contains the scripts to calculate primary production (PP) using the BIOv2 model, described [here](markdown_reports/pp_model_description.html). There are also scripts to calculate PP using the BIOv1 and Laval models, but currently (2021-12-29) they require some cleaning and reformatting since they haven't been used in awhile and changes have been made to input/output formats and the "run" script.  

Note that some model functions have been converted to C++ for the sake of speed, and are run by using the `Rcpp` package. (No changes needed in the code, just adjust variables at the top of __02_run_PP.R__ as usual).  


--------------------------------------------------------------------------------

## Data sources

### CHL, PAR

**These are stored in the PanCanadian dataset (located on _hecla_ as of 2020-11-20). If you don't have access to this, you can download them using the following steps:**  

- Download daily level-3 (L3) binned VIIRS-SNPP files from the [NASA ocean colour website](https://oceandata.sci.gsfc.nasa.gov/) (*note: you must download daily instead of 8day files for consistency with the formatting script*):  
This can be done using the *01a_download_daily_L3b.sh* from the [panCan_processing repo](https://github.com/BIO-RSG/panCan_processing). **Note** that you need an Earthdata account to download the files (see top of *01a_download_daily_L3b.sh* script for instructions).  

- Subset the files to the PanCanadian grid (panCan_processing: *02a_subset_L3b_to_panCan_grid.py*)  

- Transfer the files to the `00c_other_variables` subfolder in the primary production repo, sorted into subfolders following the convention `sensor/variable/PANCAN/year`, where sensor is MODIS or VIIRS-SNPP, variable is SST, PAR, or CHL_OCX, and year is a 4-digit number.  

- When running `01a_format_PP_input.R`, make sure to adjust the input_path variable to the `00c_other_variables` folder.  

**NOTE:** In each script, don't forget to adjust input/output paths to the appropriate directories.  


### SST, TCC (for BIOv1 and Laval only)

**SST**:  
Follow instructions for CHL and PAR above.  
*Note that SST is also required input for the lookup tables for PI and biomass parameters from the files listed below, using the functions in `01b_get_param.R`.  

**TCC**:  
[MODIS-Aqua 8day cloud fraction](https://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_E_CLD_FR)  
How to download:  

- Edit variables in __00a_TCC_download.R__ and run it to create a list of TCC files to download, written to txt files in the `00a_downloaded_TCC/interval` subfolder, where *interval* is either *8day* or *monthly*, and the txt files are in the naming format __*year*_get_tcc.txt__.  
- Navigate to the folder containing the txt files, and type on the command line: `wget --random-wait --content-disposition -i [filename here, including .txt extension]`  
- Convert TCC .png files to .csv by editing the variables in the __00b_TCC_to_csv.R__ script and running it.  


### PI curve parameters (alphaB, PmB)

`02a_LUT_pi-20141113.csv`  

### Biomass profile parameters (h, sigma, B0, zm)

`02a_LUT_cp-20141117.csv`  
*Note: this is only used if you have set uniform_profile = FALSE in 02_run_PP.R*  

<br>


--------------------------------------------------------------------------------

## INSTRUCTIONS

### STEP 00: Gather raw input data

Using the instructions in the data source section above, gather the necessary input.  

### STEP 01: Format input

Collect and format input data by editing the variables in __01_format_PP_input.R__ and running it. This merges the input data for each variable, with their bin numbers, latitudes, longitudes, and bathymetry. CHL, SST, and time of year are used as input to the function in __01b_get_param.R__, which retrieves corresponding BP (biomass profile) and PI (photosynthesis-irradiance) curve parameters, based on the look-up tables in *01a_LUT_cp-20141117.csv* and *01a_LUT_pi-20141113.csv*. All the data is then combined into a format that can be used as input to the primary production scripts, pixels with missing input data are removed, and the table is written to csv.  


### STEP 02: Calculate PP

Edit the variables in __02_run_PP.R__ and run them to create primary production files using the selected model, which will be written to the `03_PP_output/model/interval` subfolder, where *model* is either *BIOv1*, *BIOv2*, or *Laval*, and *interval* is either *daily*, *8day* or *monthly*. This "run" script take the user input and reads and organizes the formatted csv input, then calls the main function for the selected model. Two output files are created:  

- **NetCDF (.nc)**: This contains only one variable, the primary production in 4km-resolution, in the same size as the NWA (Northwest Atlantic) binned grid so that the values can be easily matched to NWA bin numbers.  

- **Rdata (.rda)**: This contains the input data, primary production, and other variables outputted by the model function (for example, surface PAR).  


### STEP 03: Examine output

You can compare the output to the other models or, for example, the [VGPM](http://orca.science.oregonstate.edu/2160.by.4320.8day.hdf.vgpm.v.chl.v.sst.php) model. Output netCDFs can be read like so:  

```
PP_filename <- "YOUR_PP_FILE.nc"
library(ncdf4)
library(raster)
library(oceancolouR)
nc <- nc_open(PP_filename)
PP <- ncvar_get(nc, "PP")
nc_close(nc)
data("nwa_bins_4km")
PP_rast <- var_to_rast(data.frame(bin=nwa_bins_4km, PP=PP, stringsAsFactors=FALSE), ext=c(-95,-42,39,82))
make_raster_map(PP_rast)
```

Or you can view the output in the `rda` file, which contains the input variables as well:  

```
load(PP_filename)
head(rda_output)

# view a map of PAR input for this primary production
library(tidyr)
par_df <- rda_output %>% dplyr::select(bin, par)
par_rast <- var_to_rast(df=par_df, ext=c(-95,-42,39,82))
make_raster_map(par_rast)
```



--------------------------------------------------------------------------------

## BIOv2 model (currently in use)

[Description of the model](markdown_reports/pp_model_description.html)  

[Comparison to other models](markdown_reports/PP_comparison_summary.html)  


--------------------------------------------------------------------------------

## BIOv1 model

Brief summary of the steps in the PP calculation in *model_BIO.R*:  

1. Get *chla* profile (uniform or non-uniform)  
2. Calculate angles and time variables  
3. Loop over time intervals/zenith angles. At each iteration:  
    a. "bird" function: Calculate direct and diffuse components of surface irradiance for all the wavelengths  
    b. Correct for seasonal variation and cloud cover  
4. Scale/correct irradiance based on satellite PAR  
5. Loop over time intervals/zenith angles. At each iteration:  
    a. Calculate reflection losses at air-sea interface  
    b. Convert values to get the right units for sam_penguin  
    c. "sam_penguin" function: Calculate subsurface irradiance down to euphotic depth (wavelength integration happens in here)  
    d. Integrate over depth  
6. Integrate over time  

Note there are two functions in *model_BIO.R*: *dwcpn()* and *dwcpn_C()*. They are identical except that dwcpn uses *bird()* and *sam_penguin()* for the surface and subsurface irradiance calculations respectively, while dwcpn_C uses *bird_C()* and *sam_penguin_C()*, which are C++ versions of the functions that are much faster.  

--------------------------------------------------------------------------------

## Laval model
