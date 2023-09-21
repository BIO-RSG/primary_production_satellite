# Primary production

*00_background* subfolder contains notes / images / research papers in pdf format for reference.  

This repo contains the scripts to calculate primary production (PP) using the BIOv2 model, described [here](markdown_reports/pp_model_description.html). There are also scripts to calculate PP using the BIOv1 and Laval models, but currently (2021-12-29) they require some cleaning and reformatting since they haven't been used in awhile and changes have been made to input/output formats and the "run" script.  

Note that some model functions have been converted to C++ for the sake of speed, and are run by using the `Rcpp` package. (No changes needed in the code, just adjust variables at the top of __02_run_PP.R__ as usual).  


--------------------------------------------------------------------------------

## Data sources

### CHL, PAR

**These are stored in the PanCanadian dataset (located on _sinumkw_ as of 2023-08-23). If you don't have access to this, you can download them using the following steps:**  

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

1. Gather raw input data from your selected sources, for whatever spatial and temporal resolution    

2. Combine input into the required format, an .Rdata file containing the following columns:  

    - **event_id**
    - **sample_id**
    - **station**
    - **latitude**
    - **longitude**
    - **doy** (day of year)
    - **year**
    - **par** (surface PAR, currently we use the satellite value)
    - **depth_data** (a list-col where each cell/element is a dataframe containing columns depth, alphab, pmb, chl)  
    - **time_data** (a list-col where each cell/element is a vector containing one column ("time"), the time(s) of day to calculate pp, **in UTC** in the format **HHMM**)

* For BIOv1, PI parameters and biomass profile parameters for satellite PP were retrieved using the lookup tables in *01a_LUT_cp-20141117.csv* and *01a_LUT_pi-20141113.csv*, based on satellite SST and Chl-a. The function to gather the data can be found in `scripts/get_pi_bp_from_LUT.R`.  

3. Calculate PP using the "01_run" scripts - either satellite PP (surface only) or in situ PP (a matrix of PP for the selected times and depths for each record)  

### OUTPUT

1. Rdata file: This contains the full results with PP for each time/depth for each record  

2. csv file: This contains a simplified table with a single PP value per record, integrated over time and depth  



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
