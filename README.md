# Primary production

This repo contains the scripts to calculate primary production (PP) using the BIOv2 model. Note that some model functions have been converted to C++ for the sake of speed, and are run by using the `Rcpp` package. (No changes needed in the code, just adjust variables at the top of __01_run_PP.R__ as usual).  

Instructions and required input formatting can be found at the top of __01_run_PP.R__.

Output:  
- Rdata file containing the input dataframe, as well as columns containing the depth_data filled to the depths where PP is calculated, PP (a matrix where row=PP at depth, col=PP at time), and scalar irradiance E0z (also a matrix with row=E0z at depth, col=E0z at time)   
- csv file containing the input columns + total PP integrated over time and depth, + total PP integrated over time and depth in the top 30m of the water column   

If you're calculating satellite PP and have a lot of pixels in the area of interest, you can use __01_run_PP_satellite_simple.R__ to run a faster version of the script with simplified output: an fst or csv file containing total PP, surface PP, and % PP down to different depths for each pixel.  
