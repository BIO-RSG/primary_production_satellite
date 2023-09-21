# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-04-18

# Generating an image for claysa_projects.pptx

library(oceancolouR)
library(dplyr)
library(data.table)
library(corrplot)
source("markdown_reports/PP_comparison_rmd_functions.R")

path <- "model_comparisons/"
date_bio <- "201807"
date_op <- "2018182"
xlims <- c(-95,-42)
ylims <- c(39,82)

file_biov1_LHTBI <- paste0(path,"BIO_v1_Carla/",date_bio,"_LHTBI.csv")
file_biov1_pi <- paste0(path,"BIO_v1_Carla/",date_bio,"_pi-20141113.csv")
file_biov1_pp <- paste0(path,"BIO_v1_Carla/",date_bio,"_pp.csv")
file_svgpm <- paste0(path,"VGPM/vgpm.viirs.",date_op,".hdf")
file_biov2_month_MPI_U <- paste0(path,"BIOv2-monthly_vs_BIOv1-monthly/BIOv2_medianPI_U_BIO_NWA_monthly_",date_bio,".rda")

df_biov1 <- dplyr::bind_cols(fread(file_biov1_LHTBI),fread(file_biov1_pi),fread(file_biov1_pp)) %>%
    dplyr::filter(cppday >= 0) %>%
    dplyr::rename(bin=Name)
rast_biov1 <- var_to_rast(df_biov1 %>% dplyr::select(bin,cppday), ext=c(xlims,ylims))
rast_svgpm <- get_op_pp(file_svgpm,rast_biov1,xlims,ylims)
load(file_biov2_month_MPI_U)
rast_biov2_month_MPI_U = var_to_rast(rda_output %>% dplyr::select(bin, PP), ext=c(xlims,ylims))

df <- data.frame(BIOv2=getValues(rast_biov2_month_MPI_U),
                 BIOv1=getValues(rast_biov1),
                 VGPM=getValues(rast_svgpm)) %>%
    tidyr::drop_na()




# ggsave(filename="model_comparison.png",
#        plot=final,
#        dpi=150,
#        units="px",
#        height=1000,
#        width=2000)
