# Stephanie.Clay@dfo-mpo.gc.ca
# 2021-12-14

# Make figures for PP_summary.

library(oceancolouR)
library(stringr)
library(data.table)
library(raster)
library(ggplot2)
library(patchwork)
source("markdown_reports/PP_comparison_rmd_functions.R")

# date system used for BIO files
dates_bio_monthly <- c("201210","201501","201704","201807")
dates_vgpm <- c("2012275","2015001","2017091","2018182")
dates_bio_8day <- list(d1=c("2012273","2012281","2012289","2012297"),
                       d2=c("2015001","2015009","2015017","2015025"),
                       d3=c("2017089","2017097","2017105","2017113"),
                       d4=c("2018185","2018193","2018201","2018209"))
# daily dates matched to the second 8day period above
dates_bio_daily_8day <- list(d1=paste0("2012",str_pad(281:288,width=3,side="left",pad="0")),
                             d2=paste0("2015",str_pad(9:16,width=3,side="left",pad="0")),
                             d3=paste0("2017",str_pad(97:104,width=3,side="left",pad="0")),
                             d4=paste0("2018",str_pad(193:200,width=3,side="left",pad="0")))
# daily dates matched to the month above
dates_bio_daily_month <- list(d1=paste0("2012",get_doys(2012,month=10)),
                              d2=paste0("2015",get_doys(2015,month=1)),
                              d3=paste0("2017",get_doys(2017,month=4)),
                              d4=paste0("2018",get_doys(2018,month=7)))

# indices of the yearmonth to use
which_dates <- 1:4
# x-axis limits of density plot for figure 2
dens_lims <- list(c(10,1e4),c(5,5000),c(30,1e4),c(80,1e4))

path <- "model_comparisons"

xlims <- c(-95,-42)
ylims <- c(39,82)



#*******************************************************************************
# FIGURE 1

pmap <- function(new,old,xlims) {
    x <- as.matrix(new)
    y <- as.matrix(old)
    diff <- old
    values(diff) <- 100*(x-y)/y
    mdens <- dens_plot(diff, xlims=xlims)
    return(diff)
}

for (ind in 1:length(which_dates)) {

    which_date <- which_dates[ind]

    output_name <- paste0("images/summary_fig01_",dates_bio_monthly[which_date],".png")

    biov1 <- lapply(paste0(path,"/BIO_v1_Carla/",dates_bio_monthly[which_date],"_",
                           c("LHTBI.csv", "pi-20141113.csv", "pp.csv")),
                    fread) %>%
        do.call(what=dplyr::bind_cols) %>%
        dplyr::filter(cppday >= 0) %>%
        dplyr::rename(bin=Name) %>%
        dplyr::select(bin,cppday) %>%
        var_to_rast(ext=c(xlims,ylims))
    vgpm <- get_op_pp(paste0(path,"/VGPM/vgpm.viirs.",dates_vgpm[which_date],".hdf"),
                      biov1,xlims,ylims)
    biov2_same <- get(load(paste0(path,"/BIOv2_vs_BIOv1/BIOv2_BIOv1-input_variablePI_NU_BIO_NWA_monthly_",dates_bio_monthly[which_date],".rda"))) %>%
        dplyr::select(bin, PP) %>%
        var_to_rast(ext=c(xlims,ylims))
    biov2_vpinu <- get(load(paste0(path,"/BIOv2-monthly_vs_BIOv1-monthly/BIOv2_variablePI_NU_BIO_NWA_monthly_",dates_bio_monthly[which_date],".rda"))) %>%
        dplyr::select(bin, PP) %>%
        var_to_rast(ext=c(xlims,ylims))
    biov2_mpinu <- get(load(paste0(path,"/BIOv2-monthly_vs_BIOv1-monthly/BIOv2_medianPI_U_BIO_NWA_monthly_",dates_bio_monthly[which_date],".rda"))) %>%
        dplyr::select(bin, PP) %>%
        var_to_rast(ext=c(xlims,ylims))
    biov2_list <- list(biov2_same, biov2_vpinu, biov2_mpinu)

    # # map the differences in chl between "same" and "vpinu", and the %difference
    # # in pp from these two options (it's relatively small)
    # r1c <- var_to_rast(biov2_same%>%dplyr::select(bin,chl)%>%dplyr::mutate(chl=log10(chl)),ext=c(xlims,ylims))
    # r2c <- var_to_rast(biov2_vpinu%>%dplyr::select(bin,chl)%>%dplyr::mutate(chl=log10(chl)),ext=c(xlims,ylims))
    # r1p <- var_to_rast(biov2_same%>%dplyr::select(bin,PP),ext=c(xlims,ylims))
    # r2p <- var_to_rast(biov2_vpinu%>%dplyr::select(bin,PP),ext=c(xlims,ylims))
    # make_raster_map(raster::stack(list(r1c,r2c)), col_limits=c(-1.5,1.5),set_extremes=TRUE)
    # rdiff <- 100*(r2p-r1p)/r1p
    # make_raster_map(rdiff, col_limits=c(-5,5),set_extremes=TRUE)

    pm_biov1 <- lapply(biov2_list, FUN=pmap, old=biov1, xlims=xlims)
    pm_vgpm <- lapply(biov2_list, FUN=pmap, old=vgpm, xlims=xlims)
    names(pm_biov1) <- names(pm_vgpm) <- c("biov1_input","new_input_VPI_NU","new_input_MPI_U")
    
    grid_biov1 <- make_raster_map(raster::stack(pm_biov1),
                                  col_limits=c(-100,100),
                                  set_extremes=TRUE,
                                  cm = colorRampPalette(c("blue","white","red"))(100)) +
        labs(y="BIOv1 % difference") +
        theme(plot.margin=unit(c(0,0,5,0),"mm"),
              axis.title.y=element_text(size=20, angle=90))
    mdens_vgpm <- wrap_plots(lapply(pm_biov1, FUN=dens_plot, xlims=xlims), nrow=1)

    grid_vgpm <- make_raster_map(raster::stack(pm_vgpm),
                                 col_limits=c(-100,100),
                                 set_extremes=TRUE,
                                 cm = colorRampPalette(c("blue","white","red"))(100)) +
        labs(y="VGPM % difference") +
        theme(plot.margin=unit(c(0,0,5,0),"mm"),
              axis.title.y=element_text(size=20, angle=90))
    mdens_vgpm <- wrap_plots(lapply(pm_vgpm, FUN=dens_plot, xlims=xlims), nrow=1)

    png(filename=output_name, width=800, height=800)
    print(mdens_vgpm / grid_biov1 / mdens_vgpm / grid_vgpm + plot_layout(heights=c(2,14,2,14)))
    dev.off()

}




stop()
#*******************************************************************************
# FIGURE 2

for (ind in 1:length(which_dates)) {
    
    which_date <- which_dates[ind]
    
    output_name <- paste0("images/summary_fig02_",dates_bio_monthly[which_date],".png")
    
    biov2_monthly <- get(load(paste0(path,"/BIOv2-monthly-8day-daily/BIOv2_medianPI_U_BIO_NWA_monthly_",dates_bio_monthly[which_date],".rda")))$PP
    
    biov2_8day <- lapply(paste0(path,"/BIOv2-monthly-8day-daily/BIOv2_medianPI_U_BIO_NWA_8day_",dates_bio_8day[[which_date]],".rda"),
                         function(x) get(load(x))$PP) %>%
        do.call(what=c)
    
    biov2_daily <- lapply(paste0(path,"/BIOv2-monthly-8day-daily/BIOv2_medianPI_U_BIO_NWA_daily_",dates_bio_daily_month[[which_date]],".rda"),
                         function(x) get(load(x))$PP) %>%
        do.call(what=c)
    
    pp_df <- data.frame(Input=c(rep("monthly",length(biov2_monthly)),
                                rep("8day",length(biov2_8day)),
                                rep("daily",length(biov2_daily))),
                        PP=c(biov2_monthly, biov2_8day, biov2_daily))
    
    p <- ggplot(pp_df) +
        geom_density(aes(PP, fill=Input), alpha=0.4) +
        scale_x_log10(limits=dens_lims[[ind]]) +
        labs(x="Primary production (daily average, mg C / m^2 / day)") +
        theme_bw() +
        theme(axis.title.y=element_blank(),
              axis.text.x=element_text(size=14),
              axis.title.x=element_text(size=16),
              legend.text=element_text(size=14),
              legend.title=element_text(size=16),
              legend.position=c(0.9, 0.7),
              legend.key.size = unit(1.5, 'cm'))
    
    png(filename=output_name, width=800, height=400)
    print(p)
    dev.off()
    
}

