# Stephanie.Clay@dfo-mpo.gc.ca
# 2023-04-18

# Generating an image for claysa_projects.pptx

library(oceancolouR)
library(dplyr)
library(ggplot2)
library(data.table)
library(raster)
library(patchwork)
library(scico)
source("markdown_reports/PP_comparison_rmd_functions.R")

path <- "model_comparisons/"
date_bio <- "201807"
date_op <- "2018182"
xlims <- c(-95,-42)
ylims <- c(39,82)
col_limits_pp <- c(0,3000)
col_limits_pdiff <- c(-100,100)
cm <- colorRampPalette(rev(scico(100, palette = 'roma')))(100)

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


mrm <- function(rast,title,cl,bh=0.6,bw=8) {
    map <- make_raster_map(rast, xlim=xlims, ylim=ylims, cm=cm, col_limits=cl,
                           set_extremes=TRUE, show_legend=TRUE, title=title, map_alpha=1) +
        guides(fill=guide_colorbar(title.hjust = 0,
                                   ticks.colour = "black",
                                   barheight = unit(bh, "cm"),
                                   barwidth = unit(bw, "cm"),
                                   frame.colour = "black")) &
            theme(legend.direction="horizontal",
                  plot.title=element_text(size=18))
    return(map)
}

p1 <- mrm(rast_biov2_month_MPI_U, "BIOv2", col_limits_pp)
p2 <- mrm(rast_biov1, "BIOv1", col_limits_pp)
p3 <- mrm(rast_svgpm, "VGPM", col_limits_pp)

pdiff <- function(v1,v2) {
    df <- (v1-v2)/v2 * 100
    df <- df[is.finite(df)]
    ggplot(data.frame(percent_difference=df)) +
        geom_vline(xintercept=0, linetype="dashed", alpha=0.6) +
        geom_density(aes(percent_difference),fill="blue",alpha=0.4) +
        theme_bw() +
        scale_x_continuous(limits=col_limits_pdiff) +
        scale_y_continuous(expand=c(0,0)) +
        theme(axis.title=element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(size=12,angle=45,hjust=1),
              axis.ticks.y=element_blank())
}

vals1 <- getValues(rast_biov2_month_MPI_U)
vals2 <- getValues(rast_biov1)
vals3 <- getValues(rast_svgpm)
p4 <- pdiff(vals1,vals2)
p5 <- pdiff(vals1,vals3)


# BIOv2 MPI_U, and % difference from BIOv1 and VGPM

final <- wrap_plots(list(p1,p2,p3,guide_area(),p4,p5)) + plot_layout(guides="collect",byrow=TRUE, heights=c(6,1))

# ggsave(filename="model_comparison.png",
#        plot=final,
#        dpi=150,
#        units="px",
#        height=1000,
#        width=2000)


stop()

pbiov2 <- p1 +
    guides(fill=guide_colorbar(title.hjust = 0,
                               ticks.colour = "black",
                               barheight = unit(8, "cm"),
                               barwidth = unit(0.6, "cm"),
                               frame.colour = "black")) +
    labs(fill = "Primary<br>production<br>(mgC /m^2 /day)<br>") +
    theme(legend.direction="vertical",
          legend.title=element_markdown(size=14,hjust=0),
          axis.text=element_text(size=12),
          legend.margin=margin(0,0,0,5),
          legend.text=element_text(size=12),
          plot.title=element_blank())

ggsave(filename="biov2.png",
       plot=pbiov2,
       dpi=150,
       units="px",
       height=800,
       width=900)



stop()

chl <- mrm(var_to_rast(df_biov1 %>% dplyr::select(bin,chlor_a) %>% dplyr::mutate(chlor_a=log10(chlor_a)), ext=c(xlims,ylims)), "log10(Chla)", c(-1.5,0.8))
pmb <- mrm(var_to_rast(df_biov1 %>% dplyr::select(bin,PmB), ext=c(xlims,ylims)), "PmB", c(1,7))
alpha <- mrm(var_to_rast(df_biov1 %>% dplyr::select(bin,alphaB), ext=c(xlims,ylims)), "alphaB", c(0.001,0.3))

final2 <- chl+pmb+alpha & theme(legend.direction="horizontal",legend.position="bottom")
ggsave(filename="model_comparison_inputs.png",
       plot=final2,
       dpi=150,
       units="px",
       height=900,
       width=2000)
