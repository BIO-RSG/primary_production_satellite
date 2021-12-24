
# function to load and format vgpm and cbpm from the ocean productivity website
get_op_pp <- function(f,new_rast,xlims,ylims) {
    rast <- raster(f)
    extent(rast) <- extent(c(xmin=-180,xmax=180,ymin=-90,ymax=90))
    rast <- crop(rast, extent(c(xlims,ylims)))
    rast[rast < 0] <- NA
    rast <- resample(rast,new_rast)
    return(rast)
}

# set extreme values of raster r to limits
set_lims <- function(r,lims) {
    r[r < lims[1]] <- lims[1]
    r[r > lims[2]] <- lims[2]
    return(r)
}

# make table of correlation coefficients with %diff values
cor_table <- function(df,var) {
    cor1 <- cor.test(df[,var],df$chl)
    cor2 <- cor.test(df[,var],df$lat)
    cor3 <- cor.test(df[,var],df$bathy)
    cor_df <- data.frame(Variable = c("chla","latitude","bathymetry"),
                         r = c(cor1$estimate,cor2$estimate,cor3$estimate),
                         p = c(cor1$p.value,cor2$p.value,cor3$p.value))
    cor_df$r <- round(cor_df$r,4)
    cor_df$p <- round(cor_df$p,4)
    p <- ggplot() +
        geom_table_npc(data=cor_df, aes(npcx=0.5, npcy=0.5), label=list(cor_df)) +
        theme_void()
    return(p)
}

# calculate percent difference between 2 matrices x and y
pdiff <- function(x,y) {
    return(100*(x-y)/y)
}

# make percent difference maps
pdiff_map <- function(r,orig_rast,xlims,col_limits_pdiff,rchl,rlat,rbathy) {
    diff <- orig_rast
    values(diff) <- pdiff(as.matrix(r),as.matrix(orig_rast))
    full_df <- data.frame(pdiff=getValues(diff),
                          chl=getValues(rchl),
                          lat=getValues(rlat),
                          bathy=getValues(rbathy))
    cor_plot <- cor_table(full_df,"pdiff")
    mdens <- dens_plot(diff, xlims=xlims)
    diff <- set_lims(diff,col_limits_pdiff)
    return(list(diff=diff, cor_plot=cor_plot, mdens=mdens))
}

arrange_maps <- function(rast_biov1, xlims, rlist, rnames, col_limits_pdiff, rast_chl, rast_lat, rast_bathy, ...) {
    pm <- lapply(rlist,
                 FUN=pdiff_map,
                 orig_rast=rast_biov1,
                 xlims=xlims,
                 col_limits_pdiff=col_limits_pdiff,
                 rchl=rast_chl,
                 rlat=rast_lat,
                 rbathy=rast_bathy)
    pp_rasts <- lapply(pm, "[[", 1)
    names(pp_rasts) <- rnames
    if (length(pp_rasts)==1) {
        pp_rasts <- pp_rasts[[1]]
    } else {
        pp_rasts <- raster::stack(pp_rasts)
    }
    cor_plots <- lapply(pm, "[[", 2)
    mrast <- make_raster_map(pp_rasts,
                             col_limits=col_limits_pdiff,
                             cm = colorRampPalette(c("blue","white","red"))(100),
                             ...) +
        theme(plot.margin=unit(c(0,0,5,0),"mm"))
    mdens <- lapply(pm, "[[", 3)
    final <- wrap_plots(mdens,nrow=1) / mrast / wrap_plots(cor_plots,nrow=1) + plot_layout(heights=c(2,14,3))
    return(final)
}

# make a percent difference density plot
dens_plot <- function(r,xlims) {
    vals <- getValues(r)
    vmean <- mean(vals,na.rm=TRUE)
    vsd <- sd(vals,na.rm=TRUE)
    xlims <- c(vmean-2*vsd,vmean+2*vsd)
    ggplot(data.frame(percent_difference=vals)) +
        geom_density(aes(percent_difference),fill="blue",alpha=0.4) +
        theme_bw() +
        scale_x_continuous(limits=xlims, position="top") +
        # ggtitle("Percent difference") +
        theme(axis.title=element_blank(),
              axis.text.y=element_blank(),
              axis.text.x=element_text(angle=45,hjust=0),
              plot.margin=unit(c(0,0,0,1), "mm"),
              axis.ticks.y=element_blank())
}
