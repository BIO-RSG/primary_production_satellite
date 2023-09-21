# plot a simple PI curve for demonstration
# here, alphab and pmb=1
library(ggplot2)
library(ggpp)
x <- seq(0,4,by=0.1)
y <- x/sqrt(1+x^2)
fr <- paste("α == {frac('ΔP', 'ΔI')}")
p <- ggplot(data.frame(x=x,y=y)) +
    geom_line(aes(x,y)) +
    geom_abline(slope=1,intercept=0,color="red",linetype="dashed") +
    geom_hline(yintercept=1,color="red",linetype="dashed") +
    geom_line(data=data.frame(x=c(0.5,0.5),y=c(0,0.5)),aes(x,y),color="red") +
    geom_line(data=data.frame(x=c(0,0.5),y=c(0,0.5)),aes(x,y),color="red") +
    geom_line(data=data.frame(x=c(0,0.5),y=c(0,0)),aes(x,y),color="red") +
    scale_x_continuous(limits=c(0,max(x)),expand=c(0,0)) +
    scale_y_continuous(limits=c(-0.001,1.1),expand=c(0,0),breaks=seq(0,1,by=0.2),labels=c(0,rep(" ",4),"PmB")) +
    theme_bw() +
    labs(x="Irradiance",y="Photosynthetic Rate") +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_text(size=16),
          axis.title.x=element_text(size=18),
          axis.title.y=element_text(size=18,margin=margin(0,-20,0,0))) +
    annotate(geom='text', x=1, y=0.25, label=fr, parse=TRUE, color='red', size=8)

ggsave(filename="PI_curve_diagram.png",
       plot=p,
       dpi=150,
       units="px",
       width=800,
       height=600)

