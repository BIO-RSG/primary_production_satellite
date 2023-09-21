# returns total, direct, and diffuse irradiance in w/m2/nm
surface_irr <- function(zendR, lam, oza, ag, aw, sco3, p0, wv, rh, am, wsm, ws, pxlvis, Fo) {

    rad = 180/pi

    # Ed(0+,lambda)
    Ed = matrix(0,23,301)
    Eddir = matrix(0,23,301)
    Eddif = matrix(0,23,301)

    # Eq(0-)
    Eqdirw = matrix(0,23,301) # direct
    Eqdifw = matrix(0,23,301) # diffus

    for (i in 1:23)
    {

        if (zendR[i] < 90.)
        {

            # above water spectral irradiance
            irr.gc = atmodd_v2(rod=0.0, ros=0.0, rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag,
                               aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                               ws = ws, vis = pxlvis, Fo = Fo)
            Ed[i,] = irr.gc$Ed
            Eddir[i,] = irr.gc$Edir
            Eddif[i,] = irr.gc$Eddif

            # Under water spectral irradiance
            res.sfc <- sfcrfl(rad,zendR[i],ws)
            irrm.gc = atmodd_v2(rod=res.sfc$rod, ros=res.sfc$ros, rad = rad, lam = lam, theta = zendR[i], oza = oza, ag = ag,
                                aw = aw, sco3 = sco3, p=p0, wv = wv, rh = rh, am = am, wsm = wsm,
                                ws = ws, vis = pxlvis, Fo = Fo)
            Eqdirw[i,] = irrm.gc$Edir
            Eqdifw[i,] = irrm.gc$Edif

        }
    }

    return(list(Ed=Ed, Eddir=Eddir, Eddif=Eddif, Eqdirw=Eqdirw, Eqdifw=Eqdifw))

}
