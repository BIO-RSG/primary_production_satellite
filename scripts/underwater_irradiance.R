underwater_irradiance <- function(E0dir, E0dif, day_hrs, nstep, chlpix, lambda, pc1f, pc2f, ratef, fracys, Sy, asw, zendw, alphaB, PBm, BW) {

    PP <- matrix(0,23,nstep)

    for (it in day_hrs) # loop on the time of day (for daylight zendR)
    {

        for (iz in 2:nstep) # loop on depth, we start at the second step as the first one is Eqdw
        {

            chlz = chlpix[iz]

            # To compute the underwater light field, we need to know total absorption and total backscattering
            # a = aw + aphy + acdom, bb = bbw + bbp, this as to be define for each wavelength

            # backscattering is modelled according to Loisel and Morel 1998
            bbtilda = (0.78+0.42*(-log10(chlz)))*0.01 # bbtilda et le rapport entre bb et b.
            if (bbtilda < 0.0005) {bbtilda = 0.0005
            } else if (bbtilda > 0.01) {bbtilda = 0.01}
            BC660 = 0.407*chlz^0.795 # ici, bc corresponds to b_chl a 660 nm (scattering by chlorophyll-a)
            BC = BC660*(660/lambda)^(-log10(chlz)) # BC spectral dependence depends on chlz,
            # this means that at chl = 1 tere is no spectral dependence.This is weird and has to be checked
            BC[BC < 0] = 0
            bbt = BC*bbtilda + BW*0.5 # Here, we get total backscattering.

            #phytoplankton absorption ### !!!!! THIS HAS TO BE UPDATED.
            aphy = pc1f*(1-exp(-ratef*chlz)) + pc2f*chlz

            # Yellow substances and detritus absorption is derived from chl
            ays = fracys * aphy[44] * exp(-Sy * (lambda -443))

            # Now we can compute total absorption
            at = asw + ays + aphy

            # Here we can start to compute attenuation at each depth for each wavelength
            Kddir = (at + bbt) / cos(zendw[it])
            Kddif = (at + bbt) / 0.83

            E0dir[it,iz,] = E0dir[it,(iz-1),] * exp(-Kddir*0.5)
            E0dif[it,iz,] = E0dif[it,(iz-1),] * exp(-Kddif*0.5)

            # Here we compute the photosynthetic action spectrum according to Kyewangala et al. 1997
            apbar = mean(aphy)
            alphaBc = alphaB * aphy / apbar

            #PIz = 1/cos(zendw[it]) * sum(alphaBc *  E0dir[it,iz,]) + 1.20 * sum(alphaBc * E0dif[it,iz,])
            PIz = sum(alphaBc * (E0dir[it,iz,] + E0dif[it,iz,]))
            PP[it,iz] = chlz * PIz/sqrt(1 + (PIz/PBm)^2)

            # Here we need to compute PP at 0 as the loop starts at 2:
            if (iz == 2)
            {
                #PIz = 1/cos(zendw[it]) * sum(alphaBc *  E0dir[it,1,]) + 1.20 * sum(alphaBc * E0dif[it,1,])
                PIz = sum(alphaBc * (E0dir[it,1,] + E0dif[it,1,]))
                PP[it,1] = chlz * PIz/sqrt(1 + (PIz/PBm)^2)
            }

        } # end loop on depth

    } # end loop on time of day (xhr)

    return(PP)

}
