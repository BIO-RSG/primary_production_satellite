# watt to quanta
Watt_perm_2_to_microMol_per_m2_per_s = function(Ed, lambda){
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed * lambda * 1e-9 / h / c / Na / 1e-6)
}

# quanta to watt
microMol_per_m2_per_s_to_watt_per_m2 = function(Ed, lambda){
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed / lambda / 1e-9 * h * c * Na * 1e-6)
}
