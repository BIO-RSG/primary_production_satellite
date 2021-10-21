# watt to quanta (lambda in nm)
Watt_perm_2_to_microMol_per_m2_per_s = function(Ed, lambda){
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed * lambda * 1e-9 / h / c / Na / 1e-6)
}

# quanta to watt (lambda in nm)
microMol_per_m2_per_s_to_watt_per_m2 = function(Ed, lambda){
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed / lambda / 1e-9 * h * c * Na * 1e-6)
}

# lambda in nm
Watt_per_m2_to_Mol_per_m2_perhour <- function(Ed, lambda) {
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed * lambda * 1e-9 * 3600 / h / c / Na)
}

# lambda in nm
Watt_per_m2_to_Mol_per_m2_persec <- function(Ed, lambda) {
    h = 6.626e-34 # in J.s
    c = 2.998e8 # in m.s-1
    Na = 6.022e23 # Avogadro number
    return(Ed * lambda * 1e-9 / h / c / Na)
}
