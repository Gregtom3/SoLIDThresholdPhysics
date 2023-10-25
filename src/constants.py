from ROOT import TGenPhaseSpace, TVector3

# Mass of Electron
# ---------------------------------------------
mE    = 0.000510999

# Mass of JPsi
# ---------------------------------------------
mJpsi = 3.0969

# Mass of Proton
# ---------------------------------------------
mP    = 0.938272

# Mass of Deuteron
# ---------------------------------------------
mD    = 1.8761358

# Fine Structure Constant
# ---------------------------------------------
alpha_em = 1.0/137.036

# Phase Space Decay
# ---------------------------------------------
GenPhase = TGenPhaseSpace()

# z-hat
# ---------------------------------------------
Zhat=TVector3(0,0,1)

# Return target configuration information
def get_target_info(parameters):
    
    target_type = parameters.get("target_type")
    
    if(target_type=="p"):
        rho     = 0.071 # Density of Target g/cm3
        X0      = 890   # Radiation Length [cm]
        factor  = 1
    elif(target_type=="d"):
        rho     = 0.169 # Density of Target g/cm3
        X0      = 769.1 # Radiation Length [cm]
        factor = 1/2
    else:
        raise ValueError("'target_type' must be either 'p' or 'd'.")
        
    return rho,X0,factor

def get_integrated_luminosity(parameters):
    
    rho,X0,factor = get_target_info(parameters)
    
    beam_current = float(parameters.get("beam_current"))
    days         = float(parameters.get("days"))
    target_length = float(parameters.get("target_length"))
    
    integrated_luminosity = factor * beam_current * 1e-6/(1.6e-19) * target_length * rho * 6.02e23 * 1.0e-24 * 1.0e-9 * 3600.0 * 24 * days 
    
    return integrated_luminosity