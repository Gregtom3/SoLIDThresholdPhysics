import numpy as np

def print_dvmp_start(parameters,fname,integrated_luminosity):
    # Load variables
    ######################################################################
    beam_energy = float(parameters.get("beam_energy"))     # GeV
    beam_current = float(parameters.get("beam_current"))   # Microamps
    target_length = float(parameters.get("target_length")) # Centimeters
    target_type = parameters.get("target_type")            # p/d
    days = float(parameters.get("days"))                   # Number of days
    process = parameters.get("process")
    num_events = int(parameters.get("num_events"))
    photon_energy_min = float(parameters.get("photon_energy_min")) # GeV
    photon_energy_max = float(parameters.get("photon_energy_max")) # GeV
    Q2_max = float(parameters.get("Q2_max")) # GeV^2
    t_min = float(parameters.get("t_min"))   # GeV^2
    
    print("--------------------------------------------------------")
    print("Filename ==>",fname)
    print("--------------------------------------------------------")
    print("Target Type: ",target_type)
    print("Target Length: ", target_length,"cm")
    if(process=="photoproduction"):
        print("Beginning Monte Carlo Simulation of e({:.2f} GeV)+{} --> gamma + {} --> J/psi(e-e+) + {}'".format(beam_energy,target_type,target_type,target_type))
        print(" *** PHOTO-PRODUCTION *** ")
    else:
        print("Beginning Monte Carlo Simulation of e({:.2f} GeV)+{} --> gamma* + {} --> e' + J/psi(e-e+) + {}'".format(beam_energy,target_type,target_type,target_type))
        print(" *** ELECTRO-PRODUCTION *** ")
    print("--------------------------------------------------------")
    print("Event Kinematics")
    print(photon_energy_min, "< Egamma <", photon_energy_max,"[GeV]")
    if(process=="electroproduction"):
        print("EPA Maximum Q2 =",Q2_max,"[GeV^2]")
    print("--------------------------------------------------------")
    print("Number of Events =",num_events)
    print("Luminosity =", np.round(integrated_luminosity/3600.0/24/days,2),"(events/nb/s)")
    print("Integrated Luminosity =",np.round(integrated_luminosity,2),"(events/nb)")
    print("--------------------------------------------------------\n\n")
    
    

def print_bh_start(parameters,fname,integrated_luminosity):
    # Load variables
    ######################################################################
    beam_energy = float(parameters.get("beam_energy"))     # GeV
    beam_current = float(parameters.get("beam_current"))   # Microamps
    target_length = float(parameters.get("target_length")) # Centimeters
    target_type = parameters.get("target_type")            # p/d
    days = float(parameters.get("days"))                   # Number of days
    num_events = int(parameters.get("num_events"))
    photon_energy_min = float(parameters.get("photon_energy_min")) # GeV
    photon_energy_max = float(parameters.get("photon_energy_max")) # GeV
    Mll2_min = float(parameters.get("Mll2_min")) # GeV^2
    Mll2_max = float(parameters.get("Mll2_max")) # GeV^2
    t_min = float(parameters.get("t_min"))   # GeV^2
    t_max = float(parameters.get("t_max"))   # GeV^2
    
    print("--------------------------------------------------------")
    print("Filename ==>",fname)
    print("--------------------------------------------------------")
    print("Target Type: ",target_type)
    print("Target Length: ", target_length,"cm")
    print("Beginning Monte Carlo Simulation of e({:.2f} GeV)+{} --> gamma + {} --> e- + e+ + {}'".format(beam_energy,target_type,target_type,target_type))
    print(" *** COHERENT BETHE-HEITLER *** ")
    print("--------------------------------------------------------")
    print("Event Kinematics")
    print(photon_energy_min, "< Egamma <", photon_energy_max,"[GeV]")
    print(np.abs(t_max), "< |t| <", np.abs(t_min),"[GeV^2]")
    print(Mll2_min, "< Mll2 <", Mll2_max,"[GeV^2]")
    print("--------------------------------------------------------")
    print("Number of Events =",num_events)
    print("Luminosity =", np.round(integrated_luminosity/3600.0/24/days,2),"(events/nb/s)")
    print("Integrated Luminosity =",np.round(integrated_luminosity,2),"(events/nb)")
    print("--------------------------------------------------------\n\n")