# Script for running the coherent DVMP generator

import sys
import os
from ROOT import TLorentzVector,TH1F
import numpy as np
from tqdm import tqdm 

sys.path.append("src")

from constants import *
from runcard_parser import *
from dvmp_tree import *
from misc import *
from photoproduction_helper import *
from electroproduction_helper import *
from dvmp_cross_section import *
from acceptance_and_smearing import *

import argparse

def main(args):
    
    # Ensure the current working directory is the base directory of the repository
    current_dir = os.getcwd()
    if not current_dir.endswith('SoLIDThresholdPhysics'):
        raise ValueError("Please run this script from the 'SoLIDThresholdPhysics' directory.")

    # Load parameters from runcard
    parameters = parse_runcard(args)
        
    # Get detector system
    set_detector(parameters)
    
    # Create output TFile and TTree and TTree for accepted events
    tfile, ttree, ttreeAcc = create_dvmp_ttree(parameters)
    
    # Calculate the integrated luminosity
    integrated_luminosity = get_integrated_luminosity(parameters)
    
    # Print out event generation info
    print_dvmp_start(parameters, tfile.GetName(),integrated_luminosity)

    
    # Event parameters
    target_mass = (mP if parameters.get("target_type")=="p" else mD)
    beam_energy = float(parameters.get("beam_energy"))     # GeV
    beam_current = float(parameters.get("beam_current"))   # Microamps
    target_length = float(parameters.get("target_length")) # Centimeters
    target_type = parameters.get("target_type")            # p/d
    if target_type not in ["p","d"]:
        raise ValueError("'target_type' must be either 'p' or 'd'.")
    days = float(parameters.get("days"))                   # Number of days
    process = parameters.get("process")
    detector = parameters.get("detector")
    model_type = parameters.get("model_type")              # '23g' or 'PomeronLQCD'
    num_events = int(parameters.get("num_events"))
    output_file_location = parameters.get("output_file_location")
    output_file_prefix = parameters.get("output_file_prefix")
    photon_energy_min = float(parameters.get("photon_energy_min")) # GeV
    photon_energy_max = float(parameters.get("photon_energy_max")) # GeV
    Q2_max = float(parameters.get("Q2_max")) # GeV^2
    t_min = float(parameters.get("t_min"))   # GeV^2

    rho,X0,_ = get_target_info(parameters)
    
    # Load DVMP production class
    DVMP = dvmpProduction(target_type, model_type)
    
    #######################################################################################
    # Event Generator
    #######################################################################################
    # Event Loop
    success=0
    for evt in tqdm(range(num_events)):
        
        
        
        # Set necessary variables to 0
        weight[0]=0.0
        decay_weight[0]=0.0
        psf[0]=0.0
        flux[0]=0.0
        acc_eOut[0]=0.0
        acc_hOut[0]=0.0
        acc_ePlus[0]=0.0
        acc_eMinus[0]=0.0
        Q2_[0]=0.0
        smear_Q2_[0]=0.0
        gammaE_[0]=0.0
        smear_gammaE_[0]=0.0
        t_[0]=0.0
        smear_t_[0]=0.0
        m_vm_[0]=0.0
        smear_m_vm_[0]=0.0
        jacobian[0] = 0.0
        smear_jacobian[0] = 0.0
        smear_eOut.SetPxPyPzE(0,0,0,0)
        smear_ePlus.SetPxPyPzE(0,0,0,0)
        smear_eMinus.SetPxPyPzE(0,0,0,0)
        smear_hOut.SetPxPyPzE(0,0,0,0)
        
        # Set initial TLorentzVectors
        eIn.SetPxPyPzE(0,0,np.sqrt(beam_energy**2-mE**2),beam_energy)
        hIn.SetXYZM(0,0,0,target_mass)
        
        # Real Photon/Virtual Photon generation
        if process == "photoproduction":
            # Do Bremmstrahlung & EPA
            flux[0], gammaE = N_Bremmstrahlung(photon_energy_min,photon_energy_max,beam_energy,target_length,X0)
            flux[0]=flux[0]+N_EquivalentPhotonApproximation(gammaE,beam_energy,Q2_max)
            q.SetPxPyPzE(0,0,gammaE,gammaE)
    
        elif process == "electroproduction":
            # Do virtual photon
            flux[0], temp = VirtualPhoton(photon_energy_min,photon_energy_max,beam_energy,target_mass,eIn,eOut,hIn)
            if(flux[0]==0):
                continue
            q.SetPxPyPzE(temp.Px(),temp.Py(),temp.Pz(),temp.E())
            gammaE=q.E()
    
    
        # Calculate Q2
        Q2 = -q*q
        
        # Set scattered electron
        eOut.SetPxPyPzE(eIn.Px()-q.Px(),
                        eIn.Py()-q.Py(),
                        eIn.Pz()-q.Pz(),
                        eIn.E()-q.E())
        
        # Check the W threshold
        W = (q+hIn).M()
        if W < mJpsi+target_mass:
            continue
            
        # Decay (gamma + p,d) --> (J/psi + p,d)
        GenPhase.SetDecay(q+hIn,2,np.array([mJpsi,target_mass]))
        GenPhase.Generate()
        temp1, temp2 = GenPhase.GetDecay(0),GenPhase.GetDecay(1)
        VM.SetPxPyPzE(temp1.Px(),temp1.Py(),temp1.Pz(),temp1.E())
        hOut.SetPxPyPzE(temp2.Px(),temp2.Py(),temp2.Pz(),temp2.E())
        psf[0] = 4 * np.pi
        
        # Get event parameters
        t = (hIn-hOut).M2()

        if(t < t_min):
            continue
    
        # Determine kinematics of first decay
        pCM_Initial = np.sqrt((W**2-hIn*hIn-q*q)**2-4*(hIn*hIn)*(q*q))/(2*W)
        pCM_Final = np.sqrt((W**2-hOut*hOut-VM*VM)**2-4*(hOut*hOut)*(VM*VM))/(2*W)
        cth = (np.sqrt(q*q+pCM_Initial**2)*np.sqrt(VM*VM+pCM_Final**2)-q*VM)/(pCM_Initial*pCM_Final)
    
        # Record jacobian
        jacobian[0] = 2*pCM_Initial*pCM_Final/(2*np.pi)
        
        # Get Cross Sectional Weight of Event --> weight[0] = dsigma_dcth 
        weight[0] = DVMP.dsigma([gammaE,t]) * jacobian[0]
        if(weight[0]<0):
            continue
        
        # Decay the J/Psi
        GenPhase.SetDecay(VM,2,np.array([mE,mE]))
        GenPhase.Generate()
        ep=GenPhase.GetDecay(0)
        em=GenPhase.GetDecay(1)
        ePlus.SetPxPyPzE(ep.Px(),ep.Py(),ep.Pz(),ep.E())
        eMinus.SetPxPyPzE(em.Px(),em.Py(),em.Pz(),em.E())
        
        # Save Kinematics
        W_[0] = W
        gammaE_[0]   = gammaE
        t_[0]        = t
        m_vm_[0]     = VM.M()
        Q2_[0]       = Q2
        
        # Calculate more kinematics
        y = (eIn-eOut).E()/eIn.E()
        Q2 = - (eIn-eOut).M2()
        gy = np.sqrt(Q2) / beam_energy
        eps = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy)
        R = (1.0 + Q2/2.164/(mJpsi**2))**(2.131)-1.0
        r = eps * R / (1.0 + eps * R)
        Ep = VM*hOut/mJpsi
        p = np.sqrt(Ep**2-target_mass**2)
        l = np.sqrt((mJpsi**2-ePlus*ePlus-eMinus*eMinus)**2-4*(ePlus*ePlus)*(eMinus*eMinus))/(2*mJpsi)
        cth_decay = (Ep * mJpsi/2 - ePlus*hOut)/(p*l)
        wth = 0.75*(1.0 + r + (1.0 - 3.0 * r) * cth_decay**2)
        branch = 0.05971
        decay_weight[0]=wth*branch
        
        # Get the acceptances of the final state particles
        acc_ePlus[0] = acc_e(ePlus)
        acc_eMinus[0]= acc_e(eMinus)
        acc_eOut[0]  = acc_e(eOut)
        acc_hOut[0]  = acc_h(hOut)
        
        # Smear outgoing scattered electron
        if acc_eOut[0]>0:
            temp = smear_particle(eOut,"e-",detector)
            smear_eOut.SetPx(temp.Px())
            smear_eOut.SetPy(temp.Py())
            smear_eOut.SetPz(temp.Pz())
            smear_eOut.SetE(temp.E())
        
        # Smear outgoing decay positron
        if acc_ePlus[0]>0:
            temp = smear_particle(ePlus,"e+",detector)
            smear_ePlus.SetPx(temp.Px())
            smear_ePlus.SetPy(temp.Py())
            smear_ePlus.SetPz(temp.Pz())
            smear_ePlus.SetE(temp.E())

        # Smear outgoing decay electron
        if acc_eMinus[0]>0:
            temp = smear_particle(eMinus,"e-",detector)
            smear_eMinus.SetPx(temp.Px())
            smear_eMinus.SetPy(temp.Py())
            smear_eMinus.SetPz(temp.Pz())
            smear_eMinus.SetE(temp.E())
        
        # Smear outgoing target
        if acc_hOut[0]>0:
            temp = smear_particle(hOut,target_type,detector)
            smear_hOut.SetPx(temp.Px())
            smear_hOut.SetPy(temp.Py())
            smear_hOut.SetPz(temp.Pz())
            smear_hOut.SetE(temp.E())
        
        if acc_ePlus[0]>0 and acc_eMinus[0]>0:
            temp = (smear_ePlus+smear_eMinus)
            smear_VM.SetPxPyPzE(temp.Px(),temp.Py(),temp.Pz(),temp.E())
            
        if acc_ePlus[0]>0 and acc_eMinus[0]>0 and acc_hOut[0]>0:
            temp = (smear_ePlus+smear_eMinus+smear_hOut-hIn)
            smear_q.SetPxPyPzE(temp.Px(),temp.Py(),temp.Pz(),temp.E())
    
        if acc_hOut[0]>0:
            smear_t_[0] = (hIn-smear_hOut).M2()
        
        smear_Q2_[0]=-smear_q*smear_q
        smear_gammaE_[0]=smear_q.E()
        smear_m_vm_[0]=(smear_VM).M()
        smear_W_[0]=(hIn+smear_q).M()
        # Record smear jacobian
        smear_pCM_Initial = np.sqrt((smear_W_[0]**2-hIn*hIn-smear_q*smear_q)**2-4*(hIn*hIn)*(smear_q*smear_q))/(2*smear_W_[0])
        smear_pCM_Final = np.sqrt((smear_W_[0]**2-smear_hOut*smear_hOut-smear_VM*smear_VM)**2-4*(smear_hOut*smear_hOut)*(smear_VM*smear_VM))/(2*smear_W_[0])
        smear_jacobian[0] = 2*smear_pCM_Initial*smear_pCM_Final/(2*np.pi)
        
        # Record if the event was successfully simulated
        if(weight[0]>0.0):
            success+=1
            ttree.Fill()
            if acc_ePlus[0]>0 and acc_eMinus[0]>0 and acc_hOut[0]>0:
                ttreeAcc.Fill()
    print(f"{success} successfully generated events out of {num_events} ({100*success/num_events:.2f}%)")
    ttree.Write()
    ttreeAcc.Write()
    # Save luminosity to the TFile
    h1_lumi = TH1F("integrated_luminosity","",1,0,1)
    h1_lumi.SetBinContent(1,integrated_luminosity)
    h1_lumi.Write()
    
    # Save Nevents to the TFile
    h1_nevents = TH1F("num_events","",1,0,1)
    h1_nevents.SetBinContent(1,num_events)
    h1_nevents.Write()
    
    # Save Beam Energy to the TFile
    h1_beamEnergy = TH1F("beam_energy","",1,0,1)
    h1_beamEnergy.SetBinContent(1,beam_energy)
    h1_beamEnergy.Write()
    
    # Save Days to the TFile
    h1_days = TH1F("days","",1,0,1)
    h1_days.SetBinContent(1,days)
    h1_days.Write()
    tfile.Close()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DVMP Generator Script")
    parser.add_argument("--runcard", type=str, required=True, help="Path to the runcard.card file.")
    parser.add_argument("--batch", type=int, default=0, required=False, help="Batch number (used for parallel running).")
    args = parser.parse_args()
    main(args)
    
    
    
# Load variables
######################################################################
# beam_energy = float(parameters.get("beam_energy"))     # GeV
# beam_current = float(parameters.get("beam_current"))   # Microamps
# target_length = float(parameters.get("target_length")) # Centimeters
# target_type = parameters.get("target_type")            # p/d
# if target_type not in ["p","d"]:
#     raise ValueError("'target_type' must be either 'p' or 'd'.")
# days = float(parameters.get("days"))                   # Number of days
# num_events = int(parameters.get("num_events"))
# output_file_location = parameters.get("output_file_location")
# output_file_prefix = parameters.get("output_file_prefix")
# photon_energy_min = float(parameters.get("photon_energy_min")) # GeV
# photon_energy_max = float(parameters.get("photon_energy_max")) # GeV
# Q2_max = float(parameters.get("Q2_max")) # GeV^2
# t_max = float(parameters.get("t_max"))   # GeV^2
