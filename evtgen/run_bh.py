# Script for running the Bethe-Heitler Event Generator

import sys
import os
from ROOT import TLorentzVector,TH1F
import numpy as np
import random
from tqdm import tqdm 

sys.path.append("src")

from constants import *
from runcard_parser import *
from bh_tree import *
from misc import *
from bh_cross_section import *
from photoproduction_helper import *
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
    
    # Create output TFile and TTree
    tfile, ttree, ttreeAcc = create_bh_ttree(parameters)
    
    # Calculate the integrated luminosity
    integrated_luminosity = get_integrated_luminosity(parameters)
    
    # Print out event generation info
    print_bh_start(parameters, tfile.GetName(),integrated_luminosity)
    
    # Event parameters
    target_mass = (mP if parameters.get("target_type")=="p" else mD)
    beam_energy = float(parameters.get("beam_energy"))     # GeV
    beam_current = float(parameters.get("beam_current"))   # Microamps
    target_length = float(parameters.get("target_length")) # Centimeters
    target_type = parameters.get("target_type")            # p/d
    if target_type not in ["p","d"]:
        raise ValueError("'target_type' must be either 'p' or 'd'.")
    days = float(parameters.get("days"))                   # Number of days
    detector = parameters.get("detector")
    num_events = int(parameters.get("num_events"))
    output_file_location = parameters.get("output_file_location")
    output_file_prefix = parameters.get("output_file_prefix")
    photon_energy_min = float(parameters.get("photon_energy_min")) # GeV
    photon_energy_max = float(parameters.get("photon_energy_max")) # GeV
    Mll2_min = float(parameters.get("Mll2_min")) # GeV^2
    Mll2_max = float(parameters.get("Mll2_max")) # GeV^2
    t_min = float(parameters.get("t_min"))   # GeV^2
    t_max = float(parameters.get("t_max"))   # GeV^2
    Q2_max = float(parameters.get("Q2_max")) # GeV^2
    
    if t_min>0 or t_max>0:
        raise ValueError("'t_min' and 't_max' by definition must be less than 0. Please edit the runcard accordingly.")
    if t_min>t_max:
        raise ValueError("'t_min'>'t_max' in the runcard.")
    rho,X0,_ = get_target_info(parameters)
    
     #######################################################################################
    # Event Generator
    #######################################################################################
    # Event Loop
    
    success=0
    for evt in tqdm(range(num_events)):
        dsigma[0]=0
        weight[0]=0
        psf[0]=0
        flux[0]=0
        flux_brem[0]=0
        flux_epa[0]=0
        acc_eOut[0]=0
        acc_ePlus[0]=0
        acc_eMinus[0]=0
        acc_hOut[0]=0
        gammaE_[0]=0
        t_[0]=0
        smear_t_[0]=0
        W_[0]=0
        smear_W_[0]=0
        Mll2_[0]=0
        Mll2_smear_[0]=0
        smear_eOut.SetPxPyPzE(0,0,0,0)
        smear_ePlus.SetPxPyPzE(0,0,0,0)
        smear_eMinus.SetPxPyPzE(0,0,0,0)
        smear_hOut.SetPxPyPzE(0,0,0,0)
        
        # Set Initial Vectors
        eIn.SetPxPyPzE(0,0,np.sqrt(beam_energy**2-mE**2),beam_energy)
        hIn.SetXYZM(0,0,0,target_mass)
        
        # Do Bremmstrahlung
        flux_brem[0], gammaE = Bremmstrahlung(photon_energy_min,photon_energy_max,beam_energy,X0,target_length)
        flux_epa[0] = N_EquivalentPhotonApproximation(gammaE,beam_energy,Q2_max)
        flux[0] = flux_brem[0]+flux_epa[0]
        
        #flux[0] += N_EquivalentPhotonApproximation(gammaE,beamE,q2max)
        q.SetPxPyPzE(0,0,gammaE,gammaE)
        
        # Get scattered electron
        eOut = eIn-q
        
        # Generate t from allowable range
        t = random.uniform(t_min,t_max)
        if t>0: # Catch in case t>0, make it negative
            t=-t
        
        # Generate Mll2 from allowable range
        Mll2_max0 = get_max_Mll2(gammaE,t,target_mass)
        if(Mll2_max<Mll2_max0):
            Mll2_max0=Mll2_max
        Mll2 = random.uniform(Mll2_min,Mll2_max0)
        
        
        # Set Lab Frame TLorentzVectors
        phi=random.uniform(0,2*np.pi)
        pHadron = get_p_d_lab(t,target_mass)
        cthHadron = get_cth_d_lab(Mll2,gammaE,t,target_mass)
        pDilepton = get_p_dilepton_lab(pHadron,cthHadron,gammaE)
        cthDilepton = get_cth_dilepton_lab(pHadron,cthHadron,gammaE)

        hOut.SetPxPyPzE(pHadron*np.sqrt(1-cthHadron**2)*np.cos(phi),
                        pHadron*np.sqrt(1-cthHadron**2)*np.sin(phi),
                        pHadron*cthHadron,
                        np.sqrt(pHadron**2+target_mass**2))

        dilepton.SetPxPyPzE(-pDilepton*np.sqrt(1-cthDilepton**2)*np.cos(phi),
                            -pDilepton*np.sqrt(1-cthDilepton**2)*np.sin(phi),
                            pDilepton*cthDilepton,
                            np.sqrt(pDilepton**2+Mll2))
        
        # Decay the Dilepton
        GenPhase.SetDecay(dilepton,2,np.array([mE,mE]))
        GenPhase.Generate()
        ep,em = GenPhase.GetDecay(0), GenPhase.GetDecay(1)
        ePlus.SetPxPyPzE(ep.Px(),ep.Py(),ep.Pz(),ep.E())
        eMinus.SetPxPyPzE(em.Px(),em.Py(),em.Pz(),em.E())
    
        # Get theta and phi cm of ePlus
    
        comBOOST=(dilepton).BoostVector()
        theta_cm,phi_cm = get_theta_phi_cm(comBOOST,ePlus,hOut,hIn)
        
        # Get event differential cross section
        dsigma[0] = get_dsigma_dt_dMll2_dcosth_dphi(t,gammaE,Mll2,theta_cm,phi_cm,target_type,target_mass)

        # Get phase space factor
        psf[0] = (t_max-t_min)*(Mll2_max0-Mll2_min)*(4*np.pi)*(photon_energy_max-photon_energy_min)
        
        # Calculate full event weight
        weight[0] = dsigma[0]*psf[0]*flux[0]
        
        # Get acceptances of final state particles
        acc_ePlus[0] = acc_e(eOut)
        acc_ePlus[0] = acc_e(ePlus)
        acc_eMinus[0]= acc_e(eMinus)
        acc_hOut[0]  = acc_h(hOut)
        
        # Store event variables
        W_[0]        = (hIn+q).M()
        gammaE_[0]   = gammaE
        t_[0]        = t
        Mll2_[0]     = Mll2

        # Smear final state particles if allowed based on acceptance
        if(acc_eOut[0]>0):
            temp = smear_particle(eOut,"e-",detector)
            smear_eOut.SetPx(temp.Px())
            smear_eOut.SetPy(temp.Py())
            smear_eOut.SetPz(temp.Pz())
            smear_eOut.SetE(temp.E())
        if(acc_ePlus[0]>0):
            temp = smear_particle(ePlus,"e+",detector)
            smear_ePlus.SetPx(temp.Px())
            smear_ePlus.SetPy(temp.Py())
            smear_ePlus.SetPz(temp.Pz())
            smear_ePlus.SetE(temp.E())
        if(acc_eMinus[0]>0):
            temp = smear_particle(eMinus,"e-",detector)
            smear_eMinus.SetPx(temp.Px())
            smear_eMinus.SetPy(temp.Py())
            smear_eMinus.SetPz(temp.Pz())
            smear_eMinus.SetE(temp.E())
        if(acc_hOut[0]>0):
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
        # Store smeared Mll2 if we reco'd the e+ e-
        if(acc_ePlus[0]>0 and acc_eMinus[0]>0):
            Mll2_smear_[0]=2*smear_ePlus.Pt()*smear_eMinus.Pt()*(np.cosh(smear_ePlus.Eta()-smear_eMinus.Eta())-np.cos(smear_ePlus.Phi()-smear_eMinus.Phi()))
        else:
            Mll2_smear_[0]=-1.0
            
        smear_W_[0]=(hIn+smear_q).M()
        # Fill the tree only when weight > 0
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
    parser = argparse.ArgumentParser(description="Bethe-Heitler Generator Script")
    parser.add_argument("--runcard", type=str, required=True, help="Path to the runcard.card file.")
    parser.add_argument("--batch", type=int, default=0, required=False, help="Batch number (used for parallel running).")
    args = parser.parse_args()
    main(args)