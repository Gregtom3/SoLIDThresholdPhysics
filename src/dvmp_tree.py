from ROOT import TLorentzVector, TFile, TTree
import numpy as np

# 4-momentum of Particles
# ---------------------------------------------
eIn, hIn, eOut, hOut, q, VM, ePlus, eMinus = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()

# 4-momentum of Smeared Particles
# ---------------------------------------------
smear_eOut, smear_hOut, smear_ePlus, smear_eMinus, smear_q, smear_VM = TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector(), TLorentzVector()

# Returned eventgen variables
# ---------------------------------------------
weight=np.array([0.0])
dsigma=np.array([0.0])
decay_weight=np.array([0.0])
psf=np.array([0.0])
flux=np.array([0.0])
flux_brem=np.array([0.0])
flux_epa=np.array([0.0])
acc_eOut=np.array([0.0])
acc_hOut=np.array([0.0])
acc_ePlus=np.array([0.0])
acc_eMinus=np.array([0.0])
Q2_=np.array([0.0])
smear_Q2_=np.array([0.0])
gammaE_=np.array([0.0])
smear_gammaE_=np.array([0.0])
t_=np.array([0.0])
smear_t_=np.array([0.0])
W_=np.array([0.0])
smear_W_=np.array([0.0])
m_vm_=np.array([0.0])
smear_m_vm_=np.array([0.0])
jacobian = np.array([0.0])
smear_jacobian = np.array([0.0])

def create_dvmp_ttree(parameters):
    
    target_type = parameters.get("target_type")            # p/d
    num_events = int(parameters.get("num_events"))
    beam_energy = float(parameters.get("beam_energy"))     # GeV
    model_type = parameters.get("model_type")
    output_file_location = parameters.get("output_file_location")
    output_file_prefix = parameters.get("output_file_prefix")
    process = parameters.get("process")
    batch = parameters.get("batch")
    fname = f"{output_file_location}/{output_file_prefix}_e{target_type}_beamE_{beam_energy:.2f}_evts_{num_events}_{process}_{model_type}_{batch}.root"
    
    outFile=TFile(fname,"RECREATE")
    # Generate TTree
    # ------------------------------
    outTree=TTree("tree","tree")
    # Create the branches for the TTree
    # ------------------------------
    outTree.Branch("eIn",eIn)
    outTree.Branch("hIn",hIn)
    outTree.Branch("eOut",eOut)
    outTree.Branch("hOut",hOut)
    outTree.Branch("ePlus",ePlus)
    outTree.Branch("eMinus",eMinus)
    outTree.Branch("smear_eOut",smear_eOut)
    outTree.Branch("smear_hOut",smear_hOut)
    outTree.Branch("smear_ePlus",smear_ePlus)
    outTree.Branch("smear_eMinus",smear_eMinus)
    outTree.Branch("q",q)
    outTree.Branch("smear_q",smear_q)
    outTree.Branch("smear_VM",smear_VM)
    outTree.Branch("dsigma",dsigma,"dsigma/D")
    outTree.Branch("weight",weight,"weight/D")
    outTree.Branch("decay_weight",decay_weight,"decay_weight/D")
    outTree.Branch("psf",psf,"psf/D")
    outTree.Branch("flux",flux,"flux/D")
    outTree.Branch("flux_brem",flux_brem,"flux_brem/D")
    outTree.Branch("flux_epa",flux_epa,"flux_epa/D")
    outTree.Branch("acc_eOut",acc_eOut,"acc_eOut/D")
    outTree.Branch("acc_hOut",acc_hOut,"acc_hOut/D")
    outTree.Branch("acc_ePlus",acc_ePlus,"acc_ePlus/D")
    outTree.Branch("acc_eMinus",acc_eMinus,"acc_eMinus/D")
    outTree.Branch("Q2",Q2_,"Q2/D")
    outTree.Branch("smear_Q2",smear_Q2_,"smear_Q2/D")
    outTree.Branch("gammaE",gammaE_,"gammaE/D")
    outTree.Branch("smear_gammaE",smear_gammaE_,"smear_gammaE/D")
    outTree.Branch("t",t_,"t/D")
    outTree.Branch("smear_t",smear_t_,"smear_t/D")
    outTree.Branch("W",W_,"W/D")
    outTree.Branch("smear_W",smear_W_,"smear_W/D")
    outTree.Branch("m_vm",m_vm_,"m_vm/D")
    outTree.Branch("smear_m_vm",smear_m_vm_,"smear_m_vm/D")
    outTree.Branch("J",jacobian,"J/D")
    outTree.Branch("smear_J",smear_jacobian,"smear_J/D")
    
    outTreeAcc=outTree.CloneTree()
    outTreeAcc.SetName("treeAcc")
    
    return outFile, outTree, outTreeAcc

