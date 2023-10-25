from ROOT import TFile, TLorentzVector
import numpy as np

# Set Detector System
# ---------------------------------------------
def set_detector(parameters):
    
    detector = parameters.get("detector")
    if detector not in ["SoLID","CLAS12"]:
        raise ValueError(f"'detector={detector}' not recognized.")
        
    global facc, acc, fres_electron, fres_positron, fres_proton, res_electron_p, res_electron_th, res_electron_phi, res_positron_p, res_positron_th, res_positron_phi, res_proton_p, res_proton_th, res_proton_phi

    if(detector=='SoLID'):
        facc = TFile("tables/acceptance_solid_JPsi_electron_target315_output.root","r")
        acc = facc.Get("acceptance_ThetaP_overall")
        fres_electron= TFile("tables/JPsi_electron_resolution_2d.root","READ")
        fres_positron = TFile("tables/JPsi_electron_resolution_2d.root","READ")
        fres_proton = TFile("tables/JPsi_proton_resolution_2d.root","READ")
        res_electron_p = fres_electron.Get("p_resolution")
        res_positron_p = fres_positron.Get("p_resolution")
        res_proton_p = fres_proton.Get("p_resolution")
        res_electron_th = fres_electron.Get("theta_resolution")
        res_positron_th = fres_positron.Get("theta_resolution")
        res_proton_th = fres_proton.Get("theta_resolution")
        res_electron_phi = fres_electron.Get("phi_resolution")
        res_positron_phi = fres_positron.Get("phi_resolution")
        res_proton_phi = fres_proton.Get("phi_resolution")
        res_electron_p.Scale(1.5)
        res_positron_p.Scale(1.5)
        res_proton_p.Scale(1.5)
        res_electron_th.Scale(1.5)
        res_positron_th.Scale(1.5)
        res_proton_th.Scale(1.5)
        res_electron_phi.Scale(1.5)
        res_positron_phi.Scale(1.5)
        res_proton_phi.Scale(1.5)
    elif(detector=='CLAS12'):
        facc = TFile("tables/clasev_acceptance.root","r")
        acc = facc.Get("acceptance_PTheta_ele")
        



def acc_e(P):
    binx = acc.GetXaxis().FindBin(P.Theta()*180.0/np.pi)
    biny = acc.GetYaxis().FindBin(P.P())
    return acc.GetBinContent(binx,biny)

def acc_h(P):
    binx = acc.GetXaxis().FindBin(P.Theta()*180.0/np.pi)
    biny = acc.GetYaxis().FindBin(P.P())
    return acc.GetBinContent(binx,biny)

def smear_particle(P,particleType,det):
    p = P.P()
    th = P.Theta()
    phi= P.Phi()
    if(det=="SoLID"):
        if(particleType=="e+" or particleType=="e-"):
            dp = res_positron_p.GetBinContent(res_positron_p.FindBin(p,th*180.0/np.pi))/100.0
            dth = res_positron_th.GetBinContent(res_positron_th.FindBin(p,th*180.0/np.pi))/100.0
            dphi = res_positron_phi.GetBinContent(res_positron_phi.FindBin(p,th*180.0/np.pi))/100.0
        elif(particleType=="p" or particleType=="d"):
            dp = res_proton_p.GetBinContent(res_proton_p.FindBin(p,th*180.0/np.pi))/100.0
            dth = res_proton_th.GetBinContent(res_proton_th.FindBin(p,th*180.0/np.pi))/100.0
            dphi = res_proton_phi.GetBinContent(res_proton_phi.FindBin(p,th*180.0/np.pi))/100.0
    elif(det=="CLAS12"):
        dp = 0.01
        dth = 0.001
        dphi = 0.004
    pprime = p*np.random.normal(1,dp)
    thprime = th*np.random.normal(1,dth)
    phiprime = phi*np.random.normal(1,dphi)
    smear_P = TLorentzVector(pprime*np.sin(thprime)*np.cos(phiprime),
                            pprime*np.sin(thprime)*np.sin(phiprime),
                            pprime*np.cos(thprime),
                            np.sqrt(pprime**2+P.M2()))
    return smear_P