import numpy as np
from constants import *
import random
from ROOT import TLorentzVector


def get_EC_ED_cth(EA,EB,MA,MB,MC,MD,t):
    # Assuming t > 0
    EC=(EA**2+2*EA*EB+EB**2+MC**2-MD**2)/(2*(EA+EB))
    ED=(EA**2+2*EA*EB+EB**2-MC**2+MD**2)/(2*(EA+EB))
    cth=(-t-MA**2-MD**2+2*EA*ED)/(2*np.sqrt(EA**2-MA**2)*np.sqrt(ED**2-MD**2))
    return EC,ED,cth

def get_tauD(t,target_mass):
    return -t/(4*target_mass**2)

def get_s(Egamma,target_mass):
    return target_mass**2+2*target_mass*Egamma

def get_cth_d_lab(Mll2,Egamma,t,target_mass):
    tauD = get_tauD(t,target_mass)
    s = get_s(Egamma,target_mass)
    return (Mll2+2*(s+target_mass**2)*tauD)/(2*(s-target_mass**2)*np.sqrt(tauD*(1+tauD)))

def get_p_d_lab(t,target_mass):
    return np.sqrt(((2*target_mass**2-t)/(2*target_mass))**2-target_mass**2)

def get_p_dilepton_lab(p_d_lab,cth_d_lab,Egamma):
    return (p_d_lab**2*np.sqrt(1-cth_d_lab**2)*np.abs(np.sqrt(1-cth_d_lab**2))+(Egamma-p_d_lab*cth_d_lab)**2)/np.sqrt(Egamma**2-2*Egamma*p_d_lab*cth_d_lab+p_d_lab**2)

def get_cth_dilepton_lab(p_d_lab,cth_d_lab,Egamma):
    return (Egamma-p_d_lab*cth_d_lab)/np.sqrt(Egamma**2-2*Egamma*p_d_lab*cth_d_lab+p_d_lab**2)

def get_beta(Mll2,target_mass):
    return np.sqrt(1-4*mE**2/Mll2)

def get_max_Mll2(Egamma,t,target_mass):
    return 2*target_mass*Egamma+t-(Egamma/target_mass)*(2*target_mass**2-t-np.sqrt(t**2-4*target_mass**2*t))

def get_delta_T(tauD,t,target_mass):
    # Approximation for small -t , large Mll2, large s
    return np.sqrt(-t*(1-tauD)-target_mass**2*tauD**2)

def get_delta_T2(s,Mll2,t,r,target_mass):
    r0 = np.sqrt((s-target_mass**2)**2)
    cth_cm = (2 * s * (t - 2 * target_mass * target_mass) + (s + target_mass * target_mass)*(s + target_mass * target_mass - Mll2)) / (r0 * r)
    sth_cm = np.sqrt(1-cth_cm**2)
    return sth_cm*r/(2*np.sqrt(s))

def get_theta_phi_cm(comBOOST,ePlus,hOut,hIn):
    # Boost into l-l+ c.m. frame
    ePlus.Boost(-comBOOST)
    hOut.Boost(-comBOOST)
    hIn.Boost(-comBOOST)
    
    # Rotate coordinate system such that dOut lies along -zaxis
    hOutVect=hOut.Vect()*(1/hOut.Vect().Mag())
    rotationVect=1/np.sqrt(2)*(hOutVect-Zhat)
    
    ePlus.Rotate(np.pi,rotationVect)
    hIn.Rotate(np.pi,rotationVect)
    
    # Rotate coordinate system such that dIn has 0 Phi
    hIn_phi = hIn.Phi()
    
    ePlus.Rotate(-hIn_phi,Zhat)
    
    # We are now in the l-l+ c.m. frame of (https://arxiv.org/abs/hep-ph/0110062) Fig5
    theta_cm = ePlus.Angle(Zhat)
    phi_cm   = ePlus.Phi()
    
    ePlus.Rotate(hIn_phi,Zhat)
    ePlus.Rotate(-np.pi,rotationVect)
    hIn.Rotate(-np.pi,rotationVect)
    
    ePlus.Boost(comBOOST)
    hOut.Boost(comBOOST)
    hIn.Boost(comBOOST)
    return theta_cm,phi_cm


def CE1(t,s,Mll2,target_mass):
    return t*(s-target_mass**2)*(s-target_mass**2-Mll2+t)*(Mll2**2+6*Mll2*t+t**2+4*mE**2*Mll2)+(Mll2-t)**2*(t**2*Mll2+target_mass**2*(Mll2+t)**2+4*mE**2*target_mass**2*Mll2)

def CE2(t,s,Mll2,target_mass):
    return -t*(s-target_mass**2)*(s-target_mass**2-Mll2+t)*(Mll2**2+t**2+4*mE**2*(Mll2+2*t-2*mE**2))+(Mll2-t)**2*(-target_mass**2*(Mll2**2+t**2)+2*mE**2*(-t**2-2*target_mass**2*Mll2+4*mE**2*target_mass**2))

def CM1(CE1,tauD,Mll2,t,target_mass):
    return CE1-2*target_mass**2*(1+tauD)*(Mll2-t)**2*(Mll2**2+t**2+4*mE**2*Mll2)

def CM2(CE2,tauD,Mll2,t,target_mass):
    return CE2+2*target_mass**2*(1+tauD)*(Mll2-t)**2*(Mll2**2+t**2+4*mE**2*(Mll2-t-2*mE**2))


conv = 0.197327
a1,a2,a3,a4,als1,als2,als3,als4 = 1.57057 * conv**2, 12.23792 * conv**2, -42.04576 * conv**2, 27.92014 * conv**2, 1.52501 * conv**2, 8.75139 * conv**2,  15.97777 * conv**2, 23.20415* conv**2
b1,b2,b3,b4,bes1,bes2,bes3,bes4 = 0.07043 * conv, 0.14443 * conv, -0.27343 * conv, 0.05856 * conv, 43.67795 * conv**2, 30.05435 * conv**2, 16.43075 * conv**2,  2.80716* conv**2
c1,c2,c3,c4,gas1,gas2,gas3,gas4 = -0.16577, 0.27557 ,  -0.05382 , -0.05598 ,  1.87055 * conv**2, 14.95683 * conv**2,  28.04312 * conv**2, 41.12940* conv**2


def get_g0_g1_g2(Q2):
    g0 = a1/(als1+Q2)+a2/(als2+Q2)+a3/(als3+Q2)+a4/(als4+Q2)
    g1 = np.sqrt(Q2)*(b1/(bes1+Q2)+b2/(bes2+Q2)+b3/(bes3+Q2)+b4/(bes4+Q2))
    g2 = Q2 * (c1/(gas1+Q2)+c2/(gas2+Q2)+c3/(gas3+Q2)+c4/(gas4+Q2))
    return g0,g1,g2

def get_GC_GQ_GM(t,target_mass):
    #t_fm = t/0.0389 # GeV^2 --> fm^-1
    tau = get_tauD(t,target_mass)
    g0,g1,g2 = get_g0_g1_g2(-t)
    GC = 1/(1-t/4/0.807)**4/(2*tau+1)*((1-2/3*tau)*g0+8/3*np.sqrt(2*tau)*g1+2/3*(2*tau-1)*g2)
    GQ = 1/(1-t/4/0.807)**4/(2*tau+1)*(-g0+np.sqrt(2/tau)*g1-(tau+1)/tau*g2)
    GM = 1/(1-t/4/0.807)**4/(2*tau+1)*(2*g0+2*(2*tau-1)/np.sqrt(2*tau)*g1-2*g2)
    return GC,GQ,GM



def get_dsigma_dt_dMll2(t,Egamma,Mll2,target_mass):
    B = get_beta(Mll2,target_mass)
    s = get_s(Egamma,target_mass)
    tauD = get_tauD(t,target_mass)
    ce1 = CE1(t,s,Mll2,target_mass)
    ce2 = CE2(t,s,Mll2,target_mass)
    cm1 = CM1(ce1,tauD,Mll2,t,target_mass)
    cm2 = CM2(ce2,tauD,Mll2,t,target_mass)
    CE = ce1 + (ce2/B)*np.log((1+B)/(1-B))
    CM = cm1 + (cm2/B)*np.log((1+B)/(1-B))
    GC,GQ,GM = get_GC_GQ_GM(t,target_mass)
    factor1 = (4*alpha_em**3*B)/((s-target_mass**2)**2*t**2*(Mll2-t)**4)
    factor2 = (CE*(GC**2+8/9*tauD**2*GQ**2)+CM*2/3*tauD*GM**2)
    return factor1*factor2*(0.389379*10**6) # nb/GeV^4

def get_dsigma_dt_dMll2_dcosth_dphi(t,Egamma,Mll2,theta,phi,target_type,target_mass):
    beta = get_beta(Mll2,target_mass)
    s = get_s(Egamma,target_mass)
    tauD = get_tauD(t,target_mass)
    r = np.sqrt((s-target_mass**2-Mll2)**2-4*target_mass**2*Mll2)
    deltaT = get_delta_T2(s,Mll2,t,r,target_mass)
    a = beta * r * np.cos(theta)
    b = beta * ((Mll2*(s-target_mass**2-Mll2)+t*(s-target_mass**2+Mll2))/r)*np.cos(theta)-beta*((2*(s-target_mass**2)*np.sqrt(Mll2)*deltaT)/r)*np.sin(theta)*np.cos(phi)
    L = ((Mll2-t)**2-b**2)/4
    A = (s-target_mass**2)**2*deltaT**2-t*a*(a+b)-target_mass**2*b**2-t*(4*target_mass**2-t)*Mll2+mE**2/L*(((Mll2-t)*(a+b)-(s-target_mass**2)*b)**2+t*(4*target_mass**2-t)*(Mll2-t)**2)
    B = (Mll2+t)**2 + b**2 + 8 * mE**2 * Mll2 - 4*mE**2*(t+2*mE**2)/L*(Mll2-t)**2
    if(target_type=="p"): 
        F1p = (1. / ((1 - t / 0.71)*(1 - t / 0.71)))*(1 / (1 - t / (4. * target_mass * target_mass)))*(1 - 2.79 * t / (4 * target_mass * target_mass))
        F2p = (1. / ((1 - t / 0.71)*(1 - t / 0.71)))*(1 / (1 - t / (4. * target_mass * target_mass)))*(2.793 - 1)
        factor1 = (alpha_em**3*beta)/(4*np.pi*(s-target_mass**2)**2*(-t)*L)
        factor2 = (F1p**2-t/4/target_mass**2*F2p**2)*A/(-t)+(F1p+F2p)**2*B/2
        return factor1*factor2*(0.389379*10**6) # nb/GeV^4
    elif(target_type=="d"):
        GC,GQ,GM = get_GC_GQ_GM(t,target_mass)
        CE = A
        CM = A + (4*target_mass**2-t)*B/2
        factor1 = (alpha_em**3*beta)/(4*np.pi*(s-target_mass**2)**2*t**2*L)
        factor2 = (CE*(GC**2+8/9*GQ**2*tauD**2)+CM*2/3*tauD*GM**2)
        return factor1*factor2*(0.389379*10**6) # nb/GeV^4