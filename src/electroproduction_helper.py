import random
import numpy as np
from ROOT import TLorentzVector
from constants import *

def VirtualPhoton(kmin,kmax,beamE,target_mass,eIn,eOut,hIn):
    cthmax,cthmin=1,0.995
    Pe  = random.uniform(beamE-kmax,beamE-kmin)
    cth = random.uniform(cthmin,cthmax)
    sth = np.sqrt(1-cth**2)
    phi = random.uniform(0,2*np.pi)
    eOut.SetXYZM(Pe * sth * np.cos(phi),Pe * sth * np.sin(phi),Pe * cth, mE)
    qq=eIn-eOut
    W2 = (qq+hIn)*(qq+hIn)
    
    if(W2 < target_mass**2):
        return 0, TLorentzVector()
    Q2 = -qq*qq
    flux = np.sqrt((hIn*qq)**2 + Q2 * target_mass**2)/np.sqrt((eIn*hIn)**2 - mE**2 * target_mass**2)
    amp = (2.0 * Q2 - 4.0 * mE**2) / (Q2**2)
    phase = eOut.P()**2/(2.0*eOut.E()*(2.0*np.pi)**3)
    volume = 2.0 * np.pi * (cthmax-cthmin)
    y = (hIn * qq)/(hIn*eIn)
    gy = hIn.M() * np.sqrt(Q2)/(hIn * eIn)
    eps = (1.0 - y - 0.25 * gy * gy) / (1.0 - y + 0.5 * y * y + 0.25 * gy * gy)
    couple = 4.0 * np.pi * alpha_em
    return couple * flux * amp * phase * volume / (1.0 - eps),qq