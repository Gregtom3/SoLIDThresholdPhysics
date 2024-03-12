import numpy as np
from constants import *
from ROOT import TF1
import random 

# bremm = TF1("bremm","(4/3-4*x/3+x*x)/x",0.01,1)
# bremm.SetNpx(1000)

# Bremmstrahlung Option A:
# 1. Generate the "y" using the PDF of dsigma_dy
# 2. Return weight of full integral over [ymin,ymax]
# def N_Bremmstrahlung(kmin,kmax,beamE,d,X0):
#     ymin = kmin/beamE
#     ymax = kmax/beamE
#     y = bremm.GetRandom(ymin,ymax)
#     # We use "0.5*d" to simulate half the target as the radiator
#     flux = 0.5* d/X0 * (4.0 / 3.0 * np.log(ymax / ymin) - 4.0 / 3.0 * (ymax - ymin) + 1.0 / 2.0 * (ymax * ymax - ymin * ymin))
#     return  2.0 * np.pi * flux, y*beamE


# Bremmstrahlung Option B:
# 1. Uniformly generate the "y" for dsigma_dy
# 2. Return "flux" which is a weight given by |ds/dy| * (ymax-ymin)
def Bremmstrahlung(kmin,kmax,beamE,target_length,X0):
    ymin = kmin/beamE
    ymax = kmax/beamE
    y = random.uniform(ymin,ymax)
    # We use "0.5*d" to simulate half the target as the radiator
    flux = 0.5*target_length/X0/(y*beamE)*(4/3 - 4/3 * y + y**2)
    return  flux, y*beamE

def N_EquivalentPhotonApproximation(gammaE,beamE,q2max):
    # (https://lss.fnal.gov/archive/other/lpc-94-35.pdf)
    # Factors of mE^2 are needed atop the fractions 1/q2min and 1/q2max for units
    # I believe this is a mistake in the derivation
    y = gammaE/beamE
    q2min = mE**2 * y**2 / (1-y)
    return (alpha_em*y)/(4*np.pi)*((2*q2min+4*mE**2)/(q2min)*np.log(q2max/q2min)-4*mE**2/q2min+4*mE**2/q2max)/beamE



