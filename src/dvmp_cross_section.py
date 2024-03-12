from scipy.interpolate import LinearNDInterpolator
import numpy as np
import pandas as pd
from constants import *
from scipy import integrate
from scipy.interpolate import make_interp_spline


class dvmpProduction:

    def __init__(self, targetType, modelType):
        self.targetType = targetType
        self.modelType = modelType
        self.m_VM = 3.096
        if self.targetType == "p":
            self.mT = 0.938272
        elif self.targetType == "d":
            self.mT = 1.8761358
        else:
            raise ValueError(f"Unknown targetType '{targetType}'.")
        
        # Determine the differential cross section & total cross section formulas
        if self.targetType == "p" and self.modelType == "23g":
            self.dsigma = self.dsigma_dt___proton_23g
            self.sigma  = self.sigma___proton_23g
            self.N2g = 6.499e3
            self.N3g = 2.894e3
            self.v = 1/(16*np.pi)
            
        elif self.targetType == "d" and self.modelType == "PomeronLQCD":
            self.dsigma = self.dsigma_dt___deuteron_PomeronLQCD
            self.sigma  = self.sigma___deuteron_PomeronLQCD
            self.df_deuteron    = pd.read_csv("/work/clas12/users/gmat/SoLIDThresholdPhysics/tables/dvmp_xsec_deuteron.csv",sep="\t").sort_values(by=['Beam Energy [GeV]', '-t [GeV**2]'])
            # Create Interpolator
            self.E_deuteron = self.df_deuteron['Beam Energy [GeV]'].values
            self.T_deuteron = self.df_deuteron['-t [GeV**2]'].values
            self.DSDT_deuteron = self.df_deuteron.iloc[:, 3].values
            self.interp_deuteron = LinearNDInterpolator(list(zip(self.E_deuteron, self.T_deuteron)), self.DSDT_deuteron)

            # Compute tmin and tmax for each energy
            self.unique_energies = np.unique(self.E_deuteron)
            self.tmtm_deuteron = []

            for energy in self.unique_energies:
                idx_of_tmin = np.where(self.E_deuteron==energy)[0][0]
                idx_of_tmax = np.where(self.E_deuteron==energy)[0][-1]
                tmin=self.T_deuteron[idx_of_tmin]
                tmax=self.T_deuteron[idx_of_tmax]
                self.tmtm_deuteron.append([energy, tmin, tmax, idx_of_tmin, idx_of_tmax])

            self.tmtm_deuteron = np.array(self.tmtm_deuteron)
    
        else:
            raise ValueError(f"For class 'dvmpProduction', targetType='{targetType}' and modelType='{modelType}' is not a supported combination.")
    
    
    # Get smallest kinematically possible tmin
    def get_tmin(self,Eg):
        mT = self.mT
        W = np.sqrt(2*mT*Eg+mT**2)
        m1_2 = 0;
        m2_2 = np.power(mT,2);
        m3_2 = np.power(self.m_VM,2);
        m4_2 = np.power(mT,2);
        p1 = np.sqrt(np.power(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
        p2 = np.sqrt(np.power(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);
        return m2_2 + m4_2 - 2 * np.sqrt(m2_2 + p1*p1) * np.sqrt(m4_2 + p2*p2) + 2 * p1 * p2;

    # Get largest kinematically possible tmax
    def get_tmax(self,Eg):
        mT = self.mT
        W = np.sqrt(2*mT*Eg+mT**2)
        m1_2 = 0;
        m2_2 = np.power(mT,2);
        m3_2 = np.power(self.m_VM,2);
        m4_2 = np.power(mT,2);
        p1 = np.sqrt(np.power(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
        p2 = np.sqrt(np.power(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);
        return m2_2 + m4_2 - 2 * np.sqrt(m2_2 + p1*p1) * np.sqrt(m4_2 + p2*p2) - 2 * p1 * p2;
    
    
    # Default function for differential cross sections
    def dsigma(self,pars):
        pass
    
    # Default function for total cross sections
    def sigma(self,pars):
        pass
    
    # ------------------------------------------------
    # Deuteron PomeronLQCD model
    #-------------------------------------------------
    # Differential cross section
    
    def dsigma_dt___deuteron_PomeronLQCD(self,pars):
        
        # Recover parameters
        photon_energy = pars[0]
        t = pars[1]
        
        # Determine index of closest energy in table
        idx_of_E = np.abs(photon_energy - self.tmtm_deuteron[:,0]).argmin()
        
        # Get table energy, tmin, tmax
        closest_photon_energy, tmin, tmax = self.tmtm_deuteron[idx_of_E,0],self.tmtm_deuteron[idx_of_E,1],self.tmtm_deuteron[idx_of_E,2]
        
        # Check t range
        if -t<tmin or -t>tmax:
            return -1
        
        # Calculate dsigma_dt
        dsigma_dt = self.interp_deuteron(closest_photon_energy,-t)
        return dsigma_dt
    
    # Total cross section
    def sigma___deuteron_PomeronLQCD(self,pars):
        
       
        gammaE = pars[0] # photon_energy as per convention
        
        # Early return for gammaE <= 6
        if gammaE <= 6:
            idx = int(self.tmtm_deuteron[self.tmtm_deuteron[:, 0] >= gammaE].min())
            DSDT_arr = self.DSDT_deuteron[int(self.tmtm_deuteron[idx, 3]):int(self.tmtm_deuteron[idx, 4])]
            T_arr = self.T_deuteron[int(self.tmtm_deuteron[idx, 3]):int(self.tmtm_deuteron[idx, 4])]
            return integrate.simpson(DSDT_arr, T_arr)
        
        # Find the nearest energy bounds from the table
        gammaE_0 = self.tmtm_deuteron[self.tmtm_deuteron[:, 0] < gammaE, 0].max()
        gammaE_1 = self.tmtm_deuteron[self.tmtm_deuteron[:, 0] > gammaE, 0].min()
        
        # Get the indices of the energies
        idx_0 = list(self.tmtm_deuteron[:, 0]).index(gammaE_0)
        idx_1 = list(self.tmtm_deuteron[:, 0]).index(gammaE_1)
        
        # Extract relevant dsigma/dt and t arrays for both energy bounds
        DSDT_arr_0 = self.DSDT_deuteron[int(self.tmtm_deuteron[idx_0, 3]):int(self.tmtm_deuteron[idx_0, 4])]
        T_arr_0 = self.T_deuteron[int(self.tmtm_deuteron[idx_0, 3]):int(self.tmtm_deuteron[idx_0, 4])]
        DSDT_arr_1 = self.DSDT_deuteron[int(self.tmtm_deuteron[idx_1, 3]):int(self.tmtm_deuteron[idx_1, 4])]
        T_arr_1 = self.T_deuteron[int(self.tmtm_deuteron[idx_1, 3]):int(self.tmtm_deuteron[idx_1, 4])]
        
        # Calculate the weighting based on how close gammaE is to gammaE_0 and gammaE_1
        diff = 1 - (gammaE - gammaE_0) / (gammaE_1 - gammaE_0)
        
        # Linearly interpolate between the two sigma values from gammaE_0 and gammaE_1
        sigma_val = diff * integrate.simpson(DSDT_arr_0, T_arr_0) + (1 - diff) * integrate.simpson(DSDT_arr_1, T_arr_1)
        
        return sigma_val
        
    def dsigma_dt___proton_23g(self,pars):
        Eg = pars[0]
        t  = pars[1]
        W = np.sqrt(2*self.mT*Eg+self.mT**2)
        x = (2*self.mT*mJpsi+mJpsi**2)/(W**2-self.mT**2)
        return  (self.N2g * self.v * (1-x)**2*np.exp(1.13*t)/(mJpsi**2)+self.N3g*self.v*np.exp(1.13*t)/(mJpsi**4))
        
    def sigma___proton_23g(self,pars):
        Eg = pars[0]
        tmin = self.get_tmin(Eg)
        tmax = self.get_tmax(Eg)
        W = np.sqrt(2*self.mT*Eg+self.mT**2)
        x = (2*self.mT*mJpsi+mJpsi**2)/(W**2-self.mT**2)
        
        def bound(t):
            return (1/1.13)*np.exp(1.13*t)*self.v*(self.N3g+mJpsi**2*self.N2g*(1-x)**2)/mJpsi**4
        
        return bound(tmax)-bound(tmin)
# N2g = 6.499e3
# N3g = 2.894e3
# v = 1/(16*np.pi)

# # -------------------------------
# # 
# # PROTON MODEL 
# # (Pomeron LQCD)
# # ds --> "d-sigma"
# # # -------------------------------
# df_proton=pd.read_csv("/work/clas12/users/gmat/SoLIDThresholdPhysics/tables/dvmp_xsec_proton_LQCD.csv",sep=",")
# dfV_proton = df_proton.values
# df_W_proton = dfV_proton[:,0]
# df_Th_proton = dfV_proton[:,1]
# df_DS_proton = dfV_proton[:,2]
# interp_proton = LinearNDInterpolator(list(zip(df_W_proton,np.cos(np.pi*df_Th_proton/180))),df_DS_proton)

# # Get dsdcth from interpolation
# def dsdcth_proton(W,cth):
#     return interp_proton(W,cth)

# # -------------------------------
# # 
# # DEUTERON MODEL 
# # Pomeron LQCD Coherent production from Harry
# #
# # -------------------------------

# # Load the Model
# df_deuteron=pd.read_csv("/work/clas12/users/gmat/SoLIDThresholdPhysics/tables/dvmp_xsec_deuteron.csv",sep="\t")
# df_deuteron=df_deuteron.sort_values(by=['Beam Energy [GeV]','-t [GeV**2]'])
# dfV_deuteron=df_deuteron.values
# df_E_deuteron=dfV_deuteron[:,1]
# df_T_deuteron=dfV_deuteron[:,2]
# df_DSDT_deuteron=dfV_deuteron[:,3]
# interp_deuteron=LinearNDInterpolator(list(zip(df_E_deuteron,df_T_deuteron)),df_DSDT_deuteron)

# # Get tmin and tmax for each energy
# # Array of [Energy, tmin, tmax, idx_of_tmin, idx_of_tmax]
# tmtm_deuteron = np.zeros((np.unique(df_E_deuteron).size,5))
# idxtmp=0
# for tmpE in np.unique(df_E_deuteron):
#     idx_of_tmin = np.where(df_E_deuteron==tmpE)[0][0]
#     idx_of_tmax = np.where(df_E_deuteron==tmpE)[0][-1]
#     tminTemp=df_T_deuteron[idx_of_tmin]
#     tmaxTemp=df_T_deuteron[idx_of_tmax]
#     tmtm_deuteron[idxtmp]=(tmpE,tminTemp,tmaxTemp,idx_of_tmin,idx_of_tmax)
#     idxtmp+=1

# #----------------------------------------------------------------------
# def get_tmtm(gammaE,targetType):
#     if targetType=="p":
#         raise ValueError("'get_tmtm()' currently not implemented for proton target.")
#     else:
#         tmtm = tmtm_deuteron
    
#     idx_of_E = np.abs(gammaE - tmtm[:,0]).argmin()
#     return tmtm[idx_of_E,0],tmtm[idx_of_E,1],tmtm[idx_of_E,2]

# # Get dsdt from interpolation
# def dsdt(gammaE,t,targetType):
#     if targetType=="p":
#         raise ValueError("'dsdt()' currently not implemented for proton target.")
#     else:
#         return interp_deuteron(gammaE,t)

# def get_xsec(targetType,mT,modelType,gammaE,W,t,cth,pCM_Initial,pCM_Final):
    
    
#     # J = "Jacobian"
#     # Since we generate gamma+(p,d) --> V + (p',d') with a uniform angular distribution (cth uniform)
#     # we must multiply cross sections dsigma/dt by the appropriate dt/dcth Jacobian
    
#     # PROTON CROSS SECTION
#     # ------------------------------------
#     if(targetType=='p'):
#         if(modelType=='PomeronLQCD'):
#             the_dsdcth = dsdcth_proton(W,cth)
#             return the_dsdcth

#         elif(modelType=='23g'):
#             x = (2*mP*self.m_VM+self.m_VM**2)/(W**2-mP**2)
#             J = 2*pCM_Initial*pCM_Final/(2*np.pi)
#             return J * (N2g * v * (1-x)**2*np.exp(1.13*t)/(self.m_VM**2)+N3g*v*np.exp(1.13*t)/(self.m_VM**4))

        
        
#     # DEUTERON CROSS SECTION
#     # ------------------------------------
#     elif(targetType=='d' and modelType=='PomeronLQCD'):
#         # Also store temporary "closest photon energy" from cross section data
#         Egamma_temp,tmin,tmax = get_tmtm(gammaE,"d")
#         if(-t<tmin or -t>tmax):
#             return -1
#         # Get dsdt of (gamma + p,d --> JPsi + p,d)
#         # We use Egamma_temp as opposed to gammaE from Monte Carlo
#         # because the interpolation for dsdt(Egamma,t) is guaranteed to work
#         # when using an Egamma pulled specifically from the crossSection data
#         the_dsdt=dsdt(Egamma_temp,-t,targetType)

#         J = 2*pCM_Initial*pCM_Final/(2*np.pi)
#         return J * the_dsdt
    
#     else:
#         raise ValueError(f"'get_xsec()' has no implementation for 'targetType={targetType}' & 'modelType={modelType}'")
        
        
# # Return cross section given photon energy
# def xsec(gammaE,targetType):
#     if targetType == "d":
#         tmtm = tmtm_deuteron
#         df_DSDT = df_DSDT_deuteron
#         df_T = df_T_deuteron
#         df_E = df_E_deuteron
#         if(gammaE<=6):
#             idx=int(tmtm[tmtm[:,0] >= gammaE].min())
#             DSDT_arr=df_DSDT[int(tmtm[idx,3]):int(tmtm[idx,4])]
#             T_arr=df_T[int(tmtm[idx,3]):int(tmtm[idx,4])]
#             return integrate.simpson(DSDT_arr,T_arr)
#         gammaE_0=tmtm[tmtm[:,0] < gammaE,0].max()
#         gammaE_1=tmtm[tmtm[:,0] > gammaE,0].min()
#         diff=(gammaE-gammaE_0)/(gammaE_1-gammaE_0)
#         idx_0=list(tmtm[:,0]).index(gammaE_0)
#         idx_1=list(tmtm[:,0]).index(gammaE_1)
#         DSDT_arr_0 = df_DSDT[int(tmtm[idx_0,3]):int(tmtm[idx_0,4])]
#         T_arr_0 = df_T[int(tmtm[idx_0,3]):int(tmtm[idx_0,4])]
#         DSDT_arr_1 = df_DSDT[int(tmtm[idx_1,3]):int(tmtm[idx_1,4])]
#         T_arr_1 = df_T[int(tmtm[idx_1,3]):int(tmtm[idx_1,4])]
#         diff=1-(gammaE-gammaE_0)/(gammaE_1-gammaE_0)
#         return diff*integrate.simpson(DSDT_arr_0,T_arr_0)+(1-diff)*integrate.simpson(DSDT_arr_1,T_arr_1)
#     elif targetType == "p":
#         raise ValueError(f"'xsec()' has no implementation for proton currently.")
        
        
# def get_tmin(Eg,targetType):
#     if targetType=="d":
#         mT = 1.875612928
#     elif targetType=="p":
#         mT = 0.938272
#     W = np.sqrt(2*mT*Eg+mT**2)
#     m1_2 = 0;
#     m2_2 = np.power(mT,2);
#     m3_2 = np.power(3.0969,2);
#     m4_2 = np.power(mT,2);
#     p1 = np.sqrt(np.power(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
#     p2 = np.sqrt(np.power(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);
#     return m2_2 + m4_2 - 2 * np.sqrt(m2_2 + p1*p1) * np.sqrt(m4_2 + p2*p2) + 2 * p1 * p2;

# def get_tmax(Eg,targetType):
#     if targetType=="d":
#         mT = 1.875612928
#     elif targetType=="p":
#         mT = 0.938272
#     W = np.sqrt(2*mT*Eg+mT**2)
#     m1_2 = 0;
#     m2_2 = np.power(mT,2);
#     m3_2 = np.power(3.0969,2);
#     m4_2 = np.power(mT,2);
#     p1 = np.sqrt(np.power(W * W - m1_2 - m2_2,2) - 4.0 * m1_2 * m2_2) / (2.0 * W);
#     p2 = np.sqrt(np.power(W * W - m3_2 - m4_2,2) - 4.0 * m3_2 * m4_2) / (2.0 * W);
#     return m2_2 + m4_2 - 2 * np.sqrt(m2_2 + p1*p1) * np.sqrt(m4_2 + p2*p2) - 2 * p1 * p2;





