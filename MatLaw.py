import numpy as np
from constants import (
    rhoi,
    rhoa,
    LCal,
    Lm,
    mH2O,
    kB,
    D0,
    ka,
    ki,
    Ca,
    Ci, 
    Tfus,
    g, 
    ka0,
    ka1,
    ka2,
    T_ref_L,
    a0,
    a1,
    a2,
    f,
    rho_ref,
    T_ref_C,
    c1,
    c2,
    c3,
    c4,
    c5,
    c6,
    R_v,
    eta_0,
    a_eta,
    b_eta,
    c_eta,
    D_rate_literature,
    pl1,
    pl2
)

############################################
###    PARAMETRIZATIONS OF RHOV_SAT      ###
############################################

###Rhov_sat following Clausius-Clapeyron
def MatLaw_RhovSatClausiusClapeyron(T):
     x = (LCal * mH2O) / (rhoi * kB)
     rho_v = rho_ref * np.exp(x * ((1 / T_ref_C) - (1 / T)))
     return rho_v
    
###Rhov_sat following Libbrecht
def MatLaw_RhovSatLibbrecht(T):
    rho_v = (np.exp(-T_ref_L / T) / (f * T) * (a0 + a1 * (T - 273) + a2 * (T - 273) ** 2))  
    return rho_v   

###Rhov_sat following Hansen
def MatLaw_RhovSatHansen(T):
    rho_v = (10.0 ** (c1 / T + c2 * np.log(T) / np.log(10) + c3 * T + c4 * T ** 2 + c5))* c6/ R_v/ T
    return rho_v 

################################################
###        CORRESPONDING dRHOV_SAT/dT        ###
################################################  

###dRhov_sat/dT following Clausius-Clapeyron
def MatLaw_dRhovSatdTClausiusClapeyron(T):   
     x = (LCal * mH2O) / (rhoi * kB)
     drho_v_dT = x / T ** 2 * rho_ref * np.exp(x * ((1 / T_ref_C) - (1 / T)))
     return drho_v_dT
    
###dRhov_sat/dT following Libbrecht
def MatLaw_dRhovSatdTLibbrecht(T):
    drho_v_dT = np.exp(-T_ref_L / T)/ (f * T ** 2) * ((a0 - a1 * 273 + a2 * 273 ** 2) * (T_ref_L / T - 1) + (a1 - a2 * 2 * 273) * T_ref_L + a2 * T ** 2 * (T_ref_L / T + 1))
    return drho_v_dT   

###dRhov_sat/dT following Hansen
def MatLaw_dRhovSatdTHansen(T):
    rho_v = (10.0 ** (c1 / T + c2 * np.log(T) / np.log(10) + c3 * T + c4 * T ** 2 + c5))* c6/ R_v/ T
    drho_v_dT = rho_v * np.log(10)* (-c1 / T ** 2 + c2 / (T * np.log(10)) + c3 + 2 * c4 * T) - rho_v / T
    return drho_v_dT  

################################################
###        PARAMETRIZATIONS OF Keff          ###
################################################ 
    
def MatLaw_keffCalonne(Phii):
    keff = ka0 + ka1 * (rhoi * Phii) + ka2 * (rhoi * Phii) ** 2
    return keff

def MatLaw_keffHansen(Phii, T):
    dSWDdT = MatLaw_dRhovSatdTLibbrecht(T)
    keff = Phii * ((1 - Phii) * ka + Phii * ki ) + (1 - Phii) * (ki*ka/(Phii*(ka + Lm * D0 * dSWDdT) + (1-Phii) * ki)) 
    return keff

################################################
###        PARAMETRIZATIONS OF Deff          ###
################################################ 
    
def MatLaw_DeffCalonne(Phii):
    x = 2 / 3 - Phii  #### Upper limit on diffusion, Implementation of Simson 2021 
#    x = 2 ### Dummy to always have b=1 for all phii when we don't want to use the simson threshold on diffusion
    b = np.heaviside(x, 1)
    Deff = D0 * (1 - 3 / 2 * Phii) * b
    return Deff

def MatLaw_DeffHansen(Phii, T):
    dSWDdT = MatLaw_dRhovSatdTLibbrecht(T)
    Deff = Phii * (1 - Phii) * D0 + (1 - Phii) * (ki*D0/(Phii*(ka + Lm * D0 * dSWDdT) + (1-Phii) * ki)) 
    return Deff

#################################################
###        PARAMETRIZATIONS OF rhoCeff        ###
################################################# 
    
def MatLaw_rhoCeffNoAir(Phii):
    rhoCeff = Phii * rhoi * Ci 
    return rhoCeff

def MatLaw_rhoCeffWithAir(Phii):
    rhoCeff = Phii * rhoi * Ci + (1 - Phii) * rhoa * Ca
    return rhoCeff

#################################################
###        PARAMETRIZATIONS OF vkin           ###
################################################# 

def MatLaw_vkin(T):
    vkin=np.sqrt((kB*T)/(2*np.pi*mH2O))
    return vkin

######################################################
###        PARAMETRIZATIONS OF VISCOSITY           ###
###################################################### 

def MatLaw_ViscosityVionnet(Phii,T): 
     eta = (eta_0* rhoi * Phii / c_eta) * np.exp(a_eta * (Tfus - T) + b_eta * rhoi * Phii)
     return eta


def MatLaw_ViscosityVionnetReduced(Phii, T): ##For enhancing settlement for scenario3 of Schurholt
    eta = (eta_0 * rhoi * Phii / c_eta) * np.exp(a_eta * (Tfus - T) + b_eta * rhoi * Phii)
    return eta/10000

def MatLaw_ViscosityVionnetTconst(Phii): 
     Tconst = 263
     eta = (eta_0* rhoi * Phii / c_eta) * np.exp(a_eta * (Tfus - Tconst) + b_eta * rhoi * Phii)
     return eta

def MatLaw_ViscosityVionnetPhiiconst(Phii, T): 
     Phiiconst = 0.1125
     restrict = np.exp(pl1 * Phii - pl2) + 1 # power law that becomes very high when Phii reaches 0.95 to render further settlement almost impossible
     eta = (eta_0* rhoi * Phiiconst / c_eta) * np.exp(a_eta * (Tfus - T) + b_eta * rhoi * Phiiconst) * restrict
     return eta

def MatLaw_ViscosityConstn1(Phii):
    Tconst = 263
    Phiiconst = 0.1125
    restrict = np.exp(pl1 * Phii - pl2) + 1 # power law that becomes very high when Phii reaches 0.95 to render further settlement almost impossible
    eta = (eta_0* rhoi * Phiiconst / c_eta) * np.exp(a_eta * (Tfus - Tconst) + b_eta * rhoi * Phiiconst) * restrict
    return eta

