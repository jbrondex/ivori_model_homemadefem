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
    g
    )

def USF_Tinit(z):
    Tinit_top = 253
    Tinit_bot = 273
    Htot = 0.5 
    
    Tinit = (Tinit_top-Tinit_bot)/Htot * z + Tinit_bot
    return Tinit 

def USF_TinitSchurholtScenario2(z):
    Tinit_top = 253
    Tinit_bot = 273
    Htot = 1 
    
    Tinit = (Tinit_top-Tinit_bot)/Htot * z + Tinit_bot
    return Tinit 

def USF_TinitSchurholtScenario3(z):
    Tinit_top = 253
    Tinit_bot = 273
    Htot = 0.02 
    
    Tinit = (Tinit_top-Tinit_bot)/Htot * z + Tinit_bot
    return Tinit


def USF_PhiiinitSimson(z):
    if 0 <= z <= 0.24:
        rho_eff = 150 ##rho_eff [kg/m^3]
    elif z >= 0.26:
        rho_eff = 75 ##rho_eff [kg/m^3]
    else:
        rho_eff = -3750 * z + 1050 ### smooth transition over 2 cm between the two layers (see Simson 2021)  
    Phii = rho_eff / rhoi 
    return Phii

def USF_PhiiinitSimsonSharpTransition(z):
    if 0 <= z < 0.25:
        rho_eff = 150 ##rho_eff [kg/m^3]
    elif z > 0.25:
        rho_eff = 75 ##rho_eff [kg/m^3]
    elif z == 0.25:
        rho_eff = (150+75)/2 ### sharp transition between the two layers (see Simson 2021)
    Phii = rho_eff / rhoi
    return Phii

def USF_PhiiinitSchurholtScenario2(z):
    if 0 <= z <= 0.08:
        Phii = 1-9.2425*z
    elif 0.08 <= z <= 0.64:
        Phii = 0.2606
    elif 0.64 <= z <= 0.72:
        Phii = 0.2606 + 4.915*(z-0.64)
    elif 0.72 <= z <= 0.75:
        Phii = 0.6538
    elif 0.75 <= z <= 0.86:
        Phii = 0.6538 - 4.915*(z-0.75335)
    elif 0.86 <= z <= 2.0:
        Phii = 0.1295895
    return Phii


def USF_PhiiinitSchurholtScenario3(z):
    sigma_sqrd = 0.0000005
    mu         = 0.01
    offset     = 0.3
    corr       = 0.00177
    pi         = 3.14159265359
    Phii  = offset + 0.2*corr/(np.sqrt(2*pi*sigma_sqrd))*np.exp(-((z - mu)**2)/(2*sigma_sqrd))
    return Phii

def USF_BCtopSchurholtscenario1(t):
    tau = 60*60*5 ### tau = 5h
    Ttop = 273 + (263-273)*(t/tau * np.heaviside(t, 1) - ((t - tau)/tau) * np.heaviside(t-tau, 1))
    return Ttop
