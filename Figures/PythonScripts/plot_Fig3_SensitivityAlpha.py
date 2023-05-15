################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt

from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

import pandas as pd ###To treat the csv files

from matplotlib import ticker
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True) 
formatter.set_powerlimits((-1,1)) 

###Defined Colors optimized for color-blind people:
Orange = [230/255, 159/255, 0/255]
SkyBlue = [86/255, 180/255, 233/255]
BluishGreen = [0/255, 158/255, 115/255]
Yellow = [240/255, 228/255, 66/255]
Blue = [0/255, 114/255, 178/255]
Vermillion = [213/255, 94/255, 0/255]
ReddishPurple= [204/255, 121/255, 167/255]
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN PART OF THE CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####             
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
#plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)

fig, axes = plt.subplots(3,1, figsize=(17,20), gridspec_kw={'height_ratios': [2,2,3]})
# pimp style
ax=axes[0]
ax.set_ylabel(r'$T$ (K)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([250,280])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)
ax.yaxis.set_major_formatter(formatter)
   
ax=axes[1]
ax.set_ylabel(r'$\rho_\mathrm{v}$ ($\mathrm{kg~m^{-3}}$)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([0.0,0.005])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)
ax.yaxis.set_major_formatter(formatter)

ax=axes[2]
ax.set_ylabel(r'$c$ ($\mathrm{kg~m^{-3}~s^{-1}}$)',fontsize=28)
ax.set_xlabel(r'$z$ (m)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([-7e-6,7e-6])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)
ax.yaxis.set_major_formatter(formatter)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig3/.')
    
    ###List all cases to plot
    AlphaCase_list = ['Alpha0', 'Alpha10minus8', 'Alpha10minus6', 'Alpha10minus4', 'Alpha10minus1', 'NoAlphaMixed']
    dt_list_MyCode = [152, 152, 152, 152, 152, 152]
    color_list_MyCode = ['grey', 'limegreen',  (42/255, 127/255, 255/255), 'orchid','darkorange',  'k']
    Mylinesstyle_list = ['-', '-', '-', '-', '-', ':']
    Mylinewidth_list = [3, 3, 3, 3, 3, 2.9]
    Mylegend_list = [r'CC_3DOF; $\alpha = 0$', r'CC_3DOF; $\alpha = 10^{-8}$',  r'CC_3DOF; $\alpha = 10^{-6}$',  r'CC_3DOF; $\alpha = 10^{-4}$',  r'CC_3DOF; $\alpha = 10^{-1}$', r'H_MF; No $\alpha$']
    
    ###Load outputs corresponding to considered cases 
    dict_Mycode_Nodal_list = [] ##Liste des données csv chargées (une par cas envisagé)
    for (AlphaCase, dt) in zip(AlphaCase_list, dt_list_MyCode):
        if AlphaCase == 'NoAlphaMixed':
            Depository_name = 'Simu_Papier_Fig3_H_MF_ProperMass_cLumped_CalPar_NoAlpha_SettlementOff_201nodes_38h_Implicit_dt15min_'
        else:
            Depository_name = 'Simu_Papier_Fig3_CC_3DOFs_ProperMass_cLumped_CalPar_{}_SettlementOff_201nodes_38h_Implicit_dt15min_NoFluxVap_'.format(AlphaCase)
        file_name_nodal =  '{}_NODAL_dt_{}.csv'.format(Depository_name,dt)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name) 
        file = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        dict_Mycode_Nodal_list.append(file)
           
    ###Plot the considered cases
    ax = axes[0]#On Subplot0 we plot Temperature
    for i in range(len(AlphaCase_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['Temperature'], label= Mylegend_list[i], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)    
        ax.legend(loc='upper right',fontsize = 18, ncol=3)

    ax = axes[1] #On Subplot1 we plot Rhov
    for i in range(len(AlphaCase_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['Rhov'], label= Mylegend_list[i], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)
        ax.legend(loc='upper right',fontsize = 18, ncol=3)

    ax = axes[2] #On Subplot2 we plot deposition rate
    for i in range(len(AlphaCase_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['DepositionRate'], label= Mylegend_list[i], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)
        ax.legend(loc='upper left', fontsize=18, ncol=3)

    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    name_output_fig = 'Fig_3_Papier_WithLegend'
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig3_Paper/.')
    fig.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()
