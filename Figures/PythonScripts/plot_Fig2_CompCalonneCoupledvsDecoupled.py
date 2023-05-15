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
formatter.set_powerlimits((-4,4))

##To solve a bug on power representation y-axis Fig 3
formatter2 = ticker.ScalarFormatter(useMathText=True)
formatter2.set_scientific(True)
formatter2.set_powerlimits((-4,4))
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
ax.set_ylim([0.0,0.006])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)
ax.yaxis.set_major_formatter(formatter)

ax=axes[2]
ax.set_ylabel(r'$c$ ($\mathrm{kg~m^{-3}~s^{-1}}$)',fontsize=28)
ax.set_xlabel(r'$z$ (m)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([-20e-6,20e-6])
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
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig2/.')
    
    ###List all cases to plot
    Case_list = ['CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CD_PC', 'CD_PC', 'CD_PC','CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CD_PC', 'CD_PC', 'CD_PC']
    # Case_list = ['CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CC_2DOFs', 'CC_2DOFs', 'CC_2DOFs','CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CC_2DOFs', 'CC_2DOFs', 'CC_2DOFs']
    TspSize_list = ['15min', '15min', '15min', '15min', '15min', '15min','5min', '5min', '5min', '5min', '5min', '5min']

    dt_list_MyCode = [8, 16, 96, 8, 16, 96, 24, 48, 288, 24, 48, 288]
    color_list_MyCode = ['limegreen', 'limegreen', 'limegreen','darkred','darkred','darkred',(42/255, 127/255, 255/255), (42/255, 127/255, 255/255), (42/255, 127/255, 255/255),'orange','orange','orange']#, 'r'] ###legend for output of my code 
    Mylinesstyle_list = [':', '--', '-', ':', '--', '-', ':', '--', '-', ':', '--', '-']
    Mylinewidth_list = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
    
    ###Load outputs corresponding to considered cases 
    dict_Mycode_Nodal_list = [] ##Liste des données csv chargées (une par cas envisagé)
    for (Case, TspSize, dt) in zip(Case_list, TspSize_list, dt_list_MyCode):
        Depository_name = 'Simu_Papier_Fig2_{}_ProperMass_cLumped_MassnotLumped_CalPar_Alpha5minus3_SettlementOff_201nodes_24h_Output2h_Implicit_dt{}_'.format(Case, TspSize)
        file_name_nodal =  '{}_NODAL_dt_{}.csv'.format(Depository_name,dt)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name) 
        file = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        dict_Mycode_Nodal_list.append(file)
           
    ###Plot the considered cases
    ax = axes[0]#On Subplot0 we plot Temperature
    for i in range(len(Case_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['Temperature'], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i])
    #### Dummy plot for legend
    ax.plot(np.NaN,np.NaN, label= 'CC_3DOF; $\Delta t =15~\mathrm{min}$', color='limegreen',linewidth= 2.5 , linestyle= '-')
    ax.plot(np.NaN,np.NaN, label= 'CD_PC; $\Delta t =15~\mathrm{min}$', color='darkred',linewidth= 2.5 , linestyle= '-')
    ax.plot(np.NaN,np.NaN, label= 'CC_3DOF; $\Delta t =5~\mathrm{min}$', color=(42/255, 127/255, 255/255),linewidth= 2.5 , linestyle= '-')
    ax.plot(np.NaN,np.NaN, label= 'CD_PC; $\Delta t =5~\mathrm{min}$', color='orange',linewidth= 2.5 , linestyle= '-')
    ax.plot(np.NaN,np.NaN, label= '$t = 2~\mathrm{h}$', color='darkgrey',linewidth= 2.5 , linestyle= ':')
    ax.plot(np.NaN,np.NaN, label= '$t = 4~\mathrm{h}$', color='darkgrey',linewidth= 2.5 , linestyle= '--')
    ax.plot(np.NaN,np.NaN, label= '$t = 24~\mathrm{h}$', color='darkgrey',linewidth= 2.5 , linestyle= '-')
    ax.legend(loc='upper right', fontsize=18, ncol=3)

    ax = axes[1] #On Subplot1 we plot Rhov
    for i in range(len(Case_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['Rhov'], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; $\Delta t =15~\mathrm{min}$', color='limegreen', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; $\Delta t =15~\mathrm{min}$', color='darkred', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; $\Delta t =5~\mathrm{min}$', color=(42 / 255, 127 / 255, 255 / 255), linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; $\Delta t =5~\mathrm{min}$', color='orange', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$t = 2~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='$t = 4~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle='--')
    ax.plot(np.NaN, np.NaN, label='$t = 24~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle='-')
    ax.legend(loc='upper right', fontsize=18, ncol=3)

    ax = axes[2] #On Subplot2 we plot deposition rate
    for i in range(len(Case_list)):
        ax.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['DepositionRate'], color=color_list_MyCode[i],linewidth=  Mylinewidth_list[i] , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; $\Delta t =15~\mathrm{min}$', color='limegreen', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; $\Delta t =15~\mathrm{min}$', color='darkred', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; $\Delta t =5~\mathrm{min}$', color=(42 / 255, 127 / 255, 255 / 255), linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; $\Delta t =5~\mathrm{min}$', color='orange', linewidth=2.5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$t = 2~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='$t = 4~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle='--')
    ax.plot(np.NaN, np.NaN, label='$t = 24~\mathrm{h}$', color='darkgrey', linewidth=2.5, linestyle='-')
    ax.legend(loc='lower center', fontsize=18, ncol=3)

    ################################################################################
    # SAVE THE FIGURE #####
    ################################################################################
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig2_Paper/.')
    name_output_fig = 'Fig_2_Papier_Main_WithLegend'
    fig.savefig(path_output_fig.joinpath(name_output_fig))

    
    ################################################################################
    # ZOOM IN FIG c  #####
    ################################################################################
    ### Set-Up the Fig 
    fig2 = plt.figure(figsize=(15,15))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
    ax2 = fig2.add_subplot(spec2[0, 0])
    ax2.set_xlim([0.03,0.1])
    ax2.set_ylim([-10e-6,10e-6])
    ax2.tick_params(labelsize=34)    # fontsize of the tick labels
    ax2.grid(True)
    ax2.yaxis.set_major_formatter(formatter)
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax2.spines[axis].set_linewidth(8)
    # increase tick width
    ax2.tick_params(width=8)
    
    for i in range(len(Case_list)):
        ax2.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['DepositionRate'], color=color_list_MyCode[i],linewidth=  8 , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)
      
    ### Set-Up the Fig 
    fig3 = plt.figure(figsize=(15,15))
    spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
    ax3 = fig3.add_subplot(spec3[0, 0])
    ax3.set_xlim([0.63,0.88])
    ax3.set_ylim([-2.8e-6,2.8e-6])
    ax3.tick_params(labelsize=34)    # fontsize of the tick labels
    ax3.grid(True)
    ax3.yaxis.set_major_formatter(formatter2)
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax3.spines[axis].set_linewidth(8)
    # increase tick width
    ax3.tick_params(width=8)
    for i in range(len(Case_list)):
        ax3.plot(dict_Mycode_Nodal_list[i]['z'], dict_Mycode_Nodal_list[i]['DepositionRate'], color=color_list_MyCode[i],linewidth=  8 , linestyle= Mylinesstyle_list[i]) #, marker=marker,markersize=markersize)

    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################  
    name_output_fig2 = 'Fig_2_Papier_ZoomIn_BottomcPeak_WithLegend'
    fig2.savefig(path_output_fig.joinpath(name_output_fig2))
    
    name_output_fig3 = 'Fig_2_Papier_ZoomIn_TopcPeak_WithLegend'
    fig3.savefig(path_output_fig.joinpath(name_output_fig3))

    plt.show()