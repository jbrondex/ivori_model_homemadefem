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
from matplotlib.ticker import ScalarFormatter
# formatter = ticker.ScalarFormatter(useMathText=True)
# formatter.set_scientific(True)
# formatter.set_powerlimits((-1,1))

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
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)

fig, axes = plt.subplots(1,2, figsize=(40,20))
# pimp style
ax=axes[0]
ax.set_ylabel(r'$\Delta E_\Omega$ ($\mathrm{J~m^{-2}}$)',fontsize=50)
ax.set_xlabel(r't (days)',fontsize=50)
ax.set_xlim([0.0,5.0])
ax.set_ylim([-302,5])
ax.tick_params(labelsize=46)    # fontsize of the tick labels
ax.grid(True)
   
ax=axes[1]
ax.set_xlabel(r't (days)',fontsize=50)
ax.set_xlim([0.0,5.0])
ax.set_ylim([-302,5])
ax.tick_params(labelsize=46)    # fontsize of the tick labels
ax.grid(True)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig4/.')
    
    ###List all cases to plot for the Calonne system
    Case_list_Calonne = ['CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CC_3DOFs', 'CD_PC', 'CD_PC', 'CD_PC', 'CD_PC']
    TspSize_list_Calonne = ['15min', '5min', '15min', '5min', '15min', '5min', '15min', '5min']
    Dep_list_Calonne = ['DepOff', 'DepOff', 'DepOn', 'DepOn', 'DepOff', 'DepOff', 'DepOn', 'DepOn']
    color_list_MyCode_Calonne = ['limegreen', 'limegreen', 'darkred','darkred', (42/255, 127/255, 255/255), (42/255, 127/255, 255/255), 'orange','orange']
    Mylinesstyle_list_Calonne = ['-', ':', '-', ':', '-', ':', '-', ':']
    Mylinewidth_list_Calonne = [5, 5, 5, 5, 5, 5, 5, 5]
    Mylegend_list_Calonne = [r'CC_3DOFs; $\Delta t = 15mn$; DepOff', r'CC_3DOFs; $\Delta t = 5mn$; DepOff', r'CC_3DOFs; $\Delta t = 15mn$; DepOn', r'CC_3DOFs; $\Delta t = 5mn$; DepOn', r'CD_PC; $\Delta t = 15mn$; DepOff $', r'CD_PC; $\Delta t = 5mn$; DepOff $', r'CD_PC; $\Delta t = 15mn$; DepOn $', r'CD_PC; $\Delta t = 5mn$; DepOn $']
    
    ###Load outputs corresponding to considered cases 
    dict_Mycode_Nodal_list_Calonne = [] ##Liste des données csv chargées (une par cas envisagé)
    for (Case, TspSize, Dep) in zip(Case_list_Calonne, TspSize_list_Calonne, Dep_list_Calonne):
        Depository_name = 'Simu_Papier_Fig4_{}_ProperMass_cLumped_MassnotLumped_CalPar_Alpha5minus3_NoFluxBCs_SettlementOff_{}_201nodes_5d_Output12h_Implicit_dt{}_'.format(Case, Dep, TspSize)  
        file_name_nodal =  '{}_ENERGY.csv'.format(Depository_name)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name) 
        file = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        dict_Mycode_Nodal_list_Calonne.append(file)
           
    ###Plot the considered cases
    ax = axes[0]#On Subplot0 we plot all Calonne Cases
    for i in range(len(Case_list_Calonne)):
        ax.plot(dict_Mycode_Nodal_list_Calonne[i]['time']/(24*3600), dict_Mycode_Nodal_list_Calonne[i]['Energy Leak from beginning'], color=color_list_MyCode_Calonne[i],linewidth=  Mylinewidth_list_Calonne[i] , linestyle= Mylinesstyle_list_Calonne[i])
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; Dep. Off ', color='limegreen', linewidth=5,linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; Dep. Off', color=(42/255, 127/255, 255/255), linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CC_3DOF; Dep. On', color='darkred', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='CD_PC; Dep. On', color='orange', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =5~\mathrm{min}$', color='darkgrey', linewidth=5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =15~\mathrm{min}$', color='darkgrey', linewidth=5, linestyle='-')
    ax.legend(loc='lower left', fontsize=26, ncol=2)

    ###List all cases to plot for the Hansen system
    Case_list_Hansen = ['H_MF', 'H_MF', 'H_MF', 'H_MF', 'H_TF', 'H_TF', 'H_TF', 'H_TF']
    TspSize_list_Hansen = ['15min', '5min', '15min', '5min', '15min', '5min', '15min', '5min']
    Dep_list_Hansen = ['DepOff', 'DepOff', 'DepOn', 'DepOn', 'DepOff', 'DepOff', 'DepOn', 'DepOn']
    color_list_MyCode_Hansen = ['limegreen', 'limegreen', 'darkred', 'darkred',(42 / 255, 127 / 255, 255 / 255), (42 / 255, 127 / 255, 255 / 255), 'orange','orange']
    Mylinesstyle_list_Hansen = ['-', ':', '-', ':', '-', ':', '-', ':']
    Mylinewidth_list_Hansen = [5, 5, 5, 5, 5, 5, 5, 5]
    Mylegend_list_Hansen = [r'H_MF; $\Delta t = 15mn$; DepOff', r'H_MF; $\Delta t = 5mn$; DepOff',r'H_MF; $\Delta t = 15mn$; DepOn', r'H_MF; $\Delta t = 5mn$; DepOn',r'H_TF; $\Delta t = 15mn$; DepOff', r'H_TF; $\Delta t = 5mn$; DepOff',r'H_TF; $\Delta t = 15mn$; DepOn', r'H_TF; $\Delta t = 5mn$; DepOn']
    ###Load outputs corresponding to considered cases
    dict_Mycode_Nodal_list_Hansen = []  ##Liste des données csv chargées (une par cas envisagé)
    for (Case, TspSize, Dep) in zip(Case_list_Hansen, TspSize_list_Hansen, Dep_list_Hansen):
        Depository_name = 'Simu_Papier_Fig4_{}_ProperMass_cLumped_MassnotLumped_CalPar_NoFluxBCs_SettlementOff_{}_201nodes_5d_Output12h_Implicit_dt{}_'.format(Case, Dep, TspSize)
        file_name_nodal = '{}_ENERGY.csv'.format(Depository_name)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name)
        file = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        dict_Mycode_Nodal_list_Hansen.append(file)

    ###Plot the considered cases
    ax = axes[1]  # On Subplot1 we plot all Hansen Cases
    for i in range(len(Case_list_Hansen)):
        ax.plot(dict_Mycode_Nodal_list_Hansen[i]['time'] / (24 * 3600),dict_Mycode_Nodal_list_Hansen[i]['Energy Leak from beginning'],color=color_list_MyCode_Hansen[i], linewidth=Mylinewidth_list_Hansen[i], linestyle=Mylinesstyle_list_Hansen[i])
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='H_MF; Dep. Off ', color='limegreen', linewidth=5,linestyle='-')
    ax.plot(np.NaN, np.NaN, label='H_TF; Dep. Off', color=(42/255, 127/255, 255/255), linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='H_MF; Dep. On', color='darkred', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='H_TF; Dep. On', color='orange', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =5~\mathrm{min}$', color='darkgrey', linewidth=5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='$\Delta t =15~\mathrm{min}$', color='darkgrey', linewidth=5, linestyle='-')
    ax.legend(loc='lower left', fontsize=26, ncol=2)
    ################################################################################
    # SAVE THE FIGURE #####
    ################################################################################
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig4_Paper')
    name_output_fig = 'Fig_4_Papier_Main_WithLegend'
    fig.savefig(path_output_fig.joinpath(name_output_fig))

    ################################################################################
    # ZOOM IN FIG CALONNE  #####
    ################################################################################
    ### Set-Up the Fig
    fig2 = plt.figure(figsize=(15, 15))
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
    ax2 = fig2.add_subplot(spec2[0, 0])
    ax2.set_xlim([0.0, 5.0])
    ax2.set_ylim([-8, 1])
    ax2.tick_params(labelsize=48)  # fontsize of the tick labels
    ax2.grid(True)
    ax2.yaxis.set_major_formatter(ScalarFormatter())
    # change all spines
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(7)
    # increase tick width
    ax2.tick_params(width=15)

    for i in range(len(Case_list_Calonne)):
        ax2.plot(dict_Mycode_Nodal_list_Calonne[i]['time'] / (24 * 3600), dict_Mycode_Nodal_list_Calonne[i]['Energy Leak from beginning'], label=Mylegend_list_Calonne[i], color=color_list_MyCode_Calonne[i], linewidth=7, linestyle=Mylinesstyle_list_Calonne[i])

    ### Set-Up the Fig
    fig3 = plt.figure(figsize=(15, 15))
    spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
    ax3 = fig3.add_subplot(spec3[0, 0])
    ax3.set_xlim([0.0, 5.0])
    ax3.set_ylim([-0.4, 0.4])
    ax3.tick_params(labelsize=48)  # fontsize of the tick labels
    ax3.grid(True)
    ax3.yaxis.set_major_formatter(ScalarFormatter())
    # change all spines
    for axis in ['top', 'bottom', 'left', 'right']:
        ax3.spines[axis].set_linewidth(7)
    # increase tick width
    ax3.tick_params(width=15)

    for i in range(len(Case_list_Hansen)):
        ax3.plot(dict_Mycode_Nodal_list_Hansen[i]['time'] / (24 * 3600),dict_Mycode_Nodal_list_Hansen[i]['Energy Leak from beginning'], label=Mylegend_list_Hansen[i],color=color_list_MyCode_Hansen[i], linewidth=7, linestyle=Mylinesstyle_list_Hansen[i])

    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    name_output_fig2 = 'Fig_4_Papier_ZoomIn_Calonne_WithLegend'
    fig2.savefig(path_output_fig.joinpath(name_output_fig2))

    name_output_fig3 = 'Fig_4_Papier_ZoomIn_Hansen_WithLegend'
    fig3.savefig(path_output_fig.joinpath(name_output_fig3))

    plt.show()