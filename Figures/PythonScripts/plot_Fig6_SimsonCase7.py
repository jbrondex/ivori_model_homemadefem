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

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)

### Set-Up the Fig for Phii, T, Rhov, c: Simson Case 7 only
fig1, axes1 = plt.subplots(2,1, figsize=(20,18))
# pimp style
ax1=axes1[0]
ax1.set_ylabel(r'$c$ ($\mathrm{kg~m^{3}~s^{-1}}$)',fontsize=32)
ax1.set_xlim([0.0,0.5])
ax1.set_ylim([-4.5e-6,4.5e-6])
ax1.tick_params(labelsize=26)    # fontsize of the tick labels
ax1.grid(True)
ax1.yaxis.set_major_formatter(formatter)

ax1=axes1[1]
ax1.set_ylabel(r'$\Phi_\mathrm{i}$ (-)',fontsize=32)
ax1.set_xlabel(r'$z$ (m)',fontsize=32)
ax1.set_xlim([0.0,0.5])
ax1.tick_params(labelsize=26)    # fontsize of the tick labels
ax1.grid(True)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~    FOR SIMSON MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
    pathroot = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig6')
    ###List all cases to plot
    Mesh = ['101']
    Case_name_list = ['AllBugsCorrected_PhiiImplicit']
    Color_Per_Case = ['k']

    Data_list=['all_c', 'all_coord', 'all_dz', 'all_phi', 'all_rho_v', 'all_T', 'all_t_passed', 'all_v']
    Key_list=['DepositionRate', 'z', 'dz', 'Phii', 'Rhov', 'Temperature', 'time', 'Velocity']

    Time_list = [24*1, 24*2, 24*5, 24*10] ### Time of plot in hours
    Legend_list = ['Simson, t=1d', 'Simson, t=2d','Simson, t=5d','Simson, t=10d']
    Linesstyle_list = [':','--','-.','-']

    dict_Simson_list = [] ### Contains a list of dict. One dict per case of Case_name_list. Each dict contains all data of data_list
    TotalMass_PieceWiseFromBotNode_list = [] ### For each case of Case_name_list, contains total mass for each tsp calculated as piecewise constant from phii of bottom node
    for (mesh,case) in zip(Mesh, Case_name_list):
        depository_name =  'Output_SimsonCode/Simu_SimsonModel_Case7_Fig6_{}nodes_dtVar_{}'.format(mesh,case)
        full_path = pathroot.joinpath(depository_name)
        dict_Simson = {}
        for (Data, Key) in zip(Data_list, Key_list):
            value = np.loadtxt(full_path.joinpath('{}_{}'.format(Data,mesh)))
            dict_Simson[Key] = value
            if Data == 'all_t_passed':
                for i in Time_list: #find index of line corresponding to time step of interest
                    print('hours of simu:',i, 'h')
                    timeh=i*3600    #i hours in seconds
                    length = len(value)
                    timeh_array = np.ones_like(length) * timeh
                    timeh_diff = np.absolute (timeh_array - value)
                    timeh_diff_list = list (timeh_diff)
                    timeh_index = timeh_diff_list.index(min(timeh_diff_list))
                    Var_name_2 = 't{}h_index'.format(i)
                    dict_Simson[Var_name_2] = timeh_index
        dict_Simson_list.append(dict_Simson)

    ################################################################################
    # MAIN FIGURES THAT ARE STARTED WITH ANNA SIMSON CODE RESULTS #####
    ################################################################################
    ###Plot the c profile (subplot 1)
    ax1 = axes1[0]  # On Subplot0 we plot c profiles
    for i in range(len(Case_name_list)):
        for idx, timeh in enumerate(Time_list):
            ax1.plot(dict_Simson_list[i]['z'][dict_Simson_list[i]['t{}h_index'.format(timeh)], :],dict_Simson_list[i]['DepositionRate'][dict_Simson_list[i]['t{}h_index'.format(timeh)], :],color=Color_Per_Case[i], linewidth=2.5, linestyle=Linesstyle_list[idx])
    #### Dummy plot for legend
    ax1.plot(np.NaN, np.NaN, label='Simson code; $c_0 = c _{N_z} = 0$', color='k', linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='H_MF; $c_0 = c _{N_z} = 0$', color='orange', linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='H_MF; No Flux BC on vap.', color=(42/255, 127/255, 255/255), linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='t = 1 d', color='darkgrey', linewidth=2.5, linestyle=':')
    ax1.plot(np.NaN, np.NaN, label='t = 2 d', color='darkgrey', linewidth=2.5, linestyle='--')
    ax1.plot(np.NaN, np.NaN, label='t = 5 d', color='darkgrey', linewidth=2.5, linestyle='-.')
    ax1.plot(np.NaN, np.NaN, label='t = 10 d', color='darkgrey', linewidth=2.5, linestyle='-')
    ax1.legend(bbox_to_anchor=(0.015, 1-0.025), loc='upper left', fontsize=22, ncol=2, borderaxespad=0)

    ###Plot the phii profile (subplot 1)
    ax1 = axes1[1]  #On Subplot1 we plot phii profiles
    for i in range(len(Case_name_list)):
        for idx,timeh in enumerate(Time_list):
            ax1.plot(dict_Simson_list[i]['z'][dict_Simson_list[i]['t{}h_index'.format(timeh)],:], dict_Simson_list[i]['Phii'][dict_Simson_list[i]['t{}h_index'.format(timeh)],:],color=Color_Per_Case[i], linewidth=2.5, linestyle= Linesstyle_list[idx])
    #### Dummy plot for legend
    ax1.plot(np.NaN, np.NaN, label='Simson code; $c_0 = c _{N_z} = 0$', color='k', linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='H_MF; $c_0 = c _{N_z} = 0$', color='orange', linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='H_MF; No Flux BC on vap.', color=(42/255, 127/255, 255/255), linewidth=2.5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='t = 1 d', color='darkgrey', linewidth=2.5, linestyle=':')
    ax1.plot(np.NaN, np.NaN, label='t = 2 d', color='darkgrey', linewidth=2.5, linestyle='--')
    ax1.plot(np.NaN, np.NaN, label='t = 5 d', color='darkgrey', linewidth=2.5, linestyle='-.')
    ax1.plot(np.NaN, np.NaN, label='t = 10 d', color='darkgrey', linewidth=2.5, linestyle='-')
    ax1.legend(loc='upper right', fontsize=22, ncol=2)
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig6/.')

    Case_list = ['BCVapNoFlux_', 'BCDepRate0_']#, 'BCNoFlux_Alpha1_']
    Time_list_MyCode = [96, 192, 480, 960]
    MyColor_list_PerCase = [ (42/255, 127/255, 255/255), 'orange']#, 'orchid']
    MyLinesstyle_list = [':','--','-.','-']
    MyLinewidth_list = [2.5, 2.5, 2.5, 2.5]

    for k, case in enumerate(Case_list):
        if k <= 1:
            Depository_name = 'Simu_Papier_Fig6_SimsonCase7_H_MF_SettlementOn_DepOn_101nodes_20d_Output1d_Implicit_dt15min_SmoothTransition_{}'.format(case)
        else:
            Depository_name = 'Simu_Papier_Fig6_SimsonCase7_CC_3DOFs_SettlementOn_DepOn_101nodes_20d_Output1d_Implicit_dt15min_SmoothTransition_{}'.format(case)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name)
        dict_Mycode_Case7_Nodal_list=[]
        dict_Mycode_Case7_Elemental_list=[]
        for time in Time_list_MyCode:
            file_name_nodal = '{}_NODAL_dt_{}.csv'.format(Depository_name, time)
            file_nodal = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
            dict_Mycode_Case7_Nodal_list.append(file_nodal)
            file_name_elemental = '{}_ELEMENTAL_dt_{}.csv'.format(Depository_name, time)
            file_elemental = pd.read_csv(full_path_mycode.joinpath(file_name_elemental))
            dict_Mycode_Case7_Elemental_list.append(file_elemental)
        ###############################################################################
        # MAIN FIGURES THAT ARE FILLED WITH MY CODE RESULTS #####
        ################################################################################
        ax1 = axes1[0]  # On Subplot0 we plot c profiles
        for idx, time in enumerate(Time_list_MyCode):
            ax1.plot(dict_Mycode_Case7_Nodal_list[idx]['z'], dict_Mycode_Case7_Nodal_list[idx]['DepositionRate'],color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[idx], linestyle=MyLinesstyle_list[idx])

        ax1 = axes1[1]#On Subplot1 we plot phii profiles
        for idx, time in enumerate(Time_list_MyCode):
            ax1.step(dict_Mycode_Case7_Elemental_list[idx]['zbot'], dict_Mycode_Case7_Elemental_list[idx]['Phii'],where='post', color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[idx],linestyle=MyLinesstyle_list[idx])

    ################################################################################
    # SAVE THE FIGURE #####
    ################################################################################
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig6_Paper/.')
    name_output_fig = 'Fig_6_Papier_WithoutCC3DOFs_WithLegend'
    fig1.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()