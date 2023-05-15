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

################################################################################
# prepare plots #####
################################################################################
# global figure pimping
#plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('axes', labelsize=20)
plt.rc('legend', fontsize=24)

### Set-Up the Fig for mass evolution and Settling Rate
fig1, axes1 = plt.subplots(1,2, figsize=(40,20))
# pimp style
ax1=axes1[0]
ax1.set_ylabel(r'Total Mass ($\mathrm{kg~m^{-2}}$)',fontsize=50)
ax1.set_xlabel(r't (days)',fontsize=50)
ax1.set_xlim([0, 20])
ax1.xaxis.set_ticks(np.linspace(0, 20, 11))
ax1.set_ylim([56.225,56.255])
ax1.tick_params(labelsize=46)    # fontsize of the tick labels
ax1.grid(True)
   
ax1=axes1[1]
ax1.set_ylabel(r'$|v_\mathrm{set}|$ ($\mathrm{cm~d^{-1}}$)',fontsize=50)
ax1.set_xlabel(r'$z$ (m)',fontsize=50)
ax1.set_xlim([0.0,0.5])
ax1.set_ylim([0,4.5])
ax1.tick_params(labelsize=46)    # fontsize of the tick labels
ax1.grid(True)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~    FOR SIMSON MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
    pathroot = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig5/.')
    ###List all cases to plot
    Mesh = ['101', '101', '51', '51', '11', '11']
    Case_name_list = ['AllBugsCorrected_PhiiImplicit', 'AllBugsCorrected_PhiiExplicit', 'AllBugsCorrected_PhiiImplicit', 'AllBugsCorrected_PhiiExplicit', 'AllBugsCorrected_PhiiImplicit', 'AllBugsCorrected_PhiiExplicit']
    Color_Per_Case = ['orange', 'orange', (42/255, 127/255, 255/255), (42/255, 127/255, 255/255), 'darkred', 'darkred']
    Data_list=['all_c', 'all_coord', 'all_dz', 'all_phi', 'all_rho_v', 'all_T', 'all_t_passed', 'all_v']
    Key_list=['DepositionRate', 'z', 'dz', 'Phii', 'Rhov', 'Temperature', 'time', 'Velocity']

    Time_list = [24*1, 24*2, 24*5, 24*20] ### Time of plot in hours
    Legend_list = ['Simson, t=1d', 'Simson, t=2d','Simson, t=5d','Simson, t=20d']
    Linesstyle_list = [':','--','-.','-']

    dict_Simson_list = [] ### Contains a list of dict. One dict per case of Case_name_list. Each dict contains all data of data_list
    TotalMass_PieceWiseFromBotNode_list = [] ### For each case of Case_name_list, contains total mass for each tsp calculated as piecewise constant from phii of bottom node
    for (mesh,case) in zip(Mesh, Case_name_list):
        depository_name =  'Output_SimsonCode/Simu_SimsonModel_Case6_Fig5_{}nodes_dt15min_{}'.format(mesh,case)
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
        ########For each of the cases, I make additionnal calculation to check mass conservation############
        Mass_per_layer_PieceWiseFromBotNode = 917*dict_Simson['Phii'][:,:-1]*dict_Simson['dz']
        TotalMass_PieceWiseFromBotNode=[]
        for i in range(np.shape(Mass_per_layer_PieceWiseFromBotNode)[0]):
            TotalMass_PieceWiseFromBotNode.append(np.sum(Mass_per_layer_PieceWiseFromBotNode[i,:]))
        TotalMass_PieceWiseFromBotNode_list.append(TotalMass_PieceWiseFromBotNode)

    ###################################################################
    # MAIN FIGURES THAT ARE STARTED WITH ANNA SIMSON CODE RESULTS #####
    ###################################################################
    ###Plot the mass evolution (subplot 0)
    ax1 = axes1[0]  #On Subplot0 we plot mass over time
    for i, case in enumerate(Case_name_list):# [1,3,5]:
        ax1.plot(dict_Simson_list[i]['time'] / (24 * 3600), TotalMass_PieceWiseFromBotNode_list[i],color=Color_Per_Case[i], linewidth=5, linestyle='-')
    ###Plot the settling velo (subplot 1)
    ax1 = axes1[1] #On Subplot1 we plot settling velo vs z
    for i in [4]: ## We do the plot only for the mesh with 11nodes with Phii implicit (i.e., Simson corrected)
        for idx,timeh in enumerate(Time_list):
            ax1.plot(dict_Simson_list[i]['z'][dict_Simson_list[i]['t{}h_index'.format(timeh)],:], abs(dict_Simson_list[i]['Velocity'][dict_Simson_list[i]['t{}h_index'.format(timeh)],:])*60*60*24*100, color=Color_Per_Case[i],linewidth=  5 , linestyle= Linesstyle_list[idx])

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig5/.')
    ###List all cases to plot
    Mesh = ['11'] ## We do the calculation of settling velocities only for the mesh with 11nodes (mass is perfectly conserved for all meshes)
    Time_list_MyCode = [96, 192, 480, 1920]
    MyColor_list_PerMesh = ['limegreen']
    MyLinesstyle_list = [':','--','-.','-']
    MyLinewidth_list = [5, 5, 5, 5]

    Total_Mass_Mycode_101nodes =[0.25*150+0.25*75] ##Initial mass
    Total_Mass_Mycode_51nodes =[0.25*150+0.25*75]  ##Initial mass
    Total_Mass_Mycode_11nodes =[0.25*150+0.25*75]  ##Initial mass
    Znodes_1d_list=[]
    Znodes_2d_list=[]
    Znodes_5d_list=[]
    Znodes_20d_list=[]
    Vset_1d_mpersec_list=[]
    Vset_2d_mpersec_list=[]
    Vset_5d_mpersec_list=[]
    Vset_20d_mpersec_list=[]
    for mesh in Mesh:
        Depository_name = 'Simu_Papier_Fig5_SimsonCase6_SettlementOnly_{}nodes_10d_Output1d_Implicit_dt15min_'.format(mesh)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name)
        for day in range(1,21):
            Time_list_MyCode_MassCons = day * 24 * 4 ## Get timestepnumber
            file_name_elemental = '{}_ELEMENTAL_dt_{}.csv'.format(Depository_name,Time_list_MyCode_MassCons)
            file_elemental = pd.read_csv(full_path_mycode.joinpath(file_name_elemental))
            if mesh == '101':
                Total_Mass_Mycode_101nodes.append(np.sum(file_elemental['Phii']*917*(file_elemental['ztop']-file_elemental['zbot'])))
            elif mesh =='51':
                Total_Mass_Mycode_51nodes.append(np.sum(file_elemental['Phii']*917*(file_elemental['ztop']-file_elemental['zbot'])))
            else:
                Total_Mass_Mycode_11nodes.append(np.sum(file_elemental['Phii']*917*(file_elemental['ztop']-file_elemental['zbot'])))
        dict_Mycode_Case6_Nodal_list=[]
        dict_Mycode_Case6_Nodal_list_prev=[]
        for time in Time_list_MyCode:
            file_name_nodal = '{}_NODAL_dt_{}.csv'.format(Depository_name, time)
            file_name_nodal_prev = '{}_NODAL_dt_{}.csv'.format(Depository_name, time-1)###We need node positions at tsp right before tsp considered to calculate settling velocity
            file_nodal = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
            file_nodal_prev = pd.read_csv(full_path_mycode.joinpath(file_name_nodal_prev))
            dict_Mycode_Case6_Nodal_list.append(file_nodal)
            dict_Mycode_Case6_Nodal_list_prev.append(file_nodal_prev)
        Znodes_1d_list.append(dict_Mycode_Case6_Nodal_list[0]['z'])
        Vset_1d_mpersec_list.append((dict_Mycode_Case6_Nodal_list[0]['z']-dict_Mycode_Case6_Nodal_list_prev[0]['z'])/(15*60)) ##we divide travelled distance per time step size in seconds (15mn*60s)
        Znodes_2d_list.append(dict_Mycode_Case6_Nodal_list[1]['z'])
        Vset_2d_mpersec_list.append((dict_Mycode_Case6_Nodal_list[1]['z']-dict_Mycode_Case6_Nodal_list_prev[1]['z'])/(15*60)) ##we divide travelled distance per time step size in seconds (15mn*60s)
        Znodes_5d_list.append(dict_Mycode_Case6_Nodal_list[2]['z'])
        Vset_5d_mpersec_list.append((dict_Mycode_Case6_Nodal_list[2]['z']-dict_Mycode_Case6_Nodal_list_prev[2]['z'])/(15*60)) ##we divide travelled distance per time step size in seconds (15mn*60s)
        Znodes_20d_list.append(dict_Mycode_Case6_Nodal_list[3]['z'])
        Vset_20d_mpersec_list.append((dict_Mycode_Case6_Nodal_list[3]['z']-dict_Mycode_Case6_Nodal_list_prev[3]['z'])/(15*60)) ##we divide travelled distance per time step size in seconds (15mn*60s)
    ###############################################################################
    # MAIN FIGURES THAT ARE FILLED WITH MY CODE RESULTS #####
    ################################################################################
    ax1 = axes1[0]#On Subplot0 we plot mass over time
    # ax1.plot(np.linspace(0,20,21), Total_Mass_Mycode_101nodes, color=MyColor_list_PerMesh[0],linewidth=1.7, linestyle='-')  # , marker=marker,markersize=markersize)
    # ax1.plot(np.linspace(0,20,21), Total_Mass_Mycode_51nodes, color=MyColor_list_PerMesh[1],linewidth=1.7, linestyle='-')
    ax1.plot(np.linspace(0,20,21), Total_Mass_Mycode_11nodes, color=MyColor_list_PerMesh[0],linewidth=5, linestyle='-')
    #### Dummy plot for legend
    ax1.plot(np.NaN, np.NaN, label='Implicit', color='limegreen', linewidth=5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='Explicit; $N_z = 11$', color='darkred', linewidth=5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='Explicit; $N_z = 51$', color=(42/255, 127/255, 255/255), linewidth=5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='Explicit; $N_z = 101$', color='orange', linewidth=5, linestyle='-')
    ax1.legend(loc='upper left', fontsize=36, ncol=2)

    ax1 = axes1[1]#On Subplot1 we plot settling velo vs z
    for i in range(len(Mesh)):
        ax1.plot(Znodes_1d_list[i], abs(Vset_1d_mpersec_list[i]) * 60 * 60 * 24 * 100, color=MyColor_list_PerMesh[i],linewidth=5, linestyle=':')
        ax1.plot(Znodes_2d_list[i], abs(Vset_2d_mpersec_list[i]) * 60 * 60 * 24 * 100, color=MyColor_list_PerMesh[i],linewidth=5, linestyle='--')
        ax1.plot(Znodes_5d_list[i], abs(Vset_5d_mpersec_list[i]) * 60 * 60 * 24 * 100, color=MyColor_list_PerMesh[i],linewidth=5, linestyle='-.')
        ax1.plot(Znodes_20d_list[i], abs(Vset_20d_mpersec_list[i]) * 60 * 60 * 24 * 100, color=MyColor_list_PerMesh[i],linewidth=5, linestyle='-')
    #### Dummy plot for legend
    ax1.plot(np.NaN, np.NaN, label='Our code', color='limegreen', linewidth=5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='Simson code', color='darkred', linewidth=5, linestyle='-')
    ax1.plot(np.NaN, np.NaN, label='t = 1 d', color='darkgrey', linewidth=5, linestyle=':')
    ax1.plot(np.NaN, np.NaN, label='t = 2 d', color='darkgrey', linewidth=5, linestyle='--')
    ax1.plot(np.NaN, np.NaN, label='t = 5 d', color='darkgrey', linewidth=5, linestyle='-.')
    ax1.plot(np.NaN, np.NaN, label='t = 20 d', color='darkgrey', linewidth=5, linestyle='-')
    ax1.legend(loc='upper left', fontsize=36, ncol=2)

    ################################################################################
    # SAVE THE FIGURE #####
    ################################################################################
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig5_Paper')
    name_output_fig = 'Fig_5_Papier_Vset11nodes_WithLegend'
    print(name_output_fig)
    fig1.savefig(path_output_fig.joinpath(name_output_fig))

    plt.show()