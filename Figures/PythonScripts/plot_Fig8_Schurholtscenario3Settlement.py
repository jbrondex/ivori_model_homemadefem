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

Orange = [230/255, 159/255, 0/255]
SkyBlue = [86/255, 180/255, 233/255]
BluishGreen = [0/255, 158/255, 115/255]
Yellow = [225/255, 190/255, 106/255]
Blue = [42/255, 127/255, 255/255]
Vermillion = [213/255, 94/255, 0/255]
ReddishPurple= [204/255, 121/255, 167/255]
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTIONS CALLED IN THE MAIN CODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####

################################################################################
# define import function for the Fenics model outputs #####
################################################################################
def getNodeArrays_Phii(vtu_filename):
    # read file
    print('Reading ' + vtu_filename + '...')
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_filename)
    reader.Update()
    data = reader.GetOutput()

    print('Number of points: ' + str(reader.GetNumberOfPoints()))
    print('Number of cells: ' + str(reader.GetNumberOfCells()))

    print('Number of point arrays: ' + str(reader.GetNumberOfPointArrays()))
    print('Number of cell arrays: ' + str(reader.GetNumberOfCellArrays()))

    for i in range(reader.GetNumberOfPointArrays()):
        print('Name of point array ' + str(i) + ': ' + str(reader.GetPointArrayName(i)))

    print('Returning numpy arrays for point data from point array 0')
    # get arrays
    points = data.GetPoints()
    x = vtk_to_numpy(points.GetData())[:, 0]

    usg = dsa.WrapDataObject(data)
    y = usg.PointData[reader.GetPointArrayName(0)]
    print('shape of x:', x.shape)
    print('shape of y:', y.shape)
    print('... done.\n')

    return x, y


def getNodeArrays_Mixed(vtu_filename, DoF):
    # read file
    print('Reading ' + vtu_filename + '...')
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_filename)
    reader.Update()
    data = reader.GetOutput()

    print('Number of points: ' + str(reader.GetNumberOfPoints()))
    print('Number of cells: ' + str(reader.GetNumberOfCells()))

    print('Number of point arrays: ' + str(reader.GetNumberOfPointArrays()))
    print('Number of cell arrays: ' + str(reader.GetNumberOfCellArrays()))

    for i in range(reader.GetNumberOfPointArrays()):
        print('Name of point array ' + str(i) + ': ' + str(reader.GetPointArrayName(i)))

    print('Returning numpy arrays for point data from point array 0')
    # get arrays
    points = data.GetPoints()
    x = vtk_to_numpy(points.GetData())[:, 0]

    usg = dsa.WrapDataObject(data)
    y = usg.PointData[reader.GetPointArrayName(0)][:, DoF]  ###DoF is 0 for vn, 1 for Rhov, 2 for T
    print('shape of x:', x.shape)
    print('shape of y:', y.shape)
    print('... done.\n')

    return x, y


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

### WE DO THREE TIMES THE SAME PLOT FOR THE THREE CONSIDERED RESOLUTION (1000, 350 and 100 nodes)
### Set-Up a list of the three fig
fig1000, axes1000 = plt.subplots(1, 2, figsize=(40, 20))
fig350, axes350 = plt.subplots(1, 2, figsize=(40, 20))
fig100, axes100 = plt.subplots(1, 2, figsize=(40, 20))

###pimp style for the three fig
for axes in [axes1000, axes350, axes100]:
    # pimp style
    ###Deposition rate for Calonne System
    ax=axes[0]
    ax.set_ylabel(r'$c$ ($\mathrm{kg~m^{3}~s^{-1}}$)',fontsize=60)
    ax.set_xlabel(r'$z$ (cm)',fontsize=60)
    ax.set_xlim([0.0,2])
    ax.locator_params(nbins=5, axis='x')
    ax.set_ylim([-2.5e-2,2.5e-2])
    ax.tick_params(labelsize=50)    # fontsize of the tick labels
    ax.grid(True)
    ###Phii for Calonne System
    ax=axes[1]
    ax.set_ylabel(r'$\Phi_\mathrm{i}$ (-)',fontsize=60)
    ax.set_xlabel(r'$z$ (cm)',fontsize=60)
    ax.set_xlim([0.0,2])
    ax.locator_params(nbins=5, axis='x')
    ax.set_ylim([0.200,0.55])
    ax.tick_params(labelsize=50)    # fontsize of the tick labels
    ax.grid(True)

#############################################################
#### Prepare also the plots for the zoom-in on oscillations
#############################################################
### For oscillations on c
fig2 = plt.figure(figsize=(15,15))
spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig2)
ax2 = fig2.add_subplot(spec2[0, 0])
ax2.set_xlim([0.0,0.07])
ax2.set_ylim([-1e-4,5e-4])
ax2.tick_params(labelsize=32)    # fontsize of the tick labels
ax2.grid(True)
ax2.yaxis.set_major_formatter(formatter)
# change all spines
for axis in ['top', 'bottom', 'left', 'right']:
    ax2.spines[axis].set_linewidth(8)
# increase tick width
ax2.tick_params(width=8)

### For oscillations on phii at bottom bc
fig3 = plt.figure(figsize=(15,15))
spec3 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig3)
ax3 = fig3.add_subplot(spec2[0, 0])
ax3.set_xlim([0.0,0.08])
ax3.set_ylim([0.325,0.55])
ax3.tick_params(labelsize=32)    # fontsize of the tick labels
ax3.grid(True)
ax3.yaxis.set_major_formatter(formatter)
# change all spines
for axis in ['top', 'bottom', 'left', 'right']:
    ax3.spines[axis].set_linewidth(8)
# increase tick width
ax3.tick_params(width=8)
################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    path_root = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig7/Output_FEniCS/.')
    #### Before all, plot the initial fields in black on the three main Figs
    for (axes,mesh) in zip([axes1000, axes350, axes100], [1000, 350,100]):
        ### Initial field of deposition rate
        x0, y0 = getNodeArrays_Mixed(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_{}_dt_60'.format(mesh), 'data_mixed_000000.vtu').absolute()), DoF=0)  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        ax = axes[0]
        ax.plot(x0 * 100, 917 * 3770 * y0, color='k', linewidth=4, linestyle='-')

        ### Initial field of phii
        x00, y00 = getNodeArrays_Phii(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_{}_dt_60'.format(mesh),'data_phi_000000.vtu').absolute()))  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        ax = axes[1]
        ax.plot(x00 * 100, y00, color='k', linewidth=4, linestyle='-')

    #### Same thing for zoom-in fig but limited to finer grid
    x0, y0 = getNodeArrays_Mixed(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_1000_dt_60','data_mixed_000000.vtu').absolute()),DoF=0)  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
    ax2.plot(x0 * 100, 917 * 3770 * y0, color='k', linewidth=3, linestyle='-')
    x00, y00 = getNodeArrays_Phii(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_1000_dt_60','data_phi_000000.vtu').absolute()))  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
    ax3.plot(x00 * 100, y00, color='k', linewidth=4, linestyle='-')

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig8')
    ###List all cases to plot for Calonne
    Case_list = ['Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3','Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3','Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3','Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3','Fig7_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'Fig8_CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3']
    Settlement_list =['SettlementOff', 'SettlementOn','SettlementOff', 'SettlementOn','SettlementOff', 'SettlementOn', 'SettlementOff', 'SettlementOn','SettlementOff', 'SettlementOn','SettlementOff', 'SettlementOn']
    Mesh_list = ['1001', '1001', '1001', '1001', '351', '351', '351', '351','101', '101', '101', '101']
    axes_list = [axes1000, axes1000, axes1000, axes1000, axes350, axes350, axes350, axes350, axes100, axes100, axes100, axes100]
    Time_list_MyCode = [2880, 2880, 1440, 1440, 2880, 2880, 1440, 1440, 2880, 2880, 1440, 1440]
    MyColor_list_PerCase = ['orange', Blue, 'orange', Blue,'orange',  Blue, 'orange', Blue, 'orange', Blue, 'orange', Blue]
    MyLinesstyle_list = ['-','-',':',':','-','-',':',':','-','-',':',':']
    MyLinewidth_list = [5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
    ### For the Calonne model
    for k, (case, settl, mesh, axes, time) in enumerate(zip(Case_list, Settlement_list, Mesh_list, axes_list, Time_list_MyCode)):
        Depository_name = 'Simu_Papier_{}_{}_DepOn_{}nodes_2d_Output3h_Implicit_dt1min_BCDirichlet_'.format(case,settl,mesh)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name)
        file_name_nodal = '{}_NODAL_dt_{}.csv'.format(Depository_name, time)
        file_nodal = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        file_name_elemental = '{}_ELEMENTAL_dt_{}.csv'.format(Depository_name, time)
        file_elemental = pd.read_csv(full_path_mycode.joinpath(file_name_elemental))
        ###############################################################################
        # MAIN FIGURES THAT ARE FILLED WITH MY CODE RESULTS #####
        ################################################################################
        ax = axes[0]  # On Subplot0 we plot c profiles
        ax.plot(file_nodal['z']*100, file_nodal['DepositionRate'],color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])
        ax = axes[1]  # On Subplot1 we plot phii profiles
        ax.plot(0.5 * (file_elemental['zbot'] + file_elemental['ztop']) * 100, file_elemental['Phii'],color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])
        ### Plot zoomed-in for finer grid only
        if mesh == '1001':
            ax2.plot(file_nodal['z']*100, file_nodal['DepositionRate'],color=MyColor_list_PerCase[k], linewidth=4, linestyle=MyLinesstyle_list[k])
            ax3.plot(0.5*(file_elemental['zbot']+file_elemental['ztop'])*100, file_elemental['Phii'], color=MyColor_list_PerCase[k], linewidth=5,linestyle=MyLinesstyle_list[k])

    ################################################################################
    # DUMMY PLOTS FOR LEGEND #####
    ################################################################################
    ax = axes1000[0]
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='Init', color='k', linewidth=4, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='Settl. Off', color='orange', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='Settl. On', color=Blue, linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='t = 1d', color='darkgrey', linewidth=5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='t = 2d', color='darkgrey', linewidth=5, linestyle='-')
    ax.legend(loc='upper right', fontsize=31, ncol=2)

    ax = axes1000[1]
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label='Init', color='k', linewidth=4, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='Settl. Off', color='orange', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='Settl. On', color=Blue, linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label='t = 1d', color='darkgrey', linewidth=5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label='t = 2d', color='darkgrey', linewidth=5, linestyle='-')
    ax.legend(loc='upper right', fontsize=31, ncol=2)

    ################################################################################
    # SAVE THE FIGURE #####
    #######################
    ### For the paper we save the fig. for 1000 nodes only
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig8_Paper')
    name_output_fig_fig1000 = 'Fig_8_Papier_Main_1000nodes_WithLegend'
    fig1000.savefig(path_output_fig.joinpath(name_output_fig_fig1000))

    plt.show()