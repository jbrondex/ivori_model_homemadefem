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


######################################################
# For the Fenics Model, FUNCTION THAT PLOT FIGURE#####
######################################################
def load_data_and_plot_figure(fig, axes, path, ne, dt, file, col, style, width, idx):
    ## load data
    path_root = Path(path)
    print('File to load is=', str(path_root.joinpath(file).absolute()))
    if 'mixed' in file:
        # vn
        x1, y1 = getNodeArrays_Mixed(str(path_root.joinpath(file).absolute()), DoF=0)  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        ## do the plots
        ax = axes[0, 0]  # On Subplot0 we plot deposition rate
        ax.plot(x1*100,  917 * 3770 * y1,  color=col, linewidth=width, linestyle=style)
        ax.yaxis.set_major_formatter(formatter)
        ax = axes[0, 1]  # On Subplot2 we plot zoom-in on deposition rate
        ax.plot(x1*100,  917 * 3770 * y1,  color=col, linewidth=width, linestyle=style)
        ax.yaxis.set_major_formatter(formatter)

    elif 'phi' in file:
        # phii
        x2, y2 = getNodeArrays_Phii(str(path_root.joinpath(file).absolute()))
        ## do the plots
        ax = axes[1, 0]  # On Subplot1 we plot phii
        ax.plot(x2*100, y2, color=col, linewidth=width, linestyle=style)
        ax = axes[1, 1]  # On Subplot3 we plot zoom in phii
        ax.plot(x2*100, y2, color=col, linewidth=width, linestyle=style)


    return fig, axes


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

### WE DO THREE TIMES THE SAME PLOT FOR THE THREE CONSIDERED RESOLUTION (1000, 350 and 100 nodes)
### Set-Up a list of the three fig
fig1000, axes1000 = plt.subplots(2, 2, figsize=(28, 20))
fig350, axes350 = plt.subplots(2, 2, figsize=(28, 20))
fig100, axes100 = plt.subplots(2, 2, figsize=(28, 20))

###pimp style for the three fig
for axes in [axes1000, axes350, axes100]:
    # pimp style
    ###Deposition rate for Calonne System
    ax=axes[0, 0]
    ax.set_ylabel(r'$c$ ($\mathrm{kg~m^{3}~s^{-1}}$)',fontsize=44)
    ax.set_xlim([0.0,2])
    ax.set_ylim([-2e-2,2e-2])
    ax.tick_params(labelsize=36)    # fontsize of the tick labels
    ax.grid(True)
    ###Phii for Calonne System
    ax=axes[1, 0]
    ax.set_ylabel(r'$\Phi_\mathrm{i}$ (-)',fontsize=44)
    ax.set_xlabel(r'$z$ (cm)',fontsize=44)
    ax.set_xlim([0.0,2])
    ax.set_ylim([0.19,0.550])
    ax.tick_params(labelsize=36)    # fontsize of the tick labels
    ax.grid(True)
    ###Deposition rate Zoom-in
    ax=axes[0, 1]
    ax.set_xlim([0.0,0.1])
    ax.set_ylim([-5e-5,4.5e-4])
    ax.tick_params(labelsize=36)    # fontsize of the tick labels
    ax.grid(True)
    ax.yaxis.set_major_formatter(formatter)
    ###Phii Zoom-in
    ax=axes[1 ,1]
    ax.set_xlabel(r'$z$ (cm)',fontsize=44)
    ax.set_xlim([0.0,0.1])
    ax.set_ylim([0.33,0.371])
    ax.tick_params(labelsize=36)    # fontsize of the tick labels
    ax.grid(True)

################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~    FOR THE FENICS MODEL: LOAD ALL DATA VIA FUNCTION AND PLOT VIA OTHER FUNCTION     ~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    path_root = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig7/Output_FEniCS/.')
    config_list = []
    Depository_list = ['MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate','MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate','MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate']
    ne_list = [1000, 1000, 350, 350, 100, 100]
    dt_list = [60, 60, 60, 60, 60, 60]
    file_name_list = ['data_mixed_000288.vtu', 'data_phi_000288.vtu', 'data_mixed_000288.vtu', 'data_phi_000288.vtu', 'data_mixed_000288.vtu', 'data_phi_000288.vtu']
    color_list = ['limegreen', 'limegreen', 'limegreen', 'limegreen', 'limegreen', 'limegreen']
    linesstyle_list = ['-', '-', '-', '-', '-', '-']
    linewidth_list = [4, 4, 4, 4, 4, 4]
    for (Depository_name, nedx, dtdx, file_name, color, style, width) in zip(Depository_list, ne_list, dt_list, file_name_list, color_list, linesstyle_list, linewidth_list):
        config = {
            'path': path_root.joinpath(Depository_name, 'ne_{}_dt_{}'.format(nedx, dtdx)),
            'ne': nedx,
            'dt': dtdx,
            'file': file_name,
            'col' : color,
            'style' : style,
            'width' : width
        }
        config_list.append(config)
    #### Before all, plot the initial fields in black on the three main Figs
    for (axes,mesh) in zip([axes1000, axes350, axes100], [1000, 350,100]):
        ### Initial field of deposition rate
        x0, y0 = getNodeArrays_Mixed(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_{}_dt_60'.format(mesh), 'data_mixed_000000.vtu').absolute()), DoF=0)  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        ax = axes[0, 0]
        ax.plot(x0 * 100, 917 * 3770 * y0, color='k', linewidth=3, linestyle='-')
        ax = axes[0, 1]
        ax.plot(x0 * 100, 917 * 3770 * y0, color='k', linewidth=3, linestyle='-')
        ### Initial field of phii
        x00, y00 = getNodeArrays_Phii(str(path_root.joinpath('MySimu_Scenario3_Implicit_Papier_Fig7_NoBCDepRate', 'ne_{}_dt_60'.format(mesh),'data_phi_000000.vtu').absolute()))  ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        ax = axes[1, 0]
        ax.plot(x00 * 100, y00, color='k', linewidth=3, linestyle='-')
        ax = axes[1, 1]
        ax.plot(x00 * 100, y00, color='k', linewidth=3, linestyle='-')

    #### For the Fenics Model, load data and plot figure
    for idx, config in enumerate(config_list):
        print('Index =', idx)
        if config['ne']== 1000:
            fig, axes = fig1000, axes1000
        elif config['ne']== 350:
            fig, axes = fig350, axes350
        elif config['ne']== 100:
            fig, axes = fig100, axes100
        fig, axes = load_data_and_plot_figure(fig=fig, axes=axes, path=config['path'], ne=config['ne'], dt=config['dt'], file=config['file'], col=config['col'], style=config['style'], width=config['width'], idx=idx)

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig7/.')
    ###List all cases to plot for Calonne
    Case_list = ['CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'CC_3DOFs_StdMass_cnotlumped_CalPar_Alpha5minus3', 'CC_3DOFs_cLumped_MassnotLumped_ConstPar_Alpha5minus3', 'CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'CC_3DOFs_StdMass_cnotlumped_CalPar_Alpha5minus3', 'CC_3DOFs_cLumped_MassnotLumped_ConstPar_Alpha5minus3','CC_3DOFs_cLumped_MassnotLumped_CalPar_Alpha5minus3', 'CC_3DOFs_StdMass_cnotlumped_CalPar_Alpha5minus3', 'CC_3DOFs_cLumped_MassnotLumped_ConstPar_Alpha5minus3']
    Mesh_list = ['1001', '1001', '1001', '351', '351', '351', '101', '101', '101']
    axes_list = [axes1000, axes1000,  axes1000, axes350, axes350, axes350, axes100, axes100, axes100]
    Time_list_MyCode = [2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880]
    MyColor_list_PerCase = [ 'orange', 'darkred', 'r', 'orange', 'darkred', 'r',  'orange', 'darkred', 'r']#, 'orchid']#, (42/255, 127/255, 255/255), 'limegreen']
    MyLinesstyle_list = ['-',':','-','-',':','-','-',':','-']
    MyLinewidth_list = [5, 5, 4, 5, 5, 4, 5, 5, 4] #,1.3]
    ### For the Calonne model
    for k, (case,mesh, axes, time) in enumerate(zip(Case_list,Mesh_list, axes_list, Time_list_MyCode)):
        Depository_name = 'Simu_Papier_Fig7_{}_SettlementOff_DepOn_{}nodes_2d_Output3h_Implicit_dt1min_BCDirichlet_'.format(case,mesh)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name)
        file_name_nodal = '{}_NODAL_dt_{}.csv'.format(Depository_name, time)
        file_nodal = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        file_name_elemental = '{}_ELEMENTAL_dt_{}.csv'.format(Depository_name, time)
        file_elemental = pd.read_csv(full_path_mycode.joinpath(file_name_elemental))
        ###############################################################################
        # MAIN FIGURES THAT ARE FILLED WITH MY CODE RESULTS #####
        ################################################################################
        ax = axes[0,0]  # On Subplot0 we plot c profiles
        ax.plot(file_nodal['z']*100, file_nodal['DepositionRate'],color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])
        ax = axes[1, 0]  # On Subplot1 we plot phii profiles
        ax.plot(0.5 * (file_elemental['zbot'] + file_elemental['ztop']) * 100, file_elemental['Phii'],color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])
        ax = axes[0, 1]  # On Subplot2 we plot zoom in on c
        ax.plot(file_nodal['z'] * 100, file_nodal['DepositionRate'], color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])
        ax = axes[1, 1]  # On Subplot3 we plot zoom in on phii
        ax.plot(0.5 * (file_elemental['zbot'] + file_elemental['ztop']) * 100, file_elemental['Phii'], color=MyColor_list_PerCase[k], linewidth=MyLinewidth_list[k], linestyle=MyLinesstyle_list[k])

    ################################################################################
    # DUMMY PLOTS FOR LEGEND #####
    ################################################################################
    ax =axes1000[0,0]
    #### Dummy plot for legend
    ax.plot(np.NaN, np.NaN, label = 'Init', color='k', linewidth=3, linestyle='-')
    ax.plot(np.NaN, np.NaN, label = 'FEniCS', color='limegreen', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label = 'CC_3DOF: No lump. / Improp.', color='darkred', linewidth=5, linestyle=':')
    ax.plot(np.NaN, np.NaN, label = 'CC_3DOF: Lump. / Prop.', color='orange', linewidth=5, linestyle='-')
    ax.plot(np.NaN, np.NaN, label = 'CC_3DOF: $D^\mathrm{eff}$ & $k^\mathrm{eff}$ cst.', color='r', linewidth=5, linestyle='-')

    lines = []
    labels = []
    for ax in fig1000.axes:
        Line, Label = ax.get_legend_handles_labels()
        # print(Label)
        lines.extend(Line)
        labels.extend(Label)
    fig1000.legend(lines, labels, loc='upper center', fontsize=34, ncol=3)

    ################################################################################
    # SAVE THE FIGURE #####
    #######################
    ### For the paper we save the fig. for 1000 nodes only
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig7_Paper')
    name_output_fig1000 = 'Fig_7_Papier_Main_1000nodes_WithLegend'
    fig1000.savefig(path_output_fig.joinpath(name_output_fig1000))


    plt.show()