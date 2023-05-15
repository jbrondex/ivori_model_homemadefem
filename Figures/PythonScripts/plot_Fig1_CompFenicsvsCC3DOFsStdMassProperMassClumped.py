################################################################################
# import #####
################################################################################
import matplotlib.pyplot as plt

from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
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
    x = vtk_to_numpy(points.GetData())[:,0]

    usg = dsa.WrapDataObject( data )
    y = usg.PointData[reader.GetPointArrayName(0)]
    print('shape of x:',x.shape)
    print('shape of y:',y.shape)
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
    x = vtk_to_numpy(points.GetData())[:,0]

    usg = dsa.WrapDataObject( data )
    y = usg.PointData[reader.GetPointArrayName(0)][:,DoF] ###DoF is 0 for vn, 1 for Rhov, 2 for T
    print('shape of x:',x.shape)
    print('shape of y:',y.shape)
    print('... done.\n')


    return x, y

################################################################################
# For the Fenics Model, FUNCTION THAT PLOT FIGURE#####
################################################################################
def load_data_and_plot_figure(fig, axes ,path, ne, dt, file, leg, idx):
     
    ## configure the properties of plot attributes
    color_list = ['k', 'limegreen']
    alpha_list = [0.95,0.95]
    linesstyle_list = ['-', '-']
    linewidth_list = [2.8, 3]
    alpha=alpha_list[idx]
    color = color_list[idx]
    line = linesstyle_list[idx]
    width =  linewidth_list[idx]
    
    ## load data
    path_root = Path(path)
    print('File to load is=', str(path_root.joinpath(file).absolute()))
    # T
    x0,y0 = getNodeArrays_Mixed(str(path_root.joinpath(file).absolute()), DoF = 2) ###T is DoF 2 of solution vector mixed = (vn, Rhov,T)
    # rhov
    x1,y1 = getNodeArrays_Mixed(str(path_root.joinpath(file).absolute()), DoF = 1) ###rhov is DoF 1 of solution vector mixed = (vn, Rhov,T)
    # vn
    x2,y2 = getNodeArrays_Mixed(str(path_root.joinpath(file).absolute()), DoF = 0) ###vn is DoF 0 of solution vector mixed = (vn, Rhov,T)
        
    ## do the plots
    ax = axes[0]#On Subplot0 we plot Temperature
    ax.plot(x0,y0, label=leg, color=color, linewidth=width ,alpha=alpha, linestyle=line)
    ax.yaxis.set_major_formatter(formatter)
     
    ax = axes[1] #On Subplot1 we plot Rhov
    ax.plot(x1, y1, label=leg, color=color, linewidth=width ,alpha=alpha, linestyle=line) 
    ax.yaxis.set_major_formatter(formatter)
    
    ax = axes[2] #On Subplot1 we plot deposition rate, i.e. rhoi*s*vn for the Fenics Model
    ax.plot(x2, 917*3770*y2, label=leg, color=color, linewidth=width ,alpha=alpha, linestyle=line) 
    ax.yaxis.set_major_formatter(formatter)
    
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
   
ax=axes[1]
ax.set_ylabel(r'$\rho_\mathrm{v}$ ($\mathrm{kg~m^{-3}}$)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([0.0,0.005])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)

ax=axes[2]
ax.set_ylabel(r'$c$ ($\mathrm{kg~m^{-3}~s^{-1}}$)',fontsize=28)
ax.set_xlabel(r'$z$ (m)',fontsize=28)
ax.set_xlim([0.0,1.0])
ax.set_ylim([-7e-6,7e-6])
ax.tick_params(labelsize=24)    # fontsize of the tick labels
ax.grid(True)


################################################################################
# Make the plots #####
################################################################################
if __name__ == "__main__":    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~    FOR THE FENICS MODEL: LOAD ALL DATA VIA FUNCTION AND PLOT VIA OTHER FUNCTION     ~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    path_root = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig1/Output_FEniCS/.')
    config_list = []
    Depository_list = ['Scenario2_FullImplicit_NoBConVn_dt15mn_ttot38h_OneOutperHour', 'Scenario2_FullImplicit_NoBConVn_dt15mn_ttot38h_OneOutperHour']
    ne_list = [200, 200]
    dt_list = [900, 900]
    file_name_list = ['data_mixed_000000.vtu', 'data_mixed_000038.vtu']
    legend_list = ['Init', 'FEniCS']
    for (Depository_name, nedx,dtdx, file_name, leg) in zip(Depository_list, ne_list,dt_list, file_name_list, legend_list):
        config = {
            'path': path_root.joinpath(Depository_name,'data','ne_{}_dt_{}'.format(nedx,dtdx)),
            'ne': nedx,
            'dt': dtdx,
            'file' : file_name,
            'leg': leg
        }
        config_list.append(config)
        
    #### For the Fenics Model, load data and plot figure
    for idx,config in enumerate(config_list):
        print('Index =', idx)
        fig, axes = load_data_and_plot_figure(fig=fig, axes=axes ,path=config['path'], ne=config['ne'], dt=config['dt'], file=config['file'], leg=config['leg'], idx=idx)
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~    FOR MY MODEL: LOAD ALL DATA AND PLOT DIRECTLY IN THIS PART OF THE CODE   ~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
    pathroot_mycode = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Output/Output_Fig1/.')
    
    ###List all cases to plot
    Depository_name_list = ['Simu_Papier_Fig1_CC_3DOFs_StdMass_cnotlumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_', 'Simu_Papier_Fig1_CC_3DOFs_ProperMass_cnotlumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_', 'Simu_Papier_Fig1_CC_3DOFs_ProperMass_cLumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_']
    # Depository_name_list = ['Simu_Papier_Fig1_CC_3DOFs_StdMass_cLumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_', 'Simu_Papier_Fig1_CC_3DOFs_ProperMass_cnotlumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_', 'Simu_Papier_Fig1_CC_3DOFs_ProperMass_cLumped_CalPar_Alpha5minus3_SettlementOff_201nodes_38h_Implicit_dt15min_']
    dt_list_MyCode = [152, 152, 152]
    color_list_MyCode = ['darkred', Blue, 'darkorange'] ###legend for output of my code 
    Mylinesstyle_list = [':', '--', '-']
    Mylinewidth_list = [3, 3.1, 3]    
    Mylegend_list = ['CC_3DOF; No lump.; Improp.', 'CC_3DOF; No Lump.; Prop.','CC_3DOF; Lump.; Prop.'] ###legend for output of my code
    
    ###Load outputs corresponding to considered cases 
    dict_Mycode_Nodal_list = [] ##Liste des données csv chargées (une par cas envisagé)
    for (Depository_name, dt) in zip(Depository_name_list, dt_list_MyCode):        
        file_name_nodal =  '{}_NODAL_dt_{}.csv'.format(Depository_name,dt)
        full_path_mycode = pathroot_mycode.joinpath(Depository_name) 
        file = pd.read_csv(full_path_mycode.joinpath(file_name_nodal))
        dict_Mycode_Nodal_list.append(file)
           
#    ###Plot the considered cases
    ax = axes[0]#On Subplot0 we plot Temperature
    ax.plot(dict_Mycode_Nodal_list[0]['z'], dict_Mycode_Nodal_list[0]['Temperature'], label= Mylegend_list[0], color=color_list_MyCode[0],linewidth=  Mylinewidth_list[0] , linestyle= Mylinesstyle_list[0]) #, marker=marker,markersize=markersize)
    ax.plot(dict_Mycode_Nodal_list[1]['z'], dict_Mycode_Nodal_list[1]['Temperature'], label=Mylegend_list[1], color=color_list_MyCode[1],linewidth= Mylinewidth_list[1] , linestyle=Mylinesstyle_list[1])
    ax.plot(dict_Mycode_Nodal_list[2]['z'], dict_Mycode_Nodal_list[2]['Temperature'], label= Mylegend_list[2], color=color_list_MyCode[2],linewidth=  Mylinewidth_list[2] , linestyle= Mylinesstyle_list[2]) #, marker=marker,markersize=markersize)
    ax.legend(loc='upper right', fontsize=18, ncol=3)

    ax = axes[1] #On Subplot1 we plot Rhov
    ax.plot(dict_Mycode_Nodal_list[0]['z'], dict_Mycode_Nodal_list[0]['Rhov'], label= Mylegend_list[0], color=color_list_MyCode[0],linewidth=  Mylinewidth_list[0] , linestyle= Mylinesstyle_list[0]) #, marker=marker,markersize=markersize)
    ax.plot(dict_Mycode_Nodal_list[1]['z'], dict_Mycode_Nodal_list[1]['Rhov'], label=Mylegend_list[1], color=color_list_MyCode[1],linewidth= Mylinewidth_list[1] , linestyle=Mylinesstyle_list[1])
    ax.plot(dict_Mycode_Nodal_list[2]['z'], dict_Mycode_Nodal_list[2]['Rhov'], label= Mylegend_list[2], color=color_list_MyCode[2],linewidth=  Mylinewidth_list[2] , linestyle= Mylinesstyle_list[2]) #, marker=marker,markersize=markersize)
    ax.legend(loc='upper right', fontsize=18, ncol=3)

    ax = axes[2] #On Subplot1 we plot deposition rate, i.e. rhoi*s*vn for the Fenics Model
    ax.plot(dict_Mycode_Nodal_list[0]['z'], dict_Mycode_Nodal_list[0]['DepositionRate'], label= Mylegend_list[0], color=color_list_MyCode[0],linewidth=  Mylinewidth_list[0] , linestyle= Mylinesstyle_list[0]) #, marker=marker,markersize=markersize)
    ax.plot(dict_Mycode_Nodal_list[1]['z'], dict_Mycode_Nodal_list[1]['DepositionRate'], label=Mylegend_list[1], color=color_list_MyCode[1],linewidth= Mylinewidth_list[1] , linestyle=Mylinesstyle_list[1])
    ax.plot(dict_Mycode_Nodal_list[2]['z'], dict_Mycode_Nodal_list[2]['DepositionRate'], label= Mylegend_list[2], color=color_list_MyCode[2],linewidth=  Mylinewidth_list[2] , linestyle= Mylinesstyle_list[2]) #, marker=marker,markersize=markersize)
    ax.legend(loc='upper left', fontsize=18, ncol=3)
    ################################################################################
    # SAVE THE FIGURES #####
    ################################################################################
    name_output_fig = 'Fig_1_Papier_WithLegend'
    path_output_fig = Path('/home/brondexj/Documents/IVORI/ivori_model_homemadefem/Figures/Fig1_Paper/.')
    fig.savefig(path_output_fig.joinpath(name_output_fig))
    plt.show()