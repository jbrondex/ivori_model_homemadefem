import copy
import csv  ## For the output
import importlib  ### For choosing which sif must be run
from pathlib import Path  ## To manage the file to produce outputs

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import FEM_classes
import USF
import MatLaw
import Solvers.CoupledCalonneHeatVap_Solver_3DOFs_cnotlumped
import Solvers.CoupledCalonneHeatVap_Solver_3DOFs_cnotlumped_StdMassMatrix
import Solvers.CoupledCalonneHeatVap_Solver_3DOFs_cLumped_StdMassMatrix
import Solvers.CoupledCalonneHeatVap_Solver_2DOFs_KreactLumped
import Solvers.CoupledCalonneHeatVap_Solver_3DOFs_cLumped
import Solvers.CalonneHeat_Solver_1DOFs_KreactLumped
import Solvers.CalonneVap_Solver_1DOFs_KreactLumped
import Solvers.CalonneVap_Solver_1DOFs_KreactLumped_SourceExplicit
import Solvers.HansenHeatVap_Solver_1DOF_TForm
import Solvers.HansenHeatVap_Solver_2DOFs_MixedForm
import Solvers.DepositionRateHansen_Solver
import Solvers.EnergyConservation_Solver
import Solvers.Stress_Solver
import Solvers.Settlement_Solver
### import sif_Simulations_Paper.sif_Fig2_CC3DOFsvsCDPC.sif_scenario2_schurholt_CC_3DOFs as sif
####/////////////////////////////####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####   CHOOSE WHICH SIF TO RUN   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
Directory_sif = 'sif_Simulations_Paper.sif_Fig2_CC3DOFsvsCDPC' ### THAT ARE THE ONLY LINES OF THE main.py THAT THE STANDARD USER MUST ADAPT
Name_sif = 'sif_scenario2_schurholt_CC_3DOFs' ### THAT ARE THE ONLY LINES OF THE main.py THAT THE STANDARD USER MUST ADAPT
sif = importlib.import_module(Directory_sif+'.'+Name_sif, package=None)

####/////////////////////####
####~~~~~~~~~~~~~~~~~~~~~####
####   INITIALIZATION    ####
####~~~~~~~~~~~~~~~~~~~~~####
####\\\\\\\\\\\\\\\\\\\\\####

###~~~~~~~~~~~~~~~~~~~~~~~~~~###
###  Initialize Simulation   ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~###
Simulation = FEM_classes.model_t(sif.Simu, sif.Snowpack, sif.Cst, sif.Mat, sif.BF, sif.BC1, sif.Solver1) ##Creates attribute Simulation of class model_t and all inner classes, creates mesh and create fields always present

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###  Initialize Variable Fields from Initial Conditions of sif  ###
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
for field in vars(sif.IC): ###loop on all fields defined in IC section of sif
    if field not in Simulation.snowpack.field_perm:
        raise TypeError('Field '+field + ' prescribed in Initial Conditions of sif does not exist in the model.')
    #if constant value
    elif type(eval('sif.IC.'+field)) == float or type(eval('sif.IC.'+field)) == int:
        Simulation.snowpack.field_perm[field].value[:] = eval('sif.IC.'+field)
        #if function
    elif type(eval('sif.IC.'+field)) == str: ### Initial Field prescribe through USF or Material Law
        Name_func = eval('sif.IC.'+field)
        liste = Name_func.split('_')
        if liste[0] not in ['USF', 'MatLaw']:
            raise TypeError('If not constant, function prescribing Initial '+ field +' is a string that must start by "USF" or "MatLaw"')
        if np.size(liste) < 3:
            raise TypeError('If not constant, function prescribing Initial '+ field +' is a string with formalism "USF/MatLaw_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
        input_fields = liste[2:] ### Get the list of variables that are the input of the function
        if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'Phii']): ### accepted are only in this list
            raise TypeError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')
        # if the field to initialize is defined at nodes
        if Simulation.snowpack.field_perm[field].NodEl == 'Nodal': 
            for i in range(Simulation.snowpack.mesh.numberofNodes):
                input_values = []
                for input_field in input_fields:
                    if Simulation.snowpack.field_perm[input_field].NodEl =='Nodal': ### If input field of function is nodal, we take the value at node
                        value = Simulation.snowpack.field_perm[input_field].value[i]
                    elif Simulation.snowpack.field_perm[input_field].NodEl =='Elemental':### If input field of function is elemental field...
                        value = 0.5*(Simulation.snowpack.field_perm[input_field].value[i-1] + Simulation.snowpack.field_perm[input_field].value[i]) ###...value at node i is arithmetic mean between values at surrounding elements
                    input_values.append(value)
                Simulation.snowpack.field_perm[field].value[i] = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
            # if the field to initialize is defined at elements
        if Simulation.snowpack.field_perm[field].NodEl == 'Elemental': 
            for i in range(Simulation.snowpack.mesh.numberofElements):
                input_values = []
                for input_field in input_fields:
                    if Simulation.snowpack.field_perm[input_field].NodEl =='Nodal': ### If input field of function is nodal...                        
                        value_below = Simulation.snowpack.field_perm[input_field].value[i]  ###...we get value at node below
                        value_above = Simulation.snowpack.field_perm[input_field].value[i+1] ###...and value at node above...
                        value = (value_below + value_above)/2 ### ... to calculate the averaged value of the input field over the element
                    elif Simulation.snowpack.field_perm[input_field].NodEl =='Elemental':### If input field of function is elemental field...
                        value = Simulation.snowpack.field_perm[input_field].value[i] ### ...we assign value of element i 
                    input_values.append(value)
                Simulation.snowpack.field_perm[field].value[i] = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values) 
    Simulation.snowpack.field_perm[field].value_prev_it[:]= Simulation.snowpack.field_perm[field].value[:] ### BE CARREFUL : The value of field at previous iteration must also be initialised     
    Simulation.snowpack.field_perm[field].value_prev_tsp[:]= Simulation.snowpack.field_perm[field].value[:] ### BE CARREFUL : The value of field at previous tsp must also be initialised (required for e.g. diagnostic of deposition rate in Hansen Solver)  
    
#### Check up for mass conservation
Mass_element_prev_tsp=[]
for i in range(Simulation.snowpack.mesh.numberofElements):
    Mass_element_prev_tsp.append(Simulation.snowpack.field_perm['Phii'].value[i]*917*Simulation.snowpack.mesh.elements[i].get_elsize())
        
Total_Mass = np.sum(Mass_element_prev_tsp)
print('Total Mass init=', Total_Mass)

#### Check up for energy conservation
Energy_Leak = 0

####/////////////////////####
####~~~~~~~~~~~~~~~~~~~~~####
####      TIME LOOP      ####
####~~~~~~~~~~~~~~~~~~~~~####
####\\\\\\\\\\\\\\\\\\\\\####
if Simulation.type == 'Transient':
    #######################################
    ###~~~~~~~ Loop on timesteps ~~~~~~~###
    #######################################
    while  Simulation.timestep <= Simulation.totaltimesteps:
        print('\n|#######################|')
        print('|~~~~~~~~~~~~~~~~~~~~~~~|')
        print('|     Tsp =',Simulation.timestep,'/',Simulation.totaltimesteps,'   |')
        print('|~~~~~~~~~~~~~~~~~~~~~~~|')
        print('|#######################|')
        #####################################
        ###~~~~~~~ Loop on Solvers ~~~~~~~###
        #####################################
        for Solver in Simulation.solver_list:   
            
            if Simulation.timestep % Solver.execsolver == 0: ## Check whether solver must be executed at considered tsp (i.e tsp can be divided by ExecSolver)...
                print('\nSolver=', Solver.name)
                eval(Solver.file+'.'+Solver.name)(Solver) ###... if yes, execute solver

                ### Mass conservation check-up
                if Solver.name == 'Settlement_Solver':
                    print('\n~~~~~~~~~~~~~~~~~~~~~~~')
                    print('~~~MASS CONSERVATION~~~')
                    print('~~~~~~~~~~~~~~~~~~~~~~~')
                    Mass_element=[]
                    for i in range(Simulation.snowpack.mesh.numberofElements):
                        Mass_element.append(Simulation.snowpack.field_perm['Phii'].value[i]*917*Simulation.snowpack.mesh.elements[i].get_elsize())
                
                    Delta_Mass_element=[]
                    for i in range(Simulation.snowpack.mesh.numberofElements):
                        Delta_Mass_element.append(Mass_element[i] - Mass_element_prev_tsp[i])
            
                    print('Total Mass prev tsp=', np.sum(Mass_element_prev_tsp))                   
                    print('Total Mass=', np.sum(Mass_element))
                    print('Delta Total Mass From Prev It=', np.sum(Mass_element) - np.sum(Mass_element_prev_tsp))
                    ### Update mass of previous timestep to current for next visit
                    Mass_element_prev_tsp =  copy.deepcopy(Mass_element)            
                
        ###~~~~~~~ End loop on Solvers ~~~~~~~###
        
        ### Energy conservation check-up
        if Simulation.timestep == 1 and Simulation.save_csv:
            ## Store energy information in dedicated file for energy budget (paper) ---> At first tsp prepare the file only
            #setup directory
            name= Simulation.name
            path= Simulation.outputpath
            path_root = Path(path).joinpath(name)
            path_root.mkdir(exist_ok=True,parents=True)
            #setup file
            fname_Energy = str(path_root.joinpath(str(name+'_ENERGY.csv')).absolute())
            header_Energy = [str('time'), str('Energy Leak from beginning')]
            with open(fname_Energy, 'w', encoding='UTF8', newline='') as f_energy:
                writer=csv.writer(f_energy)
                writer.writerow(header_Energy) ##we write the header only once
        elif Simulation.timestep > 1: ### From timestep 2, fill-up the file and print relevant information
            print('\n~~~~~~~~~~~~~~~~~~~~~~~~~')
            print('~~~ENERGY CONSERVATION~~~')
            print('~~~~~~~~~~~~~~~~~~~~~~~~~')

            Total_Energy = np.sum(Simulation.snowpack.field_perm['Energy'].value[:])          
            print('Total Energy in System at current tsp =', Total_Energy)
            Total_Energy_prev_tsp = np.sum(Simulation.snowpack.field_perm['Energy'].value_prev_tsp[:])
            Delta_Energy = Total_Energy - Total_Energy_prev_tsp        
            print('Delta Energy from prev tsp =', Delta_Energy)  
            Energy_Flux_In_Sensible = (Simulation.snowpack.field_perm['Temperature'].residual[0] + Simulation.snowpack.field_perm['Temperature'].residual[-1]) * Simulation.timestepsize * Simulation.solver_list[0].execsolver
            Energy_Flux_In_Latent = Simulation.snowpack.constant_perm['Lm'].value*(Simulation.snowpack.field_perm['Rhov'].residual[0] + Simulation.snowpack.field_perm['Rhov'].residual[-1]) * Simulation.timestepsize * Simulation.solver_list[0].execsolver
            Energy_Flux_In = Energy_Flux_In_Sensible + Energy_Flux_In_Latent
            print('Energy Flux In =', Energy_Flux_In)
            print('Latent Energy Gone Settling =', np.sum(Simulation.snowpack.field_perm['EnergyGoneSettling'].value[:]))
            print('Energy Gap  =', Delta_Energy - Energy_Flux_In + np.sum(Simulation.snowpack.field_perm['EnergyGoneSettling'].value[:]))
            Energy_Leak += Delta_Energy - Energy_Flux_In + np.sum(Simulation.snowpack.field_perm['EnergyGoneSettling'].value[:])
            print('Total Energy Leak from beginning of Simu  =', Energy_Leak)
            
            ### Fill-up the Energy file with time and energy leak from beginning
            if Simulation.save_csv:
                with open(fname_Energy, 'a', encoding='UTF8', newline='') as f_energy:
                    writer=csv.writer(f_energy)
                    writer.writerow([Simulation.time, Energy_Leak])

        ##########################################
        ###~~~~~~~ Produce the outputs ~~~~~~~####
        ##########################################
        
        
        c = 3770*5e-3*MatLaw.MatLaw_vkin(Simulation.snowpack.field_perm['Temperature'].value) *(Simulation.snowpack.field_perm['Rhov'].value - MatLaw.MatLaw_RhovSatLibbrecht(Simulation.snowpack.field_perm['Temperature'].value))
                   
        if Simulation.timestep % Simulation.outputinterval == 0: ## Check whether an output must be produced at current timestep...
            ######################################
            ###   Production of output files   ###
            ######################################
           
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~##
            ## Produce the csv outputs ### 
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~##
            if Simulation.save_csv:  ### Only if explicitly demanded in the sif file
                #setup directory
                name= Simulation.name
                path= Simulation.outputpath
                path_root = Path(path).joinpath(name)
                path_root.mkdir(exist_ok=True,parents=True)
                #setup files
                ###---> We produce 3 csv files for each desired output time : one is general information about the simu, one is nodal fields, one is elemental fields
                fname_Info = str(path_root.joinpath(str(name+'_INFO_dt_'+str(Simulation.timestep)+'.csv')).absolute())   
                fname_Nodal = str(path_root.joinpath(str(name+'_NODAL_dt_'+str(Simulation.timestep)+'.csv')).absolute()) 
                fname_Elemental = str(path_root.joinpath(str(name+'_ELEMENTAL_dt_'+str(Simulation.timestep)+'.csv')).absolute())
                ### the fields that are always present in each of the three .csv file
                header_info = [str('Timestep n°='+str(Simulation.timestep)), str('dt ='+str(Simulation.timestepsize)+' s'), str('time ='+str(Simulation.time)+' s'), str('TimeStepping Method='+Simulation.timestepping_method), str('theta='+str(Simulation.theta))]
                header_nodal = [str('Node n°'), str('z')]
                header_elemental = [str('Elt n°'), str('zbot'), str('ztop')]
                ### Loop on nodes to treat nodal fields only
                with open(fname_Nodal, 'w', encoding='UTF8', newline='') as f_node:
                    writer=csv.writer(f_node)
                    for i,node in enumerate(Simulation.snowpack.mesh.nodes):
                        list_row=[node.numberinMesh, node.pos] 
                        for field in Simulation.snowpack.field_list:
                            field = Simulation.snowpack.field_perm[field.name] ### This line is ugly but is used to go around a bug that needs to be fix later (the .value return 0.0 if accessed from field_list and true value if accessed from field_perm)
                            if field.NodEl=='Elemental' or field.name =='Coordinate':
                                continue ### we consider only nodal field for the nodal.csv file + we don't consider the field coordinate
                            if i==0: ##we fill up the header of the nodal.csv and nodal.info files only once
                                header_nodal.append(field.name) 
                                header_info.append(str(field.name+'('+field.NodEl+')'))                            
                            list_row.append(field.value[node.numberinMesh-1]) ###Fill up a given row with all value of all nodal fields at corresponding node
                        if i==0:
                            writer.writerow(header_nodal) ##we write the header only once
                        writer.writerow(list_row)
                    
                ### Loop on elements to treat elemental fields only
                with open(fname_Elemental, 'w', encoding='UTF8', newline='') as f_el:
                    writer=csv.writer(f_el)
                    for i,el in enumerate(Simulation.snowpack.mesh.elements):
                        list_row=[el.numberinMesh, el.nodes[0].pos, el.nodes[1].pos] 
                        for field in Simulation.snowpack.field_list:
                            field = Simulation.snowpack.field_perm[field.name] ### This line is ugly but is used to go around a bug that needs to be fix later (the .value return 0.0 if accessed from field_list and true value if accessed from field_perm)
                            if field.NodEl=='Nodal' or field.name =='Coordinate':
                                continue ### we consider only nelemental field for the elemental.csv file + we don't consider the field coordinate
                            if i==0: ##we fill up the header of the elemental.csv and nodal.info files only once
                                header_elemental.append(field.name) 
                                header_info.append(str(field.name+'('+field.NodEl+')'))                            
                            list_row.append(field.value[el.numberinMesh-1]) ###Fill up a given row with all value of all nodal fields at corresponding node
                        if i==0:
                            writer.writerow(header_elemental) ##we write the header only once
                        writer.writerow(list_row)                       
                     
                ### Write the info.csv info file that contains only the header
                with open(fname_Info, 'w', encoding='UTF8', newline='') as f_info: 
                    writer=csv.writer(f_info)
                    writer.writerow(header_info)
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
            ## End of output files production ### 
            ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
        ### If we are in a timestep which does not correspond to an output timestep... 
        else:
            pass ##...then don't do anything
        ###~~~~~~~ End of produce the outputs ~~~~~~~###





        ###~~~~~~~ Update timestep and time ~~~~~~~###
        Simulation.update_timestep()
        Simulation.update_time()
    

