import FEM_classes
import numpy as np
import MatLaw
import copy

def CoupledCalonneHeatVap(Solver):
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING GLOBAL MATRICES FROM LOCAL MATRICES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleBulkMatrix(Solver):
        Solver.cleanMatrices()
        # Perform Bulk Assemble
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()            
            Kloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Mloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Jloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Floc = np.zeros((Solver.DOFs*el.numberofnodes))
            for IP in range(el.numberofIP):
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)
                
                ### Get all material properties necessary to calculate Stiffness and force matrices at IP
                Deff = Solver.parent.snowpack.material_perm['Deff'].get_material(IPpos)
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
                s = Solver.parent.snowpack.bodyforce_perm['s'].get_BF(IPpos)
                alpha = Solver.parent.snowpack.bodyforce_perm['alpha'].get_BF(IPpos)
                vkin = Solver.parent.snowpack.bodyforce_perm['vkin'].get_BF(IPpos) 
                SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(IPpos) 
                dSWDdT = Solver.parent.snowpack.material_perm['dSWDdT'].get_material(IPpos)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(IPpos)
                Keff = Solver.parent.snowpack.material_perm['Keff'].get_material(IPpos)
                Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition
                T_prev= Solver.parent.snowpack.field_perm['Temperature'].get_value_prev_it(IPpos) ###Value of temperature obtained from previous iteration

                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )
                        
                        ###~~~~~~LOCAL MASS MATRIX~~~~~~###
                        Mloc[Solver.DOFs*p,Solver.DOFs*q] += basis[p] * basis[q] * IP_weight * el_size  #For DOFs 1, here Rhov
                        Mloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += basis[p] * basis[q] * IP_weight * el_size #For DOFs 2, here T
                        Mloc[Solver.DOFs*p+2,Solver.DOFs*q+2] += (Solver.execsolver*Solver.parent.timestepsize) * basis[p] * basis[q] * IP_weight * el_size ##For DOFs 3, here c
                        ### No time derivative on DOFs 3 (here c)
                        
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                        Kloc[Solver.DOFs*p,Solver.DOFs*q] += (Deff/(1-Phii)) * dbasis[q] * dbasis[p] *  IP_weight * el_size ##K_{PP}
                        Kloc[Solver.DOFs*p,Solver.DOFs*q+2] += (1/(1-Phii)) * basis[q] * basis[p] *  IP_weight * el_size ##K_{PC}                        
                        Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += (Keff/rhoCeff) * dbasis[q] * dbasis[p] *  IP_weight * el_size  ##K_{TT}
                        Kloc[Solver.DOFs*p+1,Solver.DOFs*q+2] += - (Lm/rhoCeff) * basis[q] * basis[p] *  IP_weight * el_size  ##K_{TC}
                        Kloc[Solver.DOFs*p+2,Solver.DOFs*q] += - (s*alpha*vkin) * basis[p] * basis[q] *  IP_weight * el_size  ##K_{CP}
                        Kloc[Solver.DOFs*p+2,Solver.DOFs*q+1] += (s*alpha*vkin)* dSWDdT * basis[p] * basis[q]*  IP_weight * el_size ##K_{CT}
                        ### With this formulation K_{PT} and K_{TP} are zeros, which is the default value of Kloc
                        
                        ###~~~~~~LOCAL JACOBIAN MATRIX~~~~~~### To implement later if Newton linearization needed
                        Jloc[Solver.DOFs*p,Solver.DOFs*q] += 0 *  IP_weight * el_size 
                        Jloc[Solver.DOFs*p,Solver.DOFs*q+2] += 0 *  IP_weight * el_size
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += 0 *  IP_weight * el_size 
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q+2] += 0 *  IP_weight * el_size
                        Jloc[Solver.DOFs*p+2,Solver.DOFs*q] += 0 *  IP_weight * el_size  
                        Jloc[Solver.DOFs*p+2,Solver.DOFs*q+1] += 0 *  IP_weight * el_size  
                   
                    ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Floc[Solver.DOFs*p + 2] += s*alpha*vkin * (dSWDdT * T_prev - SWD)  * basis[p] * IP_weight * el_size ##F_C
                    ### F_P and F_T are zeros for now, but might have to be changed later (e.g. SW radiation)       
                    
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1 
                    ###Assemble global mass matrice
                    Solver.Mass[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Mloc[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Mass[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Mloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Solver.Mass[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob+2] += Mloc[Solver.DOFs*p+2,Solver.DOFs*q+2]
                    
                    ###Assemble global stiffness matrice
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob+2] += Kloc[Solver.DOFs*p,Solver.DOFs*q+2]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+2] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q+2]
                    Solver.Stiffness[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p+2,Solver.DOFs*q]
                    Solver.Stiffness[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p+2,Solver.DOFs*q+1]
                    
                    ###Assemble global jacobian matrice
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob+2] += Jloc[Solver.DOFs*p,Solver.DOFs*q+2]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+2] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q+2]
                    Solver.Jac[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p+2,Solver.DOFs*q]
                    Solver.Jac[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p+2,Solver.DOFs*q+1]
                    
                ###Assemble global force matrice
                Solver.Force_bulk[Solver.DOFs*p_glob+2] += Floc[Solver.DOFs*p+2]
        
        ###Initialize Force vector to bulk force vector (force vector will be modified later if Neumann/Robin BCs are prescribed)
        Solver.Force[:] = Solver.Force_bulk[:]
      
        ###Lump mass matrix if needed
        if Solver.masslumping == True:
            for j in range(Solver.DOFs*Solver.parent.snowpack.mesh.numberofNodes):
                x = sum(Solver.Mass[j,:])
                Solver.Mass[j,:] = 0
                Solver.Mass[j,j] = x
            
            
######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######
######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END ASSEMBLING BULK MATRICES FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######
######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ######  

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING NEUMANN AND/OR ROBIN BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleNeumannRobinBC(Solver):
        ### loop on BC nodes
        for bcs in Solver.parent.snowpack.mesh.boundaries: ##loop on boundary nodes
            tag = bcs[0] ## top/bottom boundary ?
            node_number = bcs[1].numberinMesh ## node number (/!\ start from 1)
            xnode = bcs[1].pos ## node position (/!\ start from 1)
            
            ### For vapor
            if Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Given Rhov':
                pass ### If BC at considered node is prescribed as No Flux (natural) or Dirichlet we don't do anything and go see BC on T
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Rhov Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_vapor(node_number) ### get the flux from BC section of sif             
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(xnode)
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)] += Flux/(1-Phii) ### Add the neumann flux to the force vector                                 
            else:
                print('No boundary condition given for vapor on BC ' + tag + '. Setting it to no flux.')

            
            ### For heat
            if Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Adiabatic' or Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Dirichlet we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Heat Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the flux from BC section of sif  
                print('Heat Flux =', Flux)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(xnode)
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)+1] += Flux/rhoCeff                                 
            else:
                print('No boundary condition given for heat on BC ' + tag + '. Setting it to adiabatic.')
        
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING NEUMANN AND/OR ROBIN BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~   FUNCTION TRANSFORMING THE SYSTEM MdUdt + KU = F in AU = B  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleTimeStepping(Solver):
        if Solver.parent.timestepping_method == "theta":
            Solver.LHS_matrix = (Solver.Mass/(Solver.execsolver*Solver.parent.timestepsize) + Solver.parent.theta * Solver.Stiffness) + Solver.parent.theta * Solver.Jac
            Solver.RHS_vector = np.dot((Solver.Mass/(Solver.execsolver*Solver.parent.timestepsize) + (Solver.parent.theta-1)*Solver.Stiffness), Solver.solution_prev) + Solver.Force 
            Solver.RHS_vector += np.dot(Solver.parent.theta*Solver.Jac, Solver.solution_current_it)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~   END FUNCTION TRANSFORMING THE SYSTEM MdUdt + KU = F in AU = B  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING DIRICHLET BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleDirichletBC(Solver):
        ### loop on BC nodes
        for bcs in Solver.parent.snowpack.mesh.boundaries: ##loop on boundary nodes
            tag = bcs[0] ## top/bottom boundary ?
            node_number = bcs[1].numberinMesh ## node number (/!\ start from 1)
            xnode = bcs[1].pos ## node position (/!\ start from 1)
            
            ### For vapor
            if Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Rhov Flux':
                pass ### If BC at considered node is prescribed as No Flux (natural) or Neumann we don't do anything and go see BC on T
            elif Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Given Rhov':
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ### Below some ugly hard-code when the BC on rhov is a dirichlet BC that is a function of T, the latter being also a dirichlet BC (e.g. rhov_top = rhovsat(T_top))
                if (Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature') and (type(Solver.parent.snowpack.BC_tag[tag].functionRhov)==str) and ('Temperature' in Solver.parent.snowpack.BC_tag[tag].functionRhov.split('_')[2:]):
                    input_fields = Solver.parent.snowpack.BC_tag[tag].functionRhov.split('_')[2:] #### T is not necessarily in the only input of the function given rhov dirichlet at the BC !
                    if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                        raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')
                    input_values = []
                    for input_field in input_fields:
                        if input_field == 'Time':
                            input_value = Solver.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                        else:
                            if Solver.parent.snowpack.field_perm[input_field].NodEl == 'Nodal':
                                if input_field == 'Temperature': ### THE IMPORTANT SEQUENCE: PROVIDE THE DIRICHLET BC ON T AS ARGUMENT !!
                                    input_value = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number)
                                else:
                                    input_value = Solver.parent.snowpack.field_perm[input_field].value_prev_it[NodeNumber-1] ###If the input_field is nodal and is not T, we simply get its value at previous iteration at the considered BC node                                                                   
                            elif Solver.parent.snowpack.field_perm[input_field].NodEl == 'Elemental':
                                if tag == 'bot':
                                    input_value = Solver.parent.snowpack.field_perm[input_field].value_prev_it[NodeNumber-1] ##If the input field is elemental and we are at bottom node we take value of first element
                                elif tag == 'top':
                                    input_value = Solver.parent.snowpack.field_perm[input_field].value_prev_it[NodeNumber-2] ##If the input field is elemental and we are at top node we take value of last element
                        input_values.append(input_value)
                    RhovBC = eval(Solver.parent.snowpack.BC_tag[tag].functionRhov.split('_')[0]+'.'+Solver.parent.snowpack.BC_tag[tag].functionRhov.split('_')[0]+'_'+Solver.parent.snowpack.BC_tag[tag].functionRhov.split('_')[1])(*input_values)
                else:
                #### Temperature does not intervene in calculation of dirichlet rhov -> END OF UGLY HARD-CODE (we can directly use the method get_BC_vapor to recover the BC value)   
                ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                       
                    RhovBC = Solver.parent.snowpack.BC_tag[tag].get_BC_vapor(node_number) ### get the Dirichlet BC from BC section of sif             
               
                ###Modify the LHS matrix and RHS vector to add the Dirichlet condition
                Solver.RHS_vector[Solver.DOFs*(node_number-1)] = RhovBC ### Add the dirichlet value to the RHS vector
                print('Rhov at BC', tag,'=',RhovBC)
                Solver.LHS_matrix[Solver.DOFs*(node_number-1),:] = 0 ### All the terms of the line of the LHS matrix corresponding to the considered BC node are zeros...  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1), Solver.DOFs*(node_number-1)] = 1 ### ... except the one on the diagonal that is a 1
            else:
                print('No boundary condition given for vapor on BC ' + tag + '. Setting it to no flux.')

            
            ### For heat
            if Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Adiabatic' or Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Heat Flux':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Neumann we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature':
                TemperatureBC = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the T from BC section of sif  
                ###Modify the LHS matrix and RHS vector to add the Dirichlet condition
                Solver.RHS_vector[Solver.DOFs*(node_number-1)+1] = TemperatureBC ### Add the dirichlet value to the RHS vector  
                print('T at BC', tag,'=',TemperatureBC)
                Solver.LHS_matrix[Solver.DOFs*(node_number-1)+1,:] = 0 ### All the terms of the line of the LHS matrix corresponding to the considered BC node are zeros...  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1)+1, Solver.DOFs*(node_number-1)+1] = 1 ### ... except the one on the diagonal that is a 1f                                 
            else:
                print('No boundary condition given for heat on BC ' + tag + '. Setting it to adiabatic.')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING DIRICHLET BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~   FUNCTION TO CALCULATE THE RESIDUAL OF THE MAIN SYSTEM  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def GetResidual(Solver):
        if Solver.parent.timestepping_method == "theta":
            Left_matrix = (Solver.Mass / (Solver.execsolver*Solver.parent.timestepsize)) + Solver.parent.theta * Solver.Stiffness 
            Right_vector = np.dot((Solver.Mass / (Solver.execsolver*Solver.parent.timestepsize)) + (Solver.parent.theta-1) * Solver.Stiffness, Solver.solution_prev) + Solver.Force_bulk 
            Solver.residual = np.dot(Left_matrix, Solver.solution)-Right_vector
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~   END OF FUNCTION TO CALCULATE THE RESIDUAL OF THE MAIN SYSTEM  ~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####             
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    
    ### The three solution fields of this solver:
    Rhov= Solver.parent.snowpack.field_perm['Rhov']
    T = Solver.parent.snowpack.field_perm['Temperature']
    c = Solver.parent.snowpack.field_perm['DepositionRate']
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    Rhov.value_prev_tsp[:] = Rhov.value[:]
    T.value_prev_tsp[:] = T.value[:]
    c.value_prev_tsp[:] = c.value[:]   

    # Initialize PDE
    
    ###Construct solution vector from Rhov and T fields as (Rhov0, T0, Rhov1, T1,....,Rhovn, Tn)
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Solver.solution_prev[Solver.DOFs*i] = Rhov.value[i]
        Solver.solution_prev[Solver.DOFs*i+1] = T.value[i]
        Solver.solution_prev[Solver.DOFs*i+2] = 0
    
    ###Initialize solution_current_it for non-linear iterations loop
    Solver.solution_current_it[:] =  Solver.solution_prev[:]
    
    ###Initialize counter and get parameters for non-linear iterations loop
    nonlin_count = 1
    nonlin_converged = False
    nonlin_max_it = Solver.nonlin_max_it
    nonlin_thres = Solver.nonlin_thres
    nonlin_newton_after_it = Solver.nonlin_newton_after_it
    nonlin_newton_after_thres = Solver.nonlin_newton_after_thres
    
#    np.set_printoptions(precision=2)

####~~~~~~~~~~~~~~~~~~~~~~~~~ Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    while (not nonlin_converged and nonlin_count <= nonlin_max_it):
        
        AssembleBulkMatrix(Solver)
        
        AssembleNeumannRobinBC(Solver)    
      
        #Construct AU=B system (here we do not use the method AssembleTimeStepping of the class solver_t because of equation on c that is unsual with no time derivative)
        AssembleTimeStepping(Solver)

        # Add Dirichlet
        AssembleDirichletBC(Solver)
    
        # Solve the PDE system        
        Solver.linSolve()
        NormChange = Solver.NormChange()
        
        #### Store the solution vector as solution_current_it for next iteration
        Solver.solution_current_it = copy.deepcopy(Solver.solution)
        
        #### Update the value_prev_it of solution fields to get the correct values of materials/BF/fields in Assembling of bulk matrix
        for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
            Rhov.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i] 
            T.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i+1] 
            c.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i+2]

        nonlin_converged = NormChange <= nonlin_thres
        print('ITERATION N°:', nonlin_count, 'NORM CHANGE =', NormChange)
        nonlin_count += 1
####~~~~~~~~~~~~~~~~~~~~~~~~~ End of Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Calculate the residual
    GetResidual(Solver)  
#### Update the residual of solution fields for next timestep
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.residual[i] = Solver.residual[Solver.DOFs*i] 
        T.residual[i] = Solver.residual[Solver.DOFs*i+1] 
        c.residual[i] = Solver.residual[Solver.DOFs*i+2] 

#####~~~~~ After convergence, recalculate deposition rate in the weak sense but with lumped mass matrix to remove oscillations ~~~~~####
    M_CC = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
    M_CC_lumped = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
   
    for node in Solver.parent.snowpack.mesh.nodes:
        p_glob = node.numberinMesh -1
        for node_2 in Solver.parent.snowpack.mesh.nodes:
            q_glob = node_2.numberinMesh - 1 
            M_CC[p_glob,q_glob] = Solver.Mass[Solver.DOFs*p_glob+2,Solver.DOFs*q_glob+2]
    RHS_Vector = np.dot(M_CC, c.value_prev_it[:])
    ###Lump the M_CC mass matrix
    for j in range(Solver.parent.snowpack.mesh.numberofNodes):
        x = sum(M_CC[j,:])
        M_CC_lumped[j,:] = 0
        M_CC_lumped[j,j] = x
    ###Solve the system
    Solution_c_lumped = np.linalg.solve(M_CC_lumped, RHS_Vector)
#        
####~~~~~ Update the solution fields for next timestep ~~~~~#### 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.value[i] = Solver.solution[Solver.DOFs*i] 
        T.value[i] = Solver.solution[Solver.DOFs*i+1] 
        c.value[i] = Solution_c_lumped[i] 
    
#### For this solver, it is necessary to multiply residual at boundaries (they are zeros elsewhere) by resp. (1-Phii) and (rhoC)eff to get true fluxes that could be used to check energy conservation
    Rhov.residual[0] = Rhov.residual[0] *(1-Solver.parent.snowpack.field_perm['Phii'].value[0]) 
    Rhov.residual[-1] = Rhov.residual[-1] *(1-Solver.parent.snowpack.field_perm['Phii'].value[-1])
    T.residual[0] = T.residual[0] * Solver.parent.snowpack.material_perm['rhoCeff'].get_material(Solver.parent.snowpack.mesh.nodes[0].pos)
    T.residual[-1] = T.residual[-1] * Solver.parent.snowpack.material_perm['rhoCeff'].get_material(Solver.parent.snowpack.mesh.nodes[-1].pos)      

    
