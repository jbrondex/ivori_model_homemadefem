import FEM_classes
import MatLaw
import numpy as np
import copy

###This solver diagnostizes the deposition rate for the Hansen system (c = -(1-Phii)dSWDdT dTdt + d/dz(Deff dSWDdT dTdz) in the weak sense
### from T calculated by the HansenHeatVap_Solver which must be placed right before this sover, and the T from the previous timestep
### This is equivalent to solve the system MC = - M* dT/dt - K T + F 

def DepositionRateHansen(Solver):
                  
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING GLOBAL MATRICES FROM LOCAL MATRICES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleBulkMatrix(Solver):
        Solver.cleanMatrices()
        # Perform Bulk Assemble
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()
            Kloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Mloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Fstarloc = np.zeros((Solver.DOFs*el.numberofnodes)) ## /!\ Particular force vector 
            for IP in range(el.numberofIP):
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)
                
                ### Get all material properties necessary to calculate Stiffness and force matrices at IP
                Deff = Solver.parent.snowpack.material_perm['Deff'].get_material(IPpos)
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
                dSWDdT = Solver.parent.snowpack.material_perm['dSWDdT'].get_material(IPpos) 
                ### Get the Rhovsat for current and previous tsp from T of previous and current tsp
                SWD_current = Solver.parent.snowpack.material_perm['SWD'].get_material(IPpos) ##Easy for current tsp
                ### More difficult for previous tsp: would require adaptation of get_material method that I don't want to do for this unique solver -> 
                ###-> Hard-coded here for now but might have to be improved later
                function_name = Solver.parent.snowpack.material_perm['SWD'].function
                liste = function_name.split('_')
                SWD_prev_tsp = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(T.get_value_prev_tsp(IPpos))
                
                
                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to )             
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                        Kloc[Solver.DOFs*p,Solver.DOFs*q] += (Deff * dSWDdT)* dbasis[q] * dbasis[p] *  IP_weight * el_size                                
                        ###~~~~~~LOCAL MATRIX APPLYING ON C (similar to a usual mass matrix)~~~~~~###
                        Mloc[Solver.DOFs*p, Solver.DOFs*q] += basis[p] * basis[q] * IP_weight * el_size  
                    ###~~~~~~LOCAL FORCE MATRIX Corresponding to the time derivative on rhov=rhovsat(T)~~~~~~###
                    Fstarloc[Solver.DOFs*p] += - (1 - Phii) * ((SWD_current - SWD_prev_tsp)/Delta_t) * basis[p] * IP_weight * el_size 
                        
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1                    
                    ###Assemble global stiffness matrice
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p,Solver.DOFs*q]
                    ###Assemble global matrice applying on C
                    Solver.Mass[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Mloc[Solver.DOFs*p,Solver.DOFs*q]     
                ###Assemble global force matrice
                Solver.Force[Solver.DOFs*p_glob] += Fstarloc[Solver.DOFs*p] 
        
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
            
            if Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'Given Deposition Rate':
                pass ### If BC at considered node is prescribed as No Flux (natural) or Dirichlet we don't do anything and go see BC on T
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'Mass Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_depositionrate(node_number) ### get the flux from BC section of sif             
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)] += Flux ### Add the neumann flux to the force vector                                 
            else:
                print('No boundary condition given for Deposition Rate on BC ' + tag + '. Setting it to no flux.')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING NEUMANN AND/OR ROBIN BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
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

            if Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'Mass Flux':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Neumann we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_DepositionRate == 'Given Deposition Rate':
                DepositionRateBC = Solver.parent.snowpack.BC_tag[tag].get_BC_depositionrate(node_number) ### get the flux from BC section of sif  
                ###Modify the LHS matrix and RHS vector to add the Dirichlet condition
                Solver.RHS_vector[Solver.DOFs*(node_number-1)] = DepositionRateBC ### Add the dirichlet value to the RHS vector  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1),:] = 0 ### All the terms of the line of the LHS matrix corresponding to the considered BC node are zeros...  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1), Solver.DOFs*(node_number-1)] = 1 ### ... except the one on the diagonal that is a 1f                                 
            else:
                print('No boundary condition given for Deposition Rate on BC ' + tag + '. Setting it to No Flux.')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING DIRICHLET BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####    
                

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####             
        
    ### The solution field of this solver:
    c = Solver.parent.snowpack.field_perm['DepositionRate']
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    c.value_prev_tsp[:] = c.value[:]
    
    ###The time elapsed since last execution of the solver
    Delta_t = Solver.execsolver*Solver.parent.timestepsize
    ### get the T field and calculate time derivative:
    T = Solver.parent.snowpack.field_perm['Temperature']

    ### Construct the Global Matrice M,K and F
    AssembleBulkMatrix(Solver)    
    ### Construct the Force matrix to account for Neumann BC    
    AssembleNeumannRobinBC(Solver)
    
    #Construct MC=B system from system MC = -KT + F (with F = Fdrhovsatdt + F_Neumann) 
    Solver.LHS_matrix = Solver.Mass
    Solver.RHS_vector = - np.dot(Solver.Stiffness, T.value) + Solver.Force 
    
    # Add Dirichlet (not very meaningfull physically in this solver, but we leave this possibility anyway to the user with a warning message)
    AssembleDirichletBC(Solver)
        
    ##Solve the system
    Solver.linSolve()
    #### Store the solution vector also in solution_current_it for consistency (even if never used since solution is a diagnostic)
    Solver.solution_current_it = copy.deepcopy(Solver.solution)
        
    #### Update the solution field 
    c.value = copy.deepcopy(Solver.solution)  
    c.value_prev_it = copy.deepcopy(c.value)
        
