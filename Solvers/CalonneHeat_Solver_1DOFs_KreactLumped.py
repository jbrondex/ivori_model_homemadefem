import FEM_classes
import numpy as np
import MatLaw
import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def CalonneHeat(Solver):
    
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
                ### Get parameters of IPs necessary for Gaussian integration
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)
                
                ### Get all material properties necessary to calculate Stiffness and force matrices at IP
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(IPpos)
                Keff = Solver.parent.snowpack.material_perm['Keff'].get_material(IPpos)
                Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition
                c = Solver.parent.snowpack.field_perm['DepositionRate'].get_value(IPpos)

                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )
                        
                        ###~~~~~~LOCAL MASS MATRIX~~~~~~###
                        Mloc[Solver.DOFs*p,Solver.DOFs*q] += rhoCeff * basis[p] * basis[q] * IP_weight * el_size #For DOFs 1, here T
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~### ### In this solver stiffness matrice on reaction terms is disentangle from stiffness matrice on diffusion terms (to enable lumping of the former)
                        Kloc[Solver.DOFs*p,Solver.DOFs*q] += Keff * dbasis[q] * dbasis[p] *  IP_weight * el_size  ##K_{TT} for diffusion term
                       
                        ###~~~~~~LOCAL JACOBIAN MATRIX~~~~~~### To implement later if Newton linearization needed
                        Jloc[Solver.DOFs*p,Solver.DOFs*q] += 0 *  IP_weight * el_size 
                   
                    ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Floc[Solver.DOFs*p] += Lm*c* basis[p] * IP_weight * el_size ##F_T
                        
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1 
                    ###Assemble global mass matrice
                    Solver.Mass[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Mloc[Solver.DOFs*p,Solver.DOFs*q]
                    
                    ###Assemble global stiffness matrice
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p,Solver.DOFs*q]
                 
                    ###Assemble global jacobian matrice
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p,Solver.DOFs*q]
                    
                ###Assemble global force matrice
                Solver.Force_bulk[Solver.DOFs*p_glob] += Floc[Solver.DOFs*p]
         
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
            
            ### For heat
            if Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Adiabatic' or Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Dirichlet we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Heat Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the flux from BC section of sif  
                print('Heat Flux =', Flux)
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)] += Flux                                
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
            Solver.LHS_matrix = (Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * (Solver.Stiffness)) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * Solver.Jac
            Solver.RHS_vector = np.dot((Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * (Solver.parent.theta-1) * (Solver.Stiffness)), Solver.solution_prev) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.Force 
            Solver.RHS_vector += np.dot((Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta*Solver.Jac, Solver.solution_current_it)
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
            
            ### For heat
            if Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Adiabatic' or Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Heat Flux':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Neumann we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature':
                TemperatureBC = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the T from BC section of sif  
                ###Modify the LHS matrix and RHS vector to add the Dirichlet condition
                Solver.RHS_vector[Solver.DOFs*(node_number-1)] = TemperatureBC ### Add the dirichlet value to the RHS vector  
                print('T at BC', tag,'=',TemperatureBC)
                Solver.LHS_matrix[Solver.DOFs*(node_number-1),:] = 0 ### All the terms of the line of the LHS matrix corresponding to the considered BC node are zeros...  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1), Solver.DOFs*(node_number-1)] = 1 ### ... except the one on the diagonal that is a 1                                
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
            Left_matrix = Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * (Solver.Stiffness) 
            Right_vector = np.dot((Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * (Solver.parent.theta-1) * (Solver.Stiffness)), Solver.solution_prev) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.Force_bulk 
            Solver.residual = (np.dot(Left_matrix, Solver.solution)-Right_vector) /  (Solver.execsolver*Solver.parent.timestepsize) ### Don't forget to divise by time interval if matrice system AX=B was multiplied by Delta_t in AssembleTimeStepping in order to get a true flux (i.e. quantity of energy/mass per unit of time)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~   END OF FUNCTION TO CALCULATE THE RESIDUAL OF THE MAIN SYSTEM  ~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
            

####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####             
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    
    ### The only solution fields of this solver:
    T = Solver.parent.snowpack.field_perm['Temperature']
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    T.value_prev_tsp[:] = T.value[:] 

    # Initialize PDE
    
    ###Construct solution vector from T field
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Solver.solution_prev[Solver.DOFs*i] = T.value[i]
    
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
      
        #Construct AU=B system
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
            T.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i]  

        nonlin_converged = NormChange <= nonlin_thres
        print('ITERATION N°:', nonlin_count, 'NORM CHANGE =', NormChange)
        nonlin_count += 1
####~~~~~~~~~~~~~~~~~~~~~~~~~ End of Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Update the solution fields for next timestep 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        T.value[i] = Solver.solution[Solver.DOFs*i] 

#### Calculate the residual
    GetResidual(Solver)  
#### Update the residual of solution fields for next timestep
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        T.residual[i] = Solver.residual[Solver.DOFs*i] 

