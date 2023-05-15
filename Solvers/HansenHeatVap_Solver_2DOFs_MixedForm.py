import FEM_classes
import numpy as np
import copy
import MatLaw
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

###This solver solves the heat equation of the Hansen system to get the nodal T field but for purpose of energy conservation we write the heat equation in mixed form:
### accumulation term is written in terms of enthalpy and flux divergence in terms of temperature.
### As a consequence this solver has two DOFs (H,T) and 2 Eqs : Eq 1 is heat equation in mixed form; Eq 2 is relationship between enthalpy and temperature
### Eq1 -> dH/dt -d/dz[(Deff Lm drhov_sat/dT + keff)dT/dz] = 0
### Eq2 -> H - H0 = (1 -Phii)*Lm*(rhovsat(T) - rhovsat(T0)) + (rhoC)eff*(T - T0) with H0 = 0 for T0 = Tfus 
### BE CAREFUL, this formulation is valid ONLY FOR DRY SNOW

def HansenHeatVap(Solver):
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING GLOBAL MATRICES FROM LOCAL MATRICES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleBulkMatrix(Solver):
        Solver.cleanMatrices()
        ## Get the needed constants only once
        Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition
        Tfus = Solver.parent.snowpack.constant_perm['Tfus'].get_cst() ### Fusion Temperature (K)
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
                Deff = Solver.parent.snowpack.material_perm['Deff'].get_material(IPpos)
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
                SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(IPpos)
                dSWDdT = Solver.parent.snowpack.material_perm['dSWDdT'].get_material(IPpos)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(IPpos)
                Keff = Solver.parent.snowpack.material_perm['Keff'].get_material(IPpos)
                T_prev= Solver.parent.snowpack.field_perm['Temperature'].get_value_prev_it(IPpos) ###Value of temperature obtained from previous iteration      
                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )             
                        ###~~~~~~LOCAL MASS MATRIX~~~~~~###
                        Mloc[Solver.DOFs*p,Solver.DOFs*q] += basis[p] * basis[q] * IP_weight * el_size  #For DOFs 1, here H
                        ## No terms in mass matrix for Eq.2
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                        Kloc[Solver.DOFs*p,Solver.DOFs*q+1] += (Deff * Lm * dSWDdT + Keff) * dbasis[q] * dbasis[p] *  IP_weight * el_size ##K_{1T}
                        Kloc[Solver.DOFs*p+1,Solver.DOFs*q] += basis[p] * basis[q] *  IP_weight * el_size ##K_{2H}
                        Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += - ((1 - Phii) * Lm * dSWDdT + rhoCeff) * basis[p] * basis[q] *  IP_weight * el_size ##K_{2T}                        
                        ###~~~~~~LOCAL JACOBIAN MATRIX~~~~~~### To implement later if Newton linearization needed
                        Jloc[Solver.DOFs*p,Solver.DOFs*q+1] += 0 *  IP_weight * el_size 
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q] += 0 *  IP_weight * el_size 
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += 0 *  IP_weight * el_size 
                        ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    ## No source in Eq. 1
                    Floc[Solver.DOFs*p+1] += ( (1 - Phii) * Lm * (SWD - dSWDdT * T_prev) - rhoCeff * Tfus)  * basis[p] * IP_weight * el_size ##F_{2}
#                    Floc[Solver.DOFs*p+1] += (rhoCeff *(T_prev - Tfus) + Lm*(1-Phii)*(SWD - SWD_Tfus))   * basis[p] * IP_weight * el_size ##F_{2}
#                    Floc[Solver.DOFs*p+1] += -((1 - Phii) * Lm * dSWDdT + rhoCeff)*T_prev   * basis[p] * IP_weight * el_size ##F_{2}
           
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1 
                    ###Assemble global mass matrice
                    Solver.Mass[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Mloc[Solver.DOFs*p,Solver.DOFs*q]                    
                    ###Assemble global stiffness matrice
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p,Solver.DOFs*q+1]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    ###Assemble global jacobian matrice
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p,Solver.DOFs*q+1]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                ###Assemble global force matrice
                Solver.Force_bulk[Solver.DOFs*p_glob+1] += Floc[Solver.DOFs*p+1]                
        
        ###Initialize Force vector to bulk force vector (force vector will be modified later if Neumann/Robin BCs are prescribed)
        Solver.Force[:] = Solver.Force_bulk[:]
        ###Lump mass matrix if needed, which also implies to lump the part K2H of the stiffness matrix as the latter has the form of a mass matrix on H
        if Solver.masslumping == True:
            for j in range(Solver.parent.snowpack.mesh.numberofNodes):
                x = sum(Solver.Mass[Solver.DOFs*j,:])
                Solver.Mass[Solver.DOFs*j,:] = 0
                Solver.Mass[Solver.DOFs*j,Solver.DOFs*j] = x
                if j ==0:
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j] += Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j+2]
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j+2] = 0
                elif j == Solver.parent.snowpack.mesh.numberofNodes -1:
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j] += Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j-2]
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j-2] = 0
                else:
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j] += Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j-2] + Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j+2]
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j-2] = 0
                    Solver.Stiffness[Solver.DOFs*j+1,Solver.DOFs*j+2] = 0
            
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
                TemperatureBC = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the flux from BC section of sif  
                ### For this solver we also need to convert the Dirichlet BCs on T into Dirichlet BC on enthalpy
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(xnode)
                Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition 
                Tfus = Solver.parent.snowpack.constant_perm['Tfus'].get_cst() ### Fusion Temperature (K)
                SWD_TemperatureBC = eval(Solver.parent.snowpack.material_perm['SWD'].function.split('_')[0]+'.'+Solver.parent.snowpack.material_perm['SWD'].function.split('_')[0]+'_'+Solver.parent.snowpack.material_perm['SWD'].function.split('_')[1])(TemperatureBC)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(xnode)
                
                HBC = (1-Phii)*Lm*(SWD_TemperatureBC) + rhoCeff * (TemperatureBC - Tfus)
                
                ###Modify the LHS matrix and RHS vector to add the Dirichlet condition
                Solver.RHS_vector[Solver.DOFs*(node_number-1)] = HBC ### Add the dirichlet value to the RHS vector  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1),:] = 0 ### All the terms of the line of the LHS matrix corresponding to the considered BC node are zeros...  
                Solver.LHS_matrix[Solver.DOFs*(node_number-1), Solver.DOFs*(node_number-1)] = 1 ### ... except the one on the diagonal that is a 1f
            else:
                print('No boundary condition given for heat on BC ' + tag + '. Setting it to adiabatic.')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING DIRICHLET BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####  
               
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####
####///~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\\\####             
####//////////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    
    ### The solution field of this solver:
    H = Solver.parent.snowpack.field_perm['Enthalpy']
    T = Solver.parent.snowpack.field_perm['Temperature']
    ### BEFORE ANY MODIFICATION, store the field in prev_tsp for the method get_value_prev_tsp of the class field_t
    H.value_prev_tsp[:] = H.value[:]
    T.value_prev_tsp[:] = T.value[:]

    ###~~~~~~~~~ Initialization of enthalpy ~~~~~~~~~###   
    #### In addition to T, this solver requires the variable enthalpy to be initialized when this sover is visited for the first time    
    ### We write a matrice system on the form M_HH H + K_HT T = F_H
    if Solver.parent.timestep == 1:
        M_HH = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
        K_HT = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
        F_H = np.zeros((Solver.parent.snowpack.mesh.numberofNodes))
        ## Get the needed constants only once
        Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition        
        Tfus = Solver.parent.snowpack.constant_perm['Tfus'].get_cst() ### Fusion Temperature (K)
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()
            Khtloc = np.zeros((el.numberofnodes,el.numberofnodes))
            Mhloc = np.zeros((el.numberofnodes,el.numberofnodes))
            Fhloc = np.zeros((el.numberofnodes))         
            for IP in range(el.numberofIP):
                ### Get parameters of IPs necessary for Gaussian integration
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)
                
                ### Get all material properties necessary to calculate Stiffness and force matrices at IP
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
                SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(IPpos)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(IPpos)
                               
                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )              
                        Mhloc[p,q] += basis[p] * basis[q] *  IP_weight * el_size ##M_{HH}
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                        Khtloc[p,q] += - rhoCeff * basis[p] * basis[q] *  IP_weight * el_size ##K_{HT}                        
                        ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Fhloc[p] += ( (1 - Phii) * Lm * SWD - rhoCeff * Tfus )  * basis[p] * IP_weight * el_size ##F_{H}
    
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1 
                    ###Assemble global mass matrice
                    M_HH[p_glob,q_glob] += Mhloc[p,q]                    
                    ###Assemble global stiffness matrice
                    K_HT[p_glob,q_glob] += Khtloc[p,q]
                    ###Assemble global force matrice
                F_H[p_glob] += Fhloc[p]                
        ###If mass matrix is lumped in the resolution of the PDE, then it must also be lumped when initializing H for consistency
        if Solver.masslumping == True:
            for j in range(Solver.parent.snowpack.mesh.numberofNodes):
                x = sum(M_HH[j,:])
                M_HH[j,:] = 0
                M_HH[j,j] = x
        ###Now we compute Hinit by solving the system M_HH H + K_HT T = F_H
        RHS_Vector = F_H - np.dot(K_HT, T.value[:])
        H.value_prev_tsp[:] = np.linalg.solve(M_HH, RHS_Vector)
        H.value[:] = np.linalg.solve(M_HH, RHS_Vector)
                        
        ###~~~~~~~~~ End initialization of enthalpy ~~~~~~~~~###
        
    ### Although it is not the variable of the PDE, We will also uptade the Rhov field from Rhov = Rhovsat within this solver
    Rhov = Solver.parent.snowpack.field_perm['Rhov']
    ### BEFORE ANY MODIFICATION, store the field in prev_tsp for the method get_value_prev_tsp of the class field_t
    Rhov.value_prev_tsp[:] = Rhov.value[:]
    
    # Initialize PDE
    Q = 0 ## For now, no source in this model but to adapt later
    ###Construct solution vector from H and T fields  (H0, T0, H1, T1, ..., Hn, Tn)
    ### BE CAREFUL that by choice T is DOFs 2 while H is DOFs 1 although not exported as a model variable
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Solver.solution_prev[Solver.DOFs*i] = H.value[i]
        Solver.solution_prev[Solver.DOFs*i+1] = T.value[i]
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
        Solver.AssembleTimeStepping()

        # Add Dirichlet
        AssembleDirichletBC(Solver)
    
        # Solve the PDE system        
        Solver.linSolve()
        NormChange = Solver.NormChange()

        #### Store the solution vector as solution_current_it for next iteration
        Solver.solution_current_it = copy.deepcopy(Solver.solution)

        nonlin_converged = NormChange <= nonlin_thres
        #### Update the value_prev_it of solution fields to get the correct values of materials/BF/fields in Assembling of bulk matrix
        if (not nonlin_converged and nonlin_count <= nonlin_max_it):
            for i in range(Solver.parent.snowpack.mesh.numberofNodes):
                H.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i]
                T.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i+1]

        print('ITERATION N°:', nonlin_count, 'NORM CHANGE =', NormChange)
        nonlin_count += 1
####~~~~~~~~~~~~~~~~~~~~~~~~~ End of Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####


#### Update the solution field for next timestep 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes):
        H.value[i] = Solver.solution[Solver.DOFs*i]
        T.value[i] = Solver.solution[Solver.DOFs*i+1] 
        ### Update the rhov field from Rhov = Rhov_sat (assumption of the Hansen system)
        node_pos = Solver.parent.snowpack.mesh.nodes[i].pos
        Rhov.value[i] = Solver.parent.snowpack.material_perm['SWD'].get_material(node_pos)

#### Calculate the residual
    Solver.GetResidual()  
#### Update the residual of solution fields for next timestep
    for i in range(Solver.parent.snowpack.mesh.numberofNodes):
        T.residual[i] = Solver.residual[Solver.DOFs*i] ### BE CAREFUL, heat eq. is eq. 1 so residus in terms of heat are on lines Solver.DOFs*i but for consistency with Calonne solver we store them in T.residual even if T is DOFs 2 


    
