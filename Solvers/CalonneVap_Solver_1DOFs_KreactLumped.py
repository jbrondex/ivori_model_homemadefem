import FEM_classes
import numpy as np
import MatLaw
import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def CalonneVap(Solver):
    
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
                Deff = Solver.parent.snowpack.material_perm['Deff'].get_material(IPpos)
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
                cIP = Solver.parent.snowpack.field_perm['DepositionRate'].get_value(IPpos)
                
                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )
                        
                        ###~~~~~~LOCAL MASS MATRIX~~~~~~###
                        Mloc[Solver.DOFs*p,Solver.DOFs*q] += (1-Phii) * basis[p] * basis[q] * IP_weight * el_size  #For DOFs 1, here Rhov
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~### ### In this solver stiffness matrice on reaction terms is disentangle from stiffness matrice on diffusion terms (to enable lumping of the former)
                        Kloc[Solver.DOFs*p,Solver.DOFs*q] += Deff * dbasis[q] * dbasis[p] *  IP_weight * el_size ##K_{PP} for diffusion term
                        ###~~~~~~LOCAL JACOBIAN MATRIX~~~~~~### To implement later if Newton linearization needed
                        Jloc[Solver.DOFs*p,Solver.DOFs*q] += 0 *  IP_weight * el_size 
                   
                    ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Floc[Solver.DOFs*p] += - cIP * basis[p] * IP_weight * el_size ##F_P
                        
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
            
            ### For vapor
            if Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Given Rhov':
                pass ### If BC at considered node is prescribed as No Flux (natural) or Dirichlet we don't do anything and go see BC on T
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Rhov Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_vapor(node_number) ### get the flux from BC section of sif             
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(xnode)
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)] += Flux ### Add the neumann flux to the force vector                                 
            else:
                print('No boundary condition given for vapor on BC ' + tag + '. Setting it to no flux.')
        
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
            
            ### For vapor
            if Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'No Flux' or Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Rhov Flux':
                pass ### If BC at considered node is prescribed as No Flux (natural) or Neumann we don't do anything and go see BC on T
            elif Solver.parent.snowpack.BC_tag[tag].BCType_Rhov == 'Given Rhov':
                ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ### Below some ugly hard-code when the BC on rhov is a dirichlet BC that is a function of T, the latter being also a dirichlet BC (e.g. rhov_top = rhovsat(T_top)) (Normally not a problem when heat and vap decoupled)
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
    
    ### The two solution fields of this solver:
    Rhov= Solver.parent.snowpack.field_perm['Rhov']
    c = Solver.parent.snowpack.field_perm['DepositionRate']
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    Rhov.value_prev_tsp[:] = Rhov.value[:]
    c.value_prev_tsp[:] = c.value[:]   

    # Initialize PDE
    
    ###Construct solution vector from Rhov 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Solver.solution_prev[Solver.DOFs*i] = Rhov.value[i]

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

####~~~~~~~~~~~~~~~~~~~~~~~~~ Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#### NORMALLY WHEN HEAT AND VAP DECOUPLED, VAP IS LINEAR SO ONLY ONE ITERATION EXPECTED
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
            Rhov.value_prev_it[i] = Solver.solution_current_it[Solver.DOFs*i]  

        nonlin_converged = NormChange <= nonlin_thres
        print('ITERATION N°:', nonlin_count, 'NORM CHANGE =', NormChange)
        nonlin_count += 1
####~~~~~~~~~~~~~~~~~~~~~~~~~ End of Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Update the solution fields for next timestep 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.value[i] = Solver.solution[Solver.DOFs*i] 

   

#### Calculate the residual
    GetResidual(Solver)  
#### Update the residual of solution fields for next timestep
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.residual[i] = Solver.residual[Solver.DOFs*i] 


    #### Update the deposition rate 
    for node in Solver.parent.snowpack.mesh.nodes: 
        s = Solver.parent.snowpack.bodyforce_perm['s'].get_BF(node.pos)
        alpha = Solver.parent.snowpack.bodyforce_perm['alpha'].get_BF(node.pos)
        vkin = Solver.parent.snowpack.bodyforce_perm['vkin'].get_BF(node.pos) 
        SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(node.pos) 
        
        c.value[node.numberinMesh-1] = s*alpha*vkin*(Rhov.value[node.numberinMesh-1] - SWD)
    
    c.value_prev_it[:] = c.value[:]
        
#    ###~~~~~~ Diagnostic of c in the weak sense: M_CC * C = M_CP * P + F ~~~~~~###
#    M_CC = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
#    M_CP = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
#    Force = np.zeros((Solver.parent.snowpack.mesh.numberofNodes))
#    for el in Solver.parent.snowpack.mesh.elements:
#        el_size = el.get_elsize()
#        Mccloc = np.zeros((el.numberofnodes,el.numberofnodes))
#        Mcploc = np.zeros((el.numberofnodes,el.numberofnodes))
#        Floc = np.zeros((el.numberofnodes))
#        
#        for IP in range(el.numberofIP):
#            IP_weight = el.get_IPweight(IP)
#            basis, dbasis, ddbasis = el.get_bases(IP)
#            IPpos = el.get_IPpos(IP)
#                
#            ### Get all material properties necessary to calculate Stiffness and force matrices at IP
#            s = Solver.parent.snowpack.bodyforce_perm['s'].get_BF(IPpos)
#            alpha = Solver.parent.snowpack.bodyforce_perm['alpha'].get_BF(IPpos)
#            vkin = Solver.parent.snowpack.bodyforce_perm['vkin'].get_BF(IPpos) 
#            SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(IPpos) 
#
#            for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
#                for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions)
#                    ###~~~~~~LOCAL MASS MATRIX~~~~~~###
#                    Mccloc[p,q] +=  basis[p] * basis[q] * IP_weight * el_size    
#                    Mcploc[p,q] +=  s * alpha * vkin * basis[p] * basis[q] * IP_weight * el_size  
#                ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
#                Floc[p] += -(s*alpha*vkin) * SWD * basis[p] * IP_weight * el_size ##F_P
#        
#        nodes = el.nodes
#        for p in range(el.numberofnodes):
#            p_glob = nodes[p].numberinMesh -1
#            for q in range(el.numberofnodes):
#                q_glob = nodes[q].numberinMesh - 1           
#                ###Assemble global mass matrice on unknown C
#                M_CC[p_glob,q_glob] += Mccloc[p,q]      
#                ###Assemble global mass matrice on  Delta P
#                M_CP[p_glob,q_glob] += Mcploc[p,q]   
#            ###Assemble global force matrice
#            Force[p_glob] += Floc[p]            
#            
#    ###To remove oscillations, the MCC matrix must be lumped anyway
#    for j in range(Solver.parent.snowpack.mesh.numberofNodes):
#        x = sum(M_CC[j,:])
#        M_CC[j,:] = 0
#        M_CC[j,j] = x
#    
#    #### Solve the obtained system
#    RHS_Vector = np.dot(M_CP, Rhov.value[:]) + Force[:]
#    Solution_c_diag = np.linalg.solve(M_CC, RHS_Vector)
#    
#    #### Update the deposition rate 
#    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
#        c.value[i] = Solution_c_diag[i]
#        c.value_prev_it[i] = Solution_c_diag[i] 
#
#


