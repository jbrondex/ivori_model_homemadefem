import FEM_classes
import numpy as np
import MatLaw
import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def CoupledCalonneHeatVap(Solver):
    
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION ASSEMBLING GLOBAL MATRICES FROM LOCAL MATRICES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleBulkMatrix(Solver):
        Solver.cleanMatrices()
        Kglob_react = np.zeros((Solver.DOFs * Solver.parent.snowpack.mesh.numberofNodes, Solver.DOFs * Solver.parent.snowpack.mesh.numberofNodes))
        # Perform Bulk Assemble
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()
            Kloc = np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes))
            Kloc_react =  np.zeros((Solver.DOFs*el.numberofnodes,Solver.DOFs*el.numberofnodes)) ##We separate the stiffness matrix for diffusion and the one for reaction just for tests
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
                        Mloc[Solver.DOFs*p,Solver.DOFs*q] += (1-Phii) * basis[p] * basis[q] * IP_weight * el_size  #For DOFs 1, here Rhov
                        Mloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += rhoCeff * basis[p] * basis[q] * IP_weight * el_size #For DOFs 2, here T
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~### ### In this solver stiffness matrice on reaction terms is disentangle from stiffness matrice on diffusion terms (to enable lumping of the former)
                        Kloc[Solver.DOFs*p,Solver.DOFs*q] += Deff * dbasis[q] * dbasis[p] *  IP_weight * el_size ##K_{PP} for diffusion term
                        Kloc_react[Solver.DOFs*p,Solver.DOFs*q] += (s*alpha*vkin) * basis[p] * basis[q] *  IP_weight * el_size ##K_{PP} for reaction term 
                        Kloc_react[Solver.DOFs*p,Solver.DOFs*q+1] += - (s*alpha*vkin)* dSWDdT * basis[p] * basis[q] *  IP_weight * el_size ##K_{PT} for reaction term 
                        Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += Keff * dbasis[q] * dbasis[p] *  IP_weight * el_size  ##K_{TT} for diffusion term
                        Kloc_react[Solver.DOFs*p+1,Solver.DOFs*q+1] += (Lm*s*alpha*vkin) * dSWDdT * basis[p] * basis[q] *  IP_weight * el_size  ##K_{TT} for reaction term
                        Kloc_react[Solver.DOFs*p+1,Solver.DOFs*q] += - (Lm*s*alpha*vkin) * basis[p] * basis[q]*  IP_weight * el_size ##K_{TP} for reaction term
                        
                        ###~~~~~~LOCAL JACOBIAN MATRIX~~~~~~### To implement later if Newton linearization needed
                        Jloc[Solver.DOFs*p,Solver.DOFs*q] += 0 *  IP_weight * el_size 
                        Jloc[Solver.DOFs*p,Solver.DOFs*q+1] += 0 *  IP_weight * el_size
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1] += 0 *  IP_weight * el_size  
                        Jloc[Solver.DOFs*p+1,Solver.DOFs*q] += 0 *  IP_weight * el_size  
                   
                    ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Floc[Solver.DOFs*p] += ((s*alpha*vkin) * (- dSWDdT * T_prev + SWD) ) * basis[p] * IP_weight * el_size ##F_P
                    Floc[Solver.DOFs*p +1] += - ((Lm*s*alpha*vkin) * (- dSWDdT * T_prev + SWD) ) * basis[p] * IP_weight * el_size ##F_T
                        
            nodes = el.nodes
            for p in range(el.numberofnodes):
                p_glob = nodes[p].numberinMesh -1
                for q in range(el.numberofnodes):
                    q_glob = nodes[q].numberinMesh - 1 
                    ###Assemble global mass matrice
                    Solver.Mass[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Mloc[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Mass[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Mloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    
                    ###Assemble global stiffness matrice
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p,Solver.DOFs*q]
                    Kglob_react[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Kloc_react[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Stiffness[Solver.DOFs*p_glob,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p,Solver.DOFs*q+1]
                    Kglob_react[Solver.DOFs*p_glob,Solver.DOFs*q_glob+1] += Kloc_react[Solver.DOFs*p,Solver.DOFs*q+1]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Kglob_react[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Kloc_react[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Solver.Stiffness[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob] += Kloc[Solver.DOFs*p+1,Solver.DOFs*q]
                    Kglob_react[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob] += Kloc_react[Solver.DOFs*p+1,Solver.DOFs*q]
                 
                    ###Assemble global jacobian matrice
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p,Solver.DOFs*q]
                    Solver.Jac[Solver.DOFs*p_glob,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p,Solver.DOFs*q+1]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob+1] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q+1]
                    Solver.Jac[Solver.DOFs*p_glob+1,Solver.DOFs*q_glob] += Jloc[Solver.DOFs*p+1,Solver.DOFs*q]
                    
                ###Assemble global force matrice
                Solver.Force_bulk[Solver.DOFs*p_glob] += Floc[Solver.DOFs*p]
                Solver.Force_bulk[Solver.DOFs*p_glob+1] += Floc[Solver.DOFs*p+1]
         
        ###Initialize Force vector to bulk force vector (force vector will be modified later if Neumann/Robin BCs are prescribed)
        Solver.Force[:] = Solver.Force_bulk[:]
      
        ###Lump mass matrix if needed
        if Solver.masslumping == True:
            for j in range(Solver.DOFs*Solver.parent.snowpack.mesh.numberofNodes):
                x = sum(Solver.Mass[j,:])
                Solver.Mass[j,:] = 0
                Solver.Mass[j,j] = x
        ### Lump reactive stiffness matrix anyway        
        for j in range(0, Solver.parent.snowpack.mesh.numberofNodes):            
                if j==0:
                    ## For lines corresponding to vapor Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j+2]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+2] = 0
                    ## For lines corresponding to vapor Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j+3]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+3] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+2]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+2] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+3]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+3] = 0
                elif j ==  Solver.parent.snowpack.mesh.numberofNodes -1:
                    ## For lines corresponding to vapor Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j-2]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j-2] = 0
                    ## For lines corresponding to vapor Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j-1]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j-1] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-2]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-2] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-1]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-1] = 0
                else:
                    ## For lines corresponding to vapor Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j-2] + Kglob_react[Solver.DOFs*j, Solver.DOFs*j+2]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j-2] = 0
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+2] = 0
                    ## For lines corresponding to vapor Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j, Solver.DOFs*j-1] + Kglob_react[Solver.DOFs*j, Solver.DOFs*j+3]
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j-1] = 0
                    Kglob_react[Solver.DOFs*j, Solver.DOFs*j+3] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to vapor terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-2] + Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+2]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-2] = 0
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+2] = 0
                    ## For lines corresponding to Temperature Eq. and columns corresponding to temperature terms
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+1] += Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-1] + Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+3]
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j-1] = 0
                    Kglob_react[Solver.DOFs*j+1, Solver.DOFs*j+3] = 0

        return Kglob_react
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

            
            ### For heat
            if Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Adiabatic' or Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Given Temperature':
                pass ### If BC at considered node is prescribed as Adiabatic (natural) or Dirichlet we don't do anything and go see BC on Rhov
            elif  Solver.parent.snowpack.BC_tag[tag].BCType_Heat == 'Heat Flux':
                Flux = Solver.parent.snowpack.BC_tag[tag].get_BC_heat(node_number) ### get the flux from BC section of sif  
                print('Heat Flux =', Flux)
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(xnode)
                ###Modify the Force vector to add the Neumann Flux
                Solver.Force[Solver.DOFs*(node_number-1)+1] += Flux                                
            else:
                print('No boundary condition given for heat on BC ' + tag + '. Setting it to adiabatic.')
        
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ END FUNCTION ASSEMBLING NEUMANN AND/OR ROBIN BOUNDARY CONDITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~   FUNCTION TRANSFORMING THE SYSTEM MdUdt + KU = F in AU = B  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    def AssembleTimeStepping(Solver, Kglob_react):
        if Solver.parent.timestepping_method == "theta":
            Solver.LHS_matrix = (Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * (Solver.Stiffness + Kglob_react)) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * Solver.Jac
            Solver.RHS_vector = np.dot((Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * (Solver.parent.theta-1) * (Solver.Stiffness + Kglob_react)), Solver.solution_prev) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.Force 
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
    def GetResidual(Solver, Kglob_react):
        if Solver.parent.timestepping_method == "theta":
            Left_matrix = Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * Solver.parent.theta * (Solver.Stiffness + Kglob_react) 
            Right_vector = np.dot((Solver.Mass + (Solver.execsolver*Solver.parent.timestepsize) * (Solver.parent.theta-1) * (Solver.Stiffness + Kglob_react)), Solver.solution_prev) + (Solver.execsolver*Solver.parent.timestepsize) * Solver.Force_bulk 
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
        
        Kglob_react = AssembleBulkMatrix(Solver)
        
        AssembleNeumannRobinBC(Solver)    
      
        #Construct AU=B system
        AssembleTimeStepping(Solver, Kglob_react)

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

        nonlin_converged = NormChange <= nonlin_thres
        print('ITERATION N°:', nonlin_count, 'NORM CHANGE =', NormChange)
        nonlin_count += 1
####~~~~~~~~~~~~~~~~~~~~~~~~~ End of Non-linear iteration loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####

#### Update the solution fields for next timestep 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.value[i] = Solver.solution[Solver.DOFs*i] 
        T.value[i] = Solver.solution[Solver.DOFs*i+1] 

#### Calculate the residual
    GetResidual(Solver, Kglob_react)  
#### Update the residual of solution fields for next timestep
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        Rhov.residual[i] = Solver.residual[Solver.DOFs*i] 
        T.residual[i] = Solver.residual[Solver.DOFs*i+1] 

###~~~~~~ Diagnostic of c from Eq. on rhov: M_CC * C = M_CP * (P-P_prev) + K_CP * P - Flux_Rhov ~~~~~~###
    M_CC = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
    M_CP = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
    K_CP = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
    for el in Solver.parent.snowpack.mesh.elements:
        el_size = el.get_elsize()
        Mccloc = np.zeros((el.numberofnodes,el.numberofnodes))
        Mcploc = np.zeros((el.numberofnodes,el.numberofnodes))
        Kloc = np.zeros((el.numberofnodes,el.numberofnodes))
        
        for IP in range(el.numberofIP):
            IP_weight = el.get_IPweight(IP)
            basis, dbasis, ddbasis = el.get_bases(IP)
            IPpos = el.get_IPpos(IP)
                
            ### Get all material properties necessary to calculate Stiffness and force matrices at IP
            Deff = Solver.parent.snowpack.material_perm['Deff'].get_material(IPpos)
            Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)

            for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions)
                    ###~~~~~~LOCAL MASS MATRIX~~~~~~###
                    Mccloc[p,q] +=  - basis[p] * basis[q] * IP_weight * el_size    
                    Mcploc[p,q] +=  ((1-Phii)/(Solver.execsolver*Solver.parent.timestepsize)) * basis[p] * basis[q] * IP_weight * el_size  
                    ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                    Kloc[p,q] +=  ( Deff * dbasis[q] * dbasis[p] ) *  IP_weight * el_size ##K on current P
        
        nodes = el.nodes
        for p in range(el.numberofnodes):
            p_glob = nodes[p].numberinMesh -1
            for q in range(el.numberofnodes):
                q_glob = nodes[q].numberinMesh - 1           
                ###Assemble global mass matrice on unknown C
                M_CC[p_glob,q_glob] += Mccloc[p,q]      
                ###Assemble global mass matrice on  Delta P
                M_CP[p_glob,q_glob] += Mcploc[p,q] 
                ###Assemble global stiffness matrice on current P
                K_CP[p_glob,q_glob] += Kloc[p,q]    
                    
    ###To have a system that is consistant with the main PDE system, M_CP must be lumped whenever M of main system is lumped
    if Solver.masslumping == True:
        for j in range(Solver.parent.snowpack.mesh.numberofNodes):
            x = sum(M_CP[j,:])
            M_CP[j,:] = 0
            M_CP[j,j] = x
    ###To remove oscillations, the MCC matrix must be lumped anyway
    for j in range(Solver.parent.snowpack.mesh.numberofNodes):
        x = sum(M_CC[j,:])
        M_CC[j,:] = 0
        M_CC[j,j] = x
    
    #### Solve the obtained system
    RHS_Vector = np.dot(M_CP, (Rhov.value[:]-Rhov.value_prev_tsp[:])) + np.dot(K_CP, Rhov.value[:]) - Rhov.residual[:]
    Solution_c_diag = np.linalg.solve(M_CC, RHS_Vector)
    
    #### Update the deposition rate 
    for i in range(Solver.parent.snowpack.mesh.numberofNodes): 
        c.value[i] = Solution_c_diag[i]
        c.value_prev_it[i] = Solution_c_diag[i] 
   
