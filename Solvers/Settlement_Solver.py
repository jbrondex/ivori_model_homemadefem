import FEM_classes
import numpy as np
import copy
import MatLaw

###This solver contains two dependent modules : 1/ the first one updates the density based on deposition and settlement;
### 2/ the second one moves the mesh to conserve mass per element (and remove the latent energy contains in vapor that leaves each element during compaction )

def Settlement_Solver(Solver):       
               
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####    
    ### The solution fields of this solver:
    Phii = Solver.parent.snowpack.field_perm['Phii']  
    Coordinate = Solver.parent.snowpack.field_perm['Coordinate']
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    Coordinate.value_prev_tsp[:] = Coordinate.value[:]
    Phii.value_prev_tsp[:] = Phii.value[:].copy()

    ### The fields required for computation:
    c = Solver.parent.snowpack.field_perm['DepositionRate']
    Sigma = Solver.parent.snowpack.field_perm['Sigma']
    ### Some fields are also required to close the energy budget (when the snowpack settles some vapor leaves the snowpack with some latent energy)
    Rhov = Solver.parent.snowpack.field_perm['Rhov']
    T = Solver.parent.snowpack.field_perm['Temperature']
    EnergyGone = Solver.parent.snowpack.field_perm['EnergyGoneSettling'] ### Latent energy that leaves the snowpack during settling
    H = Solver.parent.snowpack.field_perm['Enthalpy']  ###Prognostic variable for the Hansen system in mixed-form
    ### Get needed constant
    rhoi = Solver.parent.snowpack.constant_perm['rhoi'].get_cst() 
    ### Some constant are also needed to ensure energy conservation
    Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst()
    ### Get material parameter
    IsCompaction = Solver.parent.snowpack.material_perm['IsCompaction'].function 
    IsDeposition = Solver.parent.snowpack.material_perm['IsDeposition'].function
    ###The time elapsed since last execution of the solver
    Delta_t = Solver.execsolver*Solver.parent.timestepsize

###############################################################################################################
####///////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####~~~~~~~~~~~~~~~~~~~~~~~~~~    MODULE THAT UPDATE OF PHII   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////####
###############################################################################################################  
       
    Deposition_element=[] ### to store the full mass deposited per element
    Deformation_element = [] ### to store the deformation of each element over Delta_t         
    Deposited_volumes = np.zeros(Solver.parent.snowpack.mesh.numberofElements)
    
    ### Loop on the elements
    for k,el in enumerate(Solver.parent.snowpack.mesh.elements):
        el_size = el.get_elsize()
        Deformation = 0 
        Deposition = 0         
        ### Use Gauss method to do the required integrations
        for IP in range(el.numberofIP):
            IP_weight = el.get_IPweight(IP)
            basis, dbasis, ddbasis = el.get_bases(IP)
            IPpos = el.get_IPpos(IP)    
            ### Get the full instantaneous deformation of the element by integrating deformation rate using gaussian integration with value of stresses deduced at IPs through basis function
            if IsCompaction:
                # Compute stress at IP
                SigmaIP = 0
                for p in range(el.numberofnodes):
                    SigmaIP += basis[p] * Sigma.value[el.nodes[p].numberinMesh - 1]
                # Calculate total deformation of the element as (l1-l0/l0) where l1-l0 is given as the integral of local deformation 
                Deformation += - ((SigmaIP / Solver.parent.snowpack.material_perm['viscosity'].get_material(IPpos))* el_size * IP_weight)/el_size
            else:
                pass
            
            ### Get the full deposition within the element by gaussian integration with value of c deduced at IPs through basis function
            if IsDeposition:
                for p in range(el.numberofnodes):
                    Deposition += c.value[el.nodes[p].numberinMesh - 1] * basis[p] * IP_weight ### BE CAREFUL: This should be the averaged deposition within element not the total mass deposited (don't multiply by elt_size !)
            else:
                pass
       
        Deformation_element.append(Deformation) 
        Deposition_element.append(Deposition)
        Deposited_volumes[k] = el_size * Deposition/rhoi

 ####~~~~~~~~~~~~~~~~~~~~~~~~~ Update Phii ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####     
    Phii.value[:] = (Phii.value_prev_tsp[:] + Delta_t* (np.array(Deposition_element[:])/(rhoi))) / (1 + Delta_t * np.array(Deformation_element[:]))
    ### DON'T FORGET to update solver variables at previous it for next timestep
    Phii.value_prev_it = copy.deepcopy(Phii.value)    
    
    
###############################################################################################################
####///////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\####
####~~~~~~~~~~~~~~~~~~~~~~~~~~     MESH DEFORMATION MODULE     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////////////////////////////////####
###############################################################################################################  
    
####~~~~~~~~~~~~~~~~~~~~~~~~~ Move the mesh nodes if compaction is activated ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~ And calculate latent energy gone to close energ budget ~~~~~~~~~~~~~~~~~~~~~~~~~~~####
    if IsCompaction:
        node_newpos=[]
        znode = 0
        node_newpos.append(znode) ### The bottom node does not move, it stays at the soil/snow interface
        
        for k,el in enumerate(Solver.parent.snowpack.mesh.elements):
            old_el_size = el.get_elsize()### We need to wait that new element size has been calculated for all elements before updating positions of nodes...
            new_el_size = old_el_size * (1 + Delta_t * Deformation_element[el.numberinMesh - 1])### ... otherwise the get_elsize will not return the old size of the element...
            znode += new_el_size
            node_newpos.append(znode)

            #### Below we compute the latent energy leaving the snowpack through settling
            Rhov_moy = 0.5 * (Rhov.value[el.nodes[0].numberinMesh - 1] + Rhov.value[el.nodes[1].numberinMesh - 1])
            EnergyGone.value[k] = Rhov_moy * Lm * (old_el_size - new_el_size)
        ### Update position of nodes and variable coordinate
        for node in Solver.parent.snowpack.mesh.nodes:
            node.update_pos(node_newpos[node.numberinMesh-1])
            Coordinate.value[node.numberinMesh-1] = node.pos   
        ### DON'T FORGET to update solver variables at previous it for next timestep
        Coordinate.value_prev_it = copy.deepcopy(Coordinate.value)        
                       
    ### Update Enthalpy field for the Hansen Solver in mixed form
    if 'HansenHeatVap' in Solver.parent.solver_perm and 'CoupledCalonneHeatVap' not in Solver.parent.solver_perm: ### ...when Heat comes from the Hansen system, it is more complicated as H is a nodal variable of the Hansen 2DOFs system
        ### it requires to calculate the H corresponding to old T field but new Phii
        H_new = np.zeros((Solver.parent.snowpack.mesh.numberofNodes))
        M_HH = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
        K_HT = np.zeros((Solver.parent.snowpack.mesh.numberofNodes, Solver.parent.snowpack.mesh.numberofNodes))
        F_H = np.zeros((Solver.parent.snowpack.mesh.numberofNodes))
        ## Get the needed constants only once
        Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst() ### Massic latent heat of sublimation/deposition        
        Tfus = Solver.parent.snowpack.constant_perm['Tfus'].get_cst() ### Fusion Temperature (K)
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize() ### That is the new elt_size
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
                T_prev= Solver.parent.snowpack.field_perm['Temperature'].get_value_prev_it(IPpos)
                dSWDdT = Solver.parent.snowpack.material_perm['dSWDdT'].get_material(IPpos)
                for p in range(el.numberofnodes): ###loop on matrice lines (correspond to test functions )
                    for q in range(el.numberofnodes): ###loop on matrices columns (correspond to trial functions )              
                        Mhloc[p,q] += basis[p] * basis[q] *  IP_weight * el_size ##M_{HH}
                        ###~~~~~~LOCAL STIFFNESS MATRIX~~~~~~###
                        Khtloc[p,q] += - ((1 - Phii) * Lm * dSWDdT + rhoCeff) * basis[p] * basis[q] *  IP_weight * el_size ##K_{HT}
                        ###~~~~~~LOCAL FORCE VECTOR~~~~~~###
                    Fhloc[p] += ( (1 - Phii) * Lm * (SWD - dSWDdT * T_prev) - rhoCeff * Tfus)  * basis[p] * IP_weight * el_size ##F_{H}
    
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
                
        ###If mass matrix is lumped in the resolution of the Hansen system PDE, then it must also be lumped here for consistency
        if Solver.parent.solver_perm['HansenHeatVap'].masslumping == True:
            for j in range(Solver.parent.snowpack.mesh.numberofNodes):
                x = sum(M_HH[j,:])
                M_HH[j,:] = 0
                M_HH[j,j] = x
        ###Now we compute the Hnew corresponding to new energy profile obtained from new Phii and old T field
        RHS_Vector = F_H - np.dot(K_HT, T.value[:])
        H_new[:] = np.linalg.solve(M_HH, RHS_Vector)
        H.value[:] = H_new[:] ### H is a prognostic field and must be updated to account for the fact that some energy has left the system during compaction
        
            

        
        

        
        
            
        
                 
