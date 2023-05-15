import FEM_classes
import numpy as np
import copy

###This solver returns the energy contained in each elements at time of execution

def EnergyConservation(Solver):           
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####             
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####    
    ### The solution fields of this solver:
    Energy = Solver.parent.snowpack.field_perm['Energy']  
    ### BEFORE ANY MODIFICATION, store the fields in prev_tsp for the method get_value_prev_tsp of class field_t
    Energy.value_prev_tsp[:] = Energy.value[:]
    
    ### Get needed constant
    Lm = Solver.parent.snowpack.constant_perm['Lm'].get_cst()
    Tfus = Solver.parent.snowpack.constant_perm['Tfus'].get_cst() ### Fusion Temperature (K)

    if ('CoupledCalonneHeatVap' in Solver.parent.solver_perm and 'HansenHeatVap' not in Solver.parent.solver_perm) or ('CalonneHeat' in Solver.parent.solver_perm and 'CalonneVap' in Solver.parent.solver_perm and 'HansenHeatVap' not in Solver.parent.solver_perm) or ('HansenHeatVap' in Solver.parent.solver_perm and Solver.parent.solver_perm['HansenHeatVap'].file == 'Solvers.HansenHeatVap_Solver_1DOF_TForm'): ### When Heat and Vap come from the Calonne system
        ### The fields required for computation:
        Rhov = Solver.parent.snowpack.field_perm['Rhov']
        T = Solver.parent.snowpack.field_perm['Temperature']
        ### Loop on the elements
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()     
            ### Use Gauss method to do the required integrations
            Energy_element_Sensible = 0
            Energy_element_Latent = 0
            for IP in range(el.numberofIP):
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)   
                ### Get all material properties necessary at IP
                rhoCeff = Solver.parent.snowpack.material_perm['rhoCeff'].get_material(IPpos)        
                Phii = Solver.parent.snowpack.field_perm['Phii'].get_value(IPpos)
            
                ### gaussian integration of energy over considered element from basis functions
                for p in range(el.numberofnodes):
                    Energy_element_Sensible += (rhoCeff * (T.value[el.nodes[p].numberinMesh - 1]-Tfus) )* basis[p] * IP_weight * el_size
                    Energy_element_Latent += ( Lm * (1 - Phii) * Rhov.value[el.nodes[p].numberinMesh - 1]  )* basis[p] * IP_weight * el_size  
            #### Update the value of energy contained by considered element (J/m^2)    
            Energy.value[el.numberinMesh -1] = Energy_element_Sensible + Energy_element_Latent 
            
    elif 'HansenHeatVap' in Solver.parent.solver_perm and Solver.parent.solver_perm['HansenHeatVap'].file == 'Solvers.HansenHeatVap_Solver_2DOFs_MixedForm' and 'CoupledCalonneHeatVap' not in Solver.parent.solver_perm: ### When Heat comes from the Hansen system (and Rhov = Rhov_sat)
        ### The fields required for computation:
        H = Solver.parent.snowpack.field_perm['Enthalpy']
        ### Loop on the elements
        for el in Solver.parent.snowpack.mesh.elements:
            el_size = el.get_elsize()     
            ### Use Gauss method to do the required integrations
            Energy_element = 0
            for IP in range(el.numberofIP):
                IP_weight = el.get_IPweight(IP)
                basis, dbasis, ddbasis = el.get_bases(IP)
                IPpos = el.get_IPpos(IP)   
                ### gaussian integration of energy over considered element from basis functions
                for p in range(el.numberofnodes):
                    Energy_element += H.value[el.nodes[p].numberinMesh - 1] * basis[p] * IP_weight * el_size
            #### Update the value of energy contained by considered element (J/m^2)    
            Energy.value[el.numberinMesh -1] = Energy_element  
    else:
        raise TypeError('There is a conflict of solver for the problem of heat and vapor diffusion. Choose only one solver among the HansenHeatVap OR CalonneHeatVap solvers.')
    
    ### DON'T FORGET to update solver variables at previous it for next timestep
    Energy.value_prev_it = copy.deepcopy(Energy.value)

        
        
            
        
                 