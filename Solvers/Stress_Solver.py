import FEM_classes
import numpy as np
import copy

###This solver calculates the value of stress at EACH NODE of the mesh from the mass of overlying element

def Stress_Solver(Solver):           
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~ MAIN SOLVER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####             
            
####~~~~~~~~~~~~~~~~~~~~~~~~~ Initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~####    
    ### The solution field of this solver:
    Sigma = Solver.parent.snowpack.field_perm['Sigma'] 
    ### BEFORE ANY MODIFICATION, store the field in prev_tsp for the method get_value_prev_tsp of the class field_t
    Sigma.value_prev_tsp[:] = Sigma.value[:]

    ### The fields required for computation:
    Phii = Solver.parent.snowpack.field_perm['Phii']   
    ### Get needed constant
    g = Solver.parent.snowpack.constant_perm['g'].get_cst()
    rhoi = Solver.parent.snowpack.constant_perm['rhoi'].get_cst() 
    
    Mass_element=[] ###to store the mass per element    
    ### Loop on the elements
    for el in Solver.parent.snowpack.mesh.elements:
#       print('\n element n°:', el.numberinMesh)
        el_size = el.get_elsize()
        Mass = Phii.value[el.numberinMesh-1]* rhoi * g * el_size
        Mass_element.append(Mass)

    if Solver.parent.snowpack.mesh.elements[0].nodes[0].pos == 0: ##The elements and nodes of the mesh are numbered from bottom to top
        tmp = np.cumsum(Mass_element[::-1]) ### cumulative sum of the element masses from top to bottom
        Sigma_element = tmp[::-1] ####  Cumulative mass associated to the elements which are again numbered from bottom to top
        Sigma.value[:-1] = Sigma_element ### The sigma at every nodes from bottom to the one before top is the cumulative mass of overlying elements
        Sigma.value[-1] = 0 ### The sigma at very top node is zero
    else:
        raise NameError('There is a problem in the structure of the mesh. The first node should be the bottom node located at the soil/snow interface with position 0')
   
    ### No non-linear iteration required in this solver (diagnostic solver) so value_prev_it is the same as value
    Sigma.value_prev_it = copy.deepcopy(Sigma.value)
    
    
#
#    ### Loop on nodes to calculate nodal deposition rate
#    for i in Solver.parent.snowpack.mesh.nodes:
#        xnode = i.pos
#        NumberinMesh = i.numberinMesh
#        # Get needed parameters at considered node
#        s = Solver.parent.snowpack.bodyforce_perm['s'].get_BF(xnode)
#        alpha = Solver.parent.snowpack.bodyforce_perm['alpha'].get_BF(xnode)
#        vkin = Solver.parent.snowpack.bodyforce_perm['vkin'].get_BF(xnode)
#        SWD = Solver.parent.snowpack.material_perm['SWD'].get_material(xnode)
#        
#        #Get Rhov at considered node
#        Rhov = Solver.parent.snowpack.field_perm['Rhov'].value[NumberinMesh-1]
#        
#        ###Update deposition rate at considered node
#        c.value[NumberinMesh-1] = s*alpha*vkin*(Rhov - SWD)
#        
#    ### No non-linear iteration required in this solver (diagnostic solver) so value_prev_it is the same as value
#    c.value_prev_it = copy.deepcopy(c.value)
#        
