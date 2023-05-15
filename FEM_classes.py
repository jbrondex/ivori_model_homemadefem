import Gauss_Integ
import MatLaw
import USF
import Ref_Elements
import numpy as np


def Norm(solution):
    return np.sqrt(np.sum(np.power(solution,2)))

class model_t:
    def __init__(self, Simu, Snowpack, Cst, Mat, BF, BC1, Solver1):
        
        ##Predefine attributes of class model_t that are set in sif file 
        self.name = Simu.Name
        self.numericalmethod = 'FEM' ##We are in FEM_classes so FEM by default but to improve later
        self.save_csv = Simu.Save_csv
        self.outputpath = Simu.Output_Path
        self.type = None
        self.timestepsize = None
        self.totaltimesteps = None
        self.outputinterval = None
        self.timestepping_method = None
        self.theta = None
       
        ##Initialize attributes that depends on simulation type
        if Simu.Type == 'Transient':
            self.type = Simu.Type
            self.timestepsize = Simu.Timestep_Size
            self.totaltimesteps = Simu.Total_TimeSteps
            self.outputinterval = Simu.Output_Interval
            
            if Simu.Timestepping_Method == 'theta':
                self.timestepping_method = Simu.Timestepping_Method
                self.theta = Simu.Theta
            elif Simu.Timestepping_Method == 'BDF':
                self.timestepping_method = Simu.Timestepping_Method
            else:
                raise NameError('TimeStepping method not valid')
        elif Simu.Type == 'Steady':
            self.type = Simu.Type
        else:
            raise NameError('Simulation Type not valid')
        
        ##Initialize attributes of class model_t that are not set in sif file 
        self.time = 0
        self.timestep = 1
        
        ##Create inner classes
        self.snowpack = snowpack_t(Snowpack, Cst, Mat, BF, BC1, self)
        
        self.solver_list=[]
        self.solver_perm={}
        for solver in Solver1.Solver_list: ###/!\ HERE THERE IS SOME HARD CODE (what happens if first instance of class Solver of sif is not called Solver1 ?) ! To improve later
            self.solver_list.append(solver_t(solver, self))
            self.solver_perm[solver.Name] = solver_t(solver, self)
              
    ##Define methods of the class 
    def update_timestep(self):
        self.timestep += 1
        
    def update_time(self):
        self.time += self.timestepsize
               
class snowpack_t: 
    def __init__(self, Snowpack, Cst, Mat, BF, BC1, parent):
        self.parent = parent
        self.mesh = mesh_t(Snowpack)
        self.field_list = []
        self.field_perm = {}
        
        # Hard coded fields that are always present with the snow models
        # 0: Temperatture
        # 1: Phii
        # 2: Rhov
        # 3: Coordinate
        list_field_always_present = ['Temperature', 'Phii', 'Rhov', 'Enthalpy', 'DepositionRate', 'Sigma', 'Energy', 'EnergyGoneSettling','Coordinate']
        list_field_type = ['Nodal', 'Elemental', 'Nodal', 'Nodal', 'Nodal', 'Nodal', 'Elemental', 'Elemental', 'Nodal']
        for i,field_name in enumerate(list_field_always_present):
            # Create and initalize field
            self.field_list.append(field_t(field_name, list_field_type[i], self))
            self.field_perm[field_name] = field_t(field_name, list_field_type[i], self)
        
        # Constants obtained from reading Cst section (i.e. class for now) of sif file
        self.constant_list=[]
        self.constant_perm={}
        for constant_name in vars(Cst): ##loop on list of constants defined in the sif
            # Create and initialize material property
            self.constant_list.append(constant_t(constant_name, eval('Cst.'+constant_name), self))
            self.constant_perm[constant_name] = constant_t(constant_name, eval('Cst.'+constant_name), self)
                
        # Materials obtained from reading Mat section (i.e. class for now) of sif file
        self.material_list = []
        self.material_perm = {}
        for material_name in vars(Mat): ##loop on list of material properties defined in the sif
            # Create and initialize material property
            self.material_list.append(material_t(material_name, eval('Mat.'+material_name), self))
            self.material_perm[material_name] = material_t(material_name, eval('Mat.'+material_name), self)
        
        # Body Forces obtained from reading BF section (i.e. class for now) of sif file
        self.bodyforce_list = []
        self.bodyforce_perm = {}
        for bodyforce_name in vars(BF): ##loop on list of material properties defined in the sif
            # Create and initialize material property
            self.bodyforce_list.append(bodyforce_t(bodyforce_name, eval('BF.'+bodyforce_name), self))
            self.bodyforce_perm[bodyforce_name] = bodyforce_t(bodyforce_name, eval('BF.'+bodyforce_name), self)
        
        # Boundary Conditions obtained from reading BC section (i.e. class for now) of sif file
        self.BC_list = []
        self.BC_tag = {}
        for BC in BC1.BC_list: ##loop on list of material properties defined in the sif
            self.BC_list.append(boundarycondition_t(BC, self))
            self.BC_tag[BC.Tag] = boundarycondition_t(BC, self)

class mesh_t:

    def __init__(self, Snowpack):
        
        #predefine attributes of mesh class
        self.nodes = []
        self.elements = []
        self.boundaries = []
        self.numberofNodes = 0
        self.numberofElements = 0
        self.numberofBoundaryNodes = 0        
        if Snowpack.Auto_Mesh:
            self.numberofNodes = Snowpack.Nnodes
            for i in range(self.numberofNodes):
                self.nodes.append(node_t(i+1,i*Snowpack.Htot/(self.numberofNodes-1))) ##First node of mesh is always bottom node as coded here
            self.numberofElements = self.numberofNodes - 1
            for i in range(self.numberofElements):
                n1, n2 = self.nodes[i], self.nodes[i+1]
                elem = element_t(n1,n2,i+1) ##First element of mesh is always bottom element as coded here
                self.elements.append(elem)
            self.numberofBoundaryNodes = 2 ### Normally always the case in 1D but might have to be improved later on 
            self.boundaries.append(['bot', self.nodes[0]])### Bottom node is first node by default but might have to be improved later on 
            self.boundaries.append(['top', self.nodes[self.numberofNodes-1]])### Top node is last node by default but might have to be improved later on
        else:
            meshfile = Snowpack.MeshFile
            with open(meshfile, 'r') as f:
                f.readline() #Skip first line
                self.numberofNodes = int(f.readline())
                for i in range(self.numberofNodes):
                    self.nodes.append(node_t(i+1, float(f.readline())))
                if self.nodes[0].pos != 0:
                    raise NameError('Error in mesh file. First node must correspond to soil/snow interface and must have position 0')                    
                f.readline()
                f.readline()
                self.numberofElements = int(f.readline())
                for i in range(self.numberofElements):
                    l = f.readline().split()
                    n1, n2 = self.nodes[int(l[1])-1], self.nodes[int(l[2])-1]
                    elem = element_t(n1,n2,i+1)
                    self.elements.append(elem)
#                
                f.readline()
                f.readline()
                self.numberofBoundaryNodes = int(f.readline())
                for i in range(self.numberofBoundaryNodes):
                    l = f.readline().split()
                    self.boundaries.append([l[0], self.nodes[int(l[1])-1]])
                
class element_t:

    def __init__(self, n1, n2, number, nIP=2):
        self.nodes = [n1,n2]
        self.numberinMesh = number
        self.numberofnodes = 2
        self.type = 1
        self.numberofIP = nIP
        
               
    def get_elsize(self):
        return abs(self.nodes[0].pos - self.nodes[1].pos)
        
    def get_bases(self,IP):
        L= self.get_elsize()
        for i in range(self.numberofnodes):
            posIP = Gauss_Integ.IPpos[self.numberofIP-1][IP]
            bases_values = Ref_Elements.evaluate_basis_ref(posIP)
            dbases_values = Ref_Elements.evaluate_dbasis_ref(posIP) * 2/L
            ddbases_values = Ref_Elements.evaluate_ddbasis_ref(posIP) * 4/(L**2)
        
        return bases_values, dbases_values, ddbases_values
        
    def get_IPweight(self,IP):
        return Gauss_Integ.IPweights[self.numberofIP-1][IP]
    
    def get_IPpos(self,IP):
        posIP = Gauss_Integ.IPpos[self.numberofIP-1][IP] ##position of IP in reference element
        xIP = 1/2*(1-posIP)*self.nodes[0].pos + 1/2*(1+posIP)*self.nodes[1].pos
        return xIP
        

class node_t:

    def __init__(self,n,xpos):
        self.numberinMesh = n
        self.pos = xpos   
        
    def update_pos(self, new_pos):
        self.pos = new_pos
                                             
class field_t:
    def __init__(self, name, NodEl, parent):
        self.parent = parent
        self.name = name
        self.NodEl = NodEl ### Nodal or Elemental value
        self.value = None
        self.value_prev_it = None
        self.value_prev_tsp = None 
        self.residual = None ### Necessary for the diagnostic of deposition rate Calonne 2DOFs model (and also later for couplings to external models)
        
        if self.NodEl == 'Nodal':
            self.value = np.zeros(self.parent.mesh.numberofNodes)
            self.value_prev_it = np.zeros(self.parent.mesh.numberofNodes)
            self.value_prev_tsp = np.zeros(self.parent.mesh.numberofNodes)
            self.residual = np.zeros(self.parent.mesh.numberofNodes)
        elif self.NodEl == 'Elemental':
            self.value = np.zeros(self.parent.mesh.numberofElements)
            self.value_prev_it = np.zeros(self.parent.mesh.numberofElements)
            self.value_prev_tsp = np.zeros(self.parent.mesh.numberofElements)
            self.residual = np.zeros(self.parent.mesh.numberofElements)
            
        if self.name == 'Coordinate':
            for i in range(self.parent.mesh.numberofNodes):
                self.value[i] = self.parent.mesh.nodes[i].pos
                self.value_prev_it[i] =  self.parent.mesh.nodes[i].pos
                self.value_prev_tsp[i] =  self.parent.mesh.nodes[i].pos
                
    def update_value(self, new_value):
        self.value = new_value 
    
    def update_value_prev_it(self, new_value):
        self.value_prev_it = new_value
    
    def get_value(self, xpos):
        ###create a list with position of all nodes of the mesh
        xpos_list=[]
        for i in range(self.parent.mesh.numberofNodes):
            xpos_list.append(self.parent.mesh.nodes[i].pos)
        
        if self.NodEl == 'Nodal': ###Nodal field
            if xpos in xpos_list: ### We are looking exactly at a node
                value = self.value[xpos_list.index(xpos)]
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                value = np.interp(xpos, xpos_list, self.value) ### ONLY LINEAR INTERPOLATION. TO BE ADAPTED FOR HIGHER ORDER ELEMENTS
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated out of domain')
        elif self.NodEl == 'Elemental': ###Elemental field
            if xpos in xpos_list: ### We are exactly at a node
                if xpos == 0: ### We are at bottom node
                    value = self.value[xpos_list.index(xpos)] ### We take the value applying to bottom element
                elif xpos == xpos_list[-1]: ### We are at top node
                    value = self.value[xpos_list.index(xpos)-1] ### We take the value applying to top element
                else: ### We are exactly at a node that is neither the bottom node nor the top node
                    value = 0.5 * (self.value[xpos_list.index(xpos)-1] + self.value[xpos_list.index(xpos)]) ###If we are exactly between two elements (i.e. exactly at a node) returned value is arithmetic mean between values at elements just above and just below
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                idx = (np.abs(np.array(xpos_list)-xpos)).argmin() ### We need to know where we are from positions of nodes
                if xpos < xpos_list[idx]: ### idx is index of node above xpos so we are in element of index idx-1 
                    value = self.value[idx-1] ### returned value is the one of the element
                else: ### idx is index of node below xpos so we are in element of index idx
                    value = self.value[idx]
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated out of domain')
        return value

    def get_value_prev_it(self, xpos):
        ###create a list with position of all nodes of the mesh (the mesh does not change between two non-linear iterations)
        xpos_list=[]
        for i in range(self.parent.mesh.numberofNodes):
            xpos_list.append(self.parent.mesh.nodes[i].pos)
        
        if self.NodEl == 'Nodal': ###Nodal field
            if xpos in xpos_list: ### We are looking exactly at a node
                value = self.value_prev_it[xpos_list.index(xpos)]
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                value = np.interp(xpos, xpos_list, self.value_prev_it) ### ONLY LINEAR INTERPOLATION. TO BE ADAPTED FOR HIGHER ORDER ELEMENTS
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated out of domain')
        elif self.NodEl == 'Elemental': ###Elemental field
            if xpos in xpos_list: ### We are exactly at a node
                if xpos == 0: ### We are at bottom node
                    value = self.value_prev_it[xpos_list.index(xpos)] ### We take the value applying to bottom element
                elif xpos == xpos_list[-1]: ### We are at top node
                    value = self.value_prev_it[xpos_list.index(xpos)-1] ### We take the value applying to top element
                else: ### We are exactly at a node that is neither the bottom node nor the top node
                    value = 0.5 * (self.value_prev_it[xpos_list.index(xpos)-1] + self.value_prev_it[xpos_list.index(xpos)]) ###If we are exactly between two elements (i.e. exactly at a node) returned value is arithmetic mean between values at elements just above and just below
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                idx = (np.abs(np.array(xpos_list)-xpos)).argmin() ### We need to know where we are from positions of nodes
                if xpos < xpos_list[idx]: ### idx is index of node above xpos so we are in element of index idx-1 
                    value = self.value_prev_it[idx-1] ### returned value is the one of the element
                else: ### idx is index of node below xpos so we are in element of index idx
                    value = self.value_prev_it[idx]
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated out of domain')
        return value
    
    def get_value_prev_tsp(self, xpos):
        ###create a list with position of all nodes of the mesh (/!\ the mesh might have changed between current and previous timestep)
        xpos_list=[]
        for i in range(len(self.parent.field_perm['Coordinate'].value_prev_tsp)):
            xpos_list.append(self.parent.field_perm['Coordinate'].value_prev_tsp[i])
        
        if self.NodEl == 'Nodal': ###Nodal field
            if xpos in xpos_list: ### We are looking exactly at a node
                value = self.value_prev_tsp[xpos_list.index(xpos)]
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                value = np.interp(xpos, xpos_list, self.value_prev_tsp) ### ONLY LINEAR INTERPOLATION. TO BE ADAPTED FOR HIGHER ORDER ELEMENTS
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated for previous timestep is out of domain as defined in previous timestep')
        elif self.NodEl == 'Elemental': ###Elemental field
            if xpos in xpos_list: ### We are exactly at a node
                if xpos == 0: ### We are at bottom node
                    value = self.value_prev_tsp[xpos_list.index(xpos)] ### We take the value applying to bottom element
                elif xpos == xpos_list[-1]: ### We are at top node
                    value = self.value_prev_tsp[xpos_list.index(xpos)-1] ### We take the value applying to top element
                else: ### We are exactly at a node that is neither the bottom node nor the top node
                    value = 0.5 * (self.value_prev_tsp[xpos_list.index(xpos)-1] + self.value_prev_tsp[xpos_list.index(xpos)]) ###If we are exactly between two elements (i.e. exactly at a node) returned value is arithmetic mean between values at elements just above and just below
            elif min(xpos_list) < xpos < max(xpos_list):                 ### We are in between two nodes
                idx = (np.abs(np.array(xpos_list)-xpos)).argmin() ### We need to know where we are from positions of nodes
                if xpos < xpos_list[idx]: ### idx is index of node above xpos so we are in element of index idx-1 
                    value = self.value_prev_tsp[idx-1] ### returned value is the one of the element
                else: ### idx is index of node below xpos so we are in element of index idx
                    value = self.value_prev_tsp[idx]
            else:  ###we are out of domain
                raise NameError('position at which field ' + self.name + ' is being evaluated for previous timestep is out of domain as defined in previous timestep')
        return value
    
class constant_t:
    def __init__(self, name, value, parent):
        self.parent = parent
        self.name = name
        self.value = value
    
    def get_cst(self):
        return self.value    
    
class material_t:
    def __init__(self, name, function, parent):
        self.parent = parent
        self.name = name
        self.function = function
    
    def get_material(self, pos):
        ## if constant
        if type(self.function) == float or type(self.function) == int:
            value = self.function
        ## if given by a function 
        elif type(self.function) == str:
            function_name = self.function
            liste = function_name.split('_')
            if liste[0] not in ['USF', 'MatLaw']:
                raise NameError('If not constant, function prescribing Material Laws '+ self.name +' is a string that must start by "USF" or "MatLaw"')
            if np.size(liste) < 3:
                raise NameError('If not constant, function prescribing Material Laws '+ self.name +' is a string with formalism "USF/MatLaw_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
            input_fields = liste[2:] ### Get the list of variables that are the input of the function
            if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')
            input_values = []
            for input_field in input_fields:
                if input_field == 'Time':
                    input_value = self.parent.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                else:
                    input_value = self.parent.field_perm[input_field].get_value_prev_it(pos) ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED (Same problem as above)
                input_values.append(input_value)
            value = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
        return value
    
class bodyforce_t:
    def __init__(self, name, function, parent):
        self.parent = parent
        self.name = name
        self.function = function
    
    def get_BF(self, pos):
        ## if constant
        if type(self.function) == float or type(self.function) == int:
            value = self.function
        ## if given by a function 
        elif type(self.function) == str:
            function_name = self.function
            liste = function_name.split('_')
            if liste[0] not in ['USF', 'MatLaw']:
                raise NameError('If not constant, function prescribing Body Force '+ self.name +' is a string that must start by "USF" or "MatLaw"')
            if np.size(liste) < 3:
                raise NameError('If not constant, function prescribing Body Force '+ self.name +' is a string with formalism "USF/MatLaw_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
            input_fields = liste[2:] ### Get the list of variables that are the input of the function
            if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')
            input_values = []
            for input_field in input_fields:
                if input_field == 'Time':
                    input_value = self.parent.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                else:
                    input_value = self.parent.field_perm[input_field].get_value_prev_it(pos) ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED (Same problem as above)
                input_values.append(input_value)
            value = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
        return value

class boundarycondition_t:
    
    def __init__(self, BC, parent):
        self.parent = parent                
        self.tag = BC.Tag
        if self.tag not in ['top', 'bot']:
            raise NameError('Tag of BCs must be set to "top" or "bot"')
        
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~     Init for Heat      ~~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        self.BCType_Heat = None
        if hasattr(BC, 'BCType_Heat'):
            self.BCType_Heat = BC.BCType_Heat
        
        self.functionHeat = None
        if self.BCType_Heat == 'Given Temperature':
            if not hasattr(BC, 'Temperature'):
                raise NameError ('There must be a key word Temperature in the BC '+ BC.Tag + ' Section of sif as key word BCType_Heat is set to '+ self.BCType_Heat )
            else:
                self.functionHeat = BC.Temperature
        elif self.BCType_Heat == 'Heat Flux': 
            if not hasattr(BC, 'HeatFlux'):
                raise NameError ('There must be a key word HeatFlux in the BC '+ BC.Tag + ' Section of sif as key word BCType_Heat is set to '+ self.BCType_Heat )
            else:
                self.functionHeat = BC.HeatFlux
        elif self.BCType_Heat == 'Adiabatic':
            self.functionHeat = 0
        elif self.BCType_Heat == 'Energy Balance': #### HERE ALL PARAMETERS RELATED TO THE ENERGY BALANCE MODEL WILL HAVE TO BE IMPLEMENTED AT SOME POINT
            self.functionHeat = None ### Could be a dictionnary with all parameters of energy balance (sigma, emissivity...) and associated functions
         
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~     Init for Rhov      ~~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###           
        self.BCType_Rhov = None
        if hasattr(BC, 'BCType_Rhov'):
            self.BCType_Rhov = BC.BCType_Rhov
        
        self.functionRhov = None
        if self.BCType_Rhov == 'Given Rhov':
            if not hasattr(BC, 'Rhov'):
                raise NameError ('There must be a key word Rhov in the BC '+ BC.Tag + ' Section of sif as key word BCType_Rhov is set to '+ self.BCType_Rhov )
            else:
                self.functionRhov = BC.Rhov
        elif self.BCType_Rhov == 'Rhov Flux': 
            if not hasattr(BC, 'RhovFlux'):
                raise NameError ('There must be a key word RhovFlux in the BC '+ BC.Tag + ' Section of sif as key word BCType_Rhov is set to '+ self.BCType_Rhov )
            else:
                self.functionRhov = BC.RhovFlux
        elif self.BCType_Rhov == 'No Flux':
            self.functionRhov = 0
        elif self.BCType_Rhov == 'Rhov Balance': #### HERE ALL PARAMETERS RELATED TO THE MASS BALANCE MODEL WILL HAVE TO BE IMPLEMENTED AT SOME POINT
            self.Dummy = None 
            
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        ###~~~     Init for DepositionRateHansen      ~~~###
        ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
        self.BCType_DepositionRate = None
        if hasattr(BC, 'BCType_DepositionRate'):
            self.BCType_DepositionRate = BC.BCType_DepositionRate
        
        self.functionDepositionRate = None
        if self.BCType_DepositionRate == 'Given Deposition Rate':
            if not hasattr(BC, 'DepositionRate'):
                raise NameError ('There must be a key word Depostion Rate in the BC '+ BC.Tag + ' Section of sif as key word BCType_DepositionRate is set to '+ self.BCType_DepositionRate )
            else:
                print('A Dirichlet BC has been prescribed for the Deposition Rate, yet only a Neumann BC is physically meaningfull. Please, double check.')
                self.functionDepositionRate = BC.DepositionRate
        elif self.BCType_DepositionRate == 'Mass Flux': 
            if not hasattr(BC, 'MassFlux'):
                raise NameError ('There must be a key word MassFlux in the BC '+ BC.Tag + ' Section of sif as key word BCType_DepositionRate is set to '+ self.BCType_DepositionRate )
            else:
                self.functionDepositionRate = BC.MassFlux
        elif self.BCType_DepositionRate == 'No Flux':
            self.functionDepositionRate = 0
        elif self.BCType_DepositionRate == 'Mass Balance': #### HERE ALL PARAMETERS RELATED TO THE VAPOR MASS BALANCE MODEL WILL HAVE TO BE IMPLEMENTED AT SOME POINT
            self.functionDepositionRate = None ### Could be a dictionnary with all parameters of energy balance (sigma, emissivity...) and associated functions
        
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~     Method for Heat      ~~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~### 
    def get_BC_heat(self, NodeNumber): ## get value of BC for heat from associated function at considered boundary (i.e. top or bot) . BE CAREFUL NodeNumber start from 1.
        
    ## if constant
        if type(self.functionHeat) == float or type(self.functionHeat) == int:
            value = self.functionHeat
        ## if given by a function 
        elif type(self.functionHeat) == str:
            function_name = self.functionHeat
            liste = function_name.split('_')
            if liste[0] not in ['USF', 'MatLaw']:
                raise NameError('If not constant, functions used to calculate BC on Heat '+ self.tag +' must be string starting with "USF" or "MatLaw"')
            if np.size(liste) < 3:
                raise NameError('If not constant, functions used to calculate BC on Heat '+ self.tag +' must be string with formalism "USF_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
            input_fields = liste[2:] ### Get the list of variables that are the input of the function
            if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')          
            input_values = []
            for input_field in input_fields:
                if input_field == 'Time':
                    input_value = self.parent.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                else:
                    if self.parent.field_perm[input_field].NodEl == 'Nodal':
                        input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ###If the input_field is nodal we simply get its value at previous iteration at the considered BC node
                    elif  self.parent.field_perm[input_field].NodEl == 'Elemental':
                        if self.tag == 'bot':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ##If the input field is elemental and we are at bottom node we take value of first element
                        elif self.tag =='top':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-2] ##If the input field is elemental and we are at top node we take value of last element
                input_values.append(input_value)
            value = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
        return value
        
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~     Method for Rhov      ~~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~### 
    def get_BC_vapor(self, NodeNumber): ## get value of BC for heat from associated function at considered boundary (i.e. top or bot). BE CAREFUL NodeNumber start from 1. 
        
    ## if constant
        if type(self.functionRhov) == float or type(self.functionRhov) == int:
            value = self.functionRhov
        ## if given by a function 
        elif type(self.functionRhov) == str:
            function_name = self.functionRhov
            liste = function_name.split('_')
            if liste[0] not in ['USF', 'MatLaw']:
                raise NameError('If not constant, functions used to calculate BC on Vapor '+ self.tag +' must be string starting with "USF" or "MatLaw"')
            if np.size(liste) < 3:
                raise NameError('If not constant, functions used to calculate BC on Vapor '+ self.tag +' must be string with formalism "USF_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
            input_fields = liste[2:] ### Get the list of variables that are the input of the function
            if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')          
            input_values = []
            for input_field in input_fields:
                if input_field == 'Time':
                    input_value = self.parent.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                else:
                    if self.parent.field_perm[input_field].NodEl == 'Nodal':
                        input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ###If the input_field is nodal we simply get its value at previous iteration at the considered BC node
                    elif  self.parent.field_perm[input_field].NodEl == 'Elemental':
                        if self.tag == 'bot':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ##If the input field is elemental and we are at bottom node we take value of first element
                        elif self.tag == 'top':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-2] ##If the input field is elemental and we are at top node we take value of last element
                input_values.append(input_value)
            value = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
        return value
    
    
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
    ###~~~     Method for Deposition Rate      ~~~###
    ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~### 
    def get_BC_depositionrate(self, NodeNumber): ## get value of BC for deposition rate from associated function at considered boundary (i.e. top or bot) . BE CAREFUL NodeNumber start from 1.
        
    ## if constant
        if type(self.functionDepositionRate) == float or type(self.functionDepositionRate) == int:
            value = self.functionDepositionRate
        ## if given by a function 
        elif type(self.functionDepositionRate) == str:
            function_name = self.functionDepositionRate
            liste = function_name.split('_')
            if liste[0] not in ['USF', 'MatLaw']:
                raise NameError('If not constant, functions used to calculate BC on Deposition Rate '+ self.tag +' must be string starting with "USF" or "MatLaw"')
            if np.size(liste) < 3:
                raise NameError('If not constant, functions used to calculate BC on Deposition Rate '+ self.tag +' must be string with formalism "USF_Name_Var1_Var2_Varn". There must be at least one variable prescribed.')
            input_fields = liste[2:] ### Get the list of variables that are the input of the function
            if not set(input_fields).issubset(['Coordinate', 'Temperature', 'Rhov', 'DepositionRate', 'Sigma', 'Phii', 'Time']): ### accepted are only in this list
                raise NameError('Variable in '+liste[0]+'_'+liste[1]+' not defined.')          
            input_values = []
            for input_field in input_fields:
                if input_field == 'Time':
                    input_value = self.parent.parent.time ###/!\ HERE MIGHT NOT WORK ANYMORE IF ARCHITECTURE IS CHANGED
                else:
                    if self.parent.field_perm[input_field].NodEl == 'Nodal':
                        input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ###If the input_field is nodal we simply get its value at previous iteration at the considered BC node
                    elif  self.parent.field_perm[input_field].NodEl == 'Elemental':
                        if self.tag == 'bot':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-1] ##If the input field is elemental and we are at bottom node we take value of first element
                        elif self.tag =='top':
                            input_value = self.parent.field_perm[input_field].value_prev_it[NodeNumber-2] ##If the input field is elemental and we are at top node we take value of last element
                input_values.append(input_value)
            value = eval(liste[0]+'.'+liste[0]+'_'+liste[1])(*input_values)
        return value
             
class solver_t:
    def __init__(self, Solver, parent):
        self.parent = parent
        
        self.name = Solver.Name 
        self.file = Solver.File
        
        self.execsolver = 1 ### Solver is executed at each timestep by default
        if hasattr(Solver, 'ExecSolver'):
            self.execsolver = Solver.ExecSolver
         
        self.DOFs = 1
        if hasattr(Solver, 'DOFs'):
            self.DOFs = Solver.DOFs
        else:
            print('Solver '+ Solver.Name +' has no integer attribute DOFs. Setting to default value DOFs=1.')
        
        self.masslumping = False #### No mass lumping by default
        if hasattr(Solver, 'MassLumping'):
            self.masslumping = Solver.MassLumping
        else:
            print('Solver '+ Solver.Name +' has no boolean attribute MassLumping. Setting to default value False.')

        
        self.solution = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.solution_current_it = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.solution_prev = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.residual = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        
        self.nonlin_max_it = None
        if hasattr(Solver, 'Nonlin_Max_It'):
            self.nonlin_max_it = Solver.Nonlin_Max_It
        
        self.nonlin_thres = None
        if hasattr(Solver, 'Nonlin_Thres'):
            self.nonlin_thres = Solver.Nonlin_Thres
        
        self.nonlin_newton_after_it = None
        if hasattr(Solver, 'Nonlin_Newton_After_It'):  
            self.nonlin_newton_after_it = Solver.Nonlin_Newton_After_It
        
        self.nonlin_newton_after_thres = None
        if hasattr(Solver, 'Nonlin_Newton_After_Thres'): 
            self.nonlin_newton_after_thres = Solver.Nonlin_Newton_After_Thres
        
        self.relaxation_factor = None
        if hasattr(Solver, 'Relaxation_Factor'):
            self.relaxation_factor = Solver.Relaxation_Factor 
        
        self.Stiffness = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes,self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.Mass = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes,self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.Force = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.Force_bulk = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes)) ##Necessary to compute the residual
        self.Jac = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes,self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.RHS_vector = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes))
        self.LHS_matrix = np.zeros((self.DOFs * self.parent.snowpack.mesh.numberofNodes,self.DOFs * self.parent.snowpack.mesh.numberofNodes))

    def AssembleTimeStepping(self):
        if self.parent.timestepping_method == "theta":
            self.LHS_matrix = (self.Mass + (self.execsolver*self.parent.timestepsize) * self.parent.theta * self.Stiffness) + (self.execsolver*self.parent.timestepsize)*self.parent.theta*self.Jac
            self.RHS_vector = np.dot((self.Mass + (self.execsolver*self.parent.timestepsize)*(self.parent.theta-1)*self.Stiffness), self.solution_prev) + (self.execsolver*self.parent.timestepsize)*self.Force 
            self.RHS_vector += np.dot((self.execsolver*self.parent.timestepsize)*self.parent.theta*self.Jac, self.solution_current_it)
    
    def GetResidual(self):
        if self.parent.timestepping_method == "theta":
            LHS_matrix = (self.Mass + (self.execsolver*self.parent.timestepsize) * self.parent.theta * self.Stiffness)
            RHS_vector = np.dot((self.Mass + (self.execsolver*self.parent.timestepsize)*(self.parent.theta-1)*self.Stiffness), self.solution_prev) + (self.execsolver*self.parent.timestepsize)*self.Force_bulk 
            self.residual = (np.dot(LHS_matrix, self.solution)-RHS_vector) / (self.execsolver*self.parent.timestepsize)
    
    def cleanMatrices(self):
        self.Stiffness[:] = 0
        self.Mass[:] = 0
        self.Jac[:] = 0
        self.Force[:] = 0
        self.Force_bulk[:] = 0
    
    def NormChange(self):
        return (2 *abs(Norm(self.solution) - Norm(self.solution_current_it))) / (Norm(self.solution_current_it) + Norm(self.solution))  ### Same definition as in Elmer
    
    def linSolve(self):
        self.solution = np.linalg.solve(self.LHS_matrix,self.RHS_vector)
                        
    def DirichletBC(self):
        for bcs in self.parent.snowpack.mesh.boundaries:
            tag = bcs[0]
            bc_node = bcs[1]    
            if tag == 'top':
                self.RHS_vector[bc_node.numberinMesh-1] = 1  ##MODIF 1 BY BCs in SIF
                self.LHS_matrix[bc_node.numberinMesh-1,:] = 0
                self.LHS_matrix[bc_node.numberinMesh-1,bc_node.numberinMesh-1] = 1
                