import numpy as np
import model_classes
import FEM_classes
from constants import (
    rhoi,
    rhoa,
    LCal,
    Lm,
    mH2O,
    kB,
    D0,
    ka,
    ki,
    Ca,
    Ci, 
    Tfus,
    g
    )

####~~~~~~~~~~~~####
#### Simulation ####
####~~~~~~~~~~~~####
Simu = model_classes.Simulation() #Create object Simu of class Simulation
Simu.Name = 'Simu_Papier_Fig5_SimsonCase6_SettlementOnly_101nodes_10d_Output1d_Implicit_dt15min_' #Prescribed by user
Simu.NumericalMethod = 'FEM' #FEM, FVM or FDM
Simu.Type = 'Transient' #Steady or Transient
    
# Only for transient simulations:
Simu.Timestepping_Method = 'theta' #Only theta for now
Simu.Theta = 1.0#0<=theta<=1 -> in particular, 0: full explicit, 0.5: Cranck-Nicholson, 1: full implicit
Simu.Timestep_Size = 15*60 # Delta_t in seconds (?)
Simu.Total_TimeSteps = 191#4*24*20#4*24*20#60*38 #60*24 # Total number of timesteps
Simu.Output_Interval = 191#4*24#12*38 # One output every X timesteps

Simu.Save_csv = True ### Do we want to save csv files ? True or False
Simu.Output_Path = "./Output/Output_Fig5"  #path towards where the output will be stored

####~~~~~~~~~~~~~~~~~~~~~####
#### Constants of nature ####
####~~~~~~~~~~~~~~~~~~~~~####
Cst = model_classes.Constants() #Create object Simu of class Simulation
Cst.rhoi = rhoi  # [kgm-3] ice density
Cst.rhoa = rhoa # [kgm-3] dry air density
Cst.LCal = LCal # [Jm-3] latent heat of sublimation Calonne et al. (2014)
Cst.Lm = Lm   # [Jkg-1] LCal/rhoi
Cst.mH2O = mH2O  # [kg] mass h2o molecule
Cst.kB = kB  # [JK-1] Boltzmann constant
Cst.D0 = D0# [m2s-1] Diffusion coefficient in air
Cst.ka = ka  # [Wm-1K-1] thermal conductivity air
Cst.ki = ki  # [Wm-1K-1] thermal conductivity ice
Cst.Ca = Ca # [Jkg-1K-1] heat capacity of dry air
Cst.Ci = Ci  # [Jkg-1K-1] heat capaity of ice
Cst.Tfus = Tfus  # [K] Melting temperature of water
Cst.g = g  # [ms-2] gravitational constant
    
####~~~~~~~~~~####
#### Snowpack ####
####~~~~~~~~~~####
Snowpack = model_classes.Geometry() #Create object Snowpack of class Geometry
Snowpack.Auto_Mesh = True # True if mesh must be generated automatically, false if it must be downloaded from external file
# if Snowpack_AutoMesh = True:
Snowpack.Htot = 0.5 #0.1 #Total Snowpack height [m]
Snowpack.Nnodes = 101 #201 #Total number of nodes
# else if Snowpack_AutoMesh = False:
#Snowpack.MeshFile = 'mesh_1D.mesh' #Name of mesh file to download

####~~~~~~~~~~####
#### Material ####
####~~~~~~~~~~####
Mat = model_classes.Material() #Create object Mat of class Material
#For the settlement solver
Mat.viscosity = 'MatLaw_ViscosityVionnetTconst_Phii' #Parameterization of viscosity
Mat.IsCompaction = True# Boolean activating/deactivating compaction when updating Phii with settlement solver
Mat.IsDeposition = False # Boolean activating/deactivating deposition when updating Phii with settlement solver

####~~~~~~~~~~~~####       
#### Body Force ####
####~~~~~~~~~~~~####
BF = model_classes.BodyForce() #Create object BF of class BodyForce
    
####~~~~~~~~~~~~~~~~~~~####
#### Initial Condition ####
####~~~~~~~~~~~~~~~~~~~####
IC = model_classes.InitialCondition() #Create object IC of class InitialCondition
#For the settlement solver
IC.Phii =  'USF_PhiiinitSimson_Coordinate' #'USF_Phiiinit_Coordinate' #Constant or defined through USF

####~~~~~~~~~~~~~~~~~~~~~####
#### Boundary Conditions ####
####~~~~~~~~~~~~~~~~~~~~~####
## BC TOP ##
BC1 = model_classes.BoundaryCondition() #Create object BC_Top of class BoundaryCondition
BC1.Tag = "top" #To which boundary does BC1 apply ?

# BC BOTTOM ##
BC2 = model_classes.BoundaryCondition() #Create object BC_Bot of class BoundaryCondition
BC2.Tag = "bot" #To which boundary does BC2 apply ?

###~~~~~~~~~####
### SOLVER  ####
###~~~~~~~~~####
####Solver 1
Solver1 = model_classes.Solver() #Create object Solver2 of class Solver
Solver1.ExecSolver = 1 # Execute this solver every N time steps
Solver1.DOFs = 1
Solver1.File = "Solvers.Stress_Solver" #In which file the solver is stored
Solver1.Name = "Stress_Solver" #Solver Name

#####Solver 2
Solver2 = model_classes.Solver() #Create object Solver2 of class Solver
Solver2.ExecSolver = 1 # Execute this solver every N time steps
Solver2.DOFs = 1
Solver2.File = "Solvers.Settlement_Solver" #In which file the solver is stored
Solver2.Name = "Settlement_Solver" #Solver Name

