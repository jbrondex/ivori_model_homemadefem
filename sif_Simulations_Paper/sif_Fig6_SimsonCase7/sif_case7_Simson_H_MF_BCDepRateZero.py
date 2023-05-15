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
Simu.Name = 'Simu_Papier_Fig6_SimsonCase7_H_MF_SettlementOn_DepOn_101nodes_20d_Output1d_Implicit_dt15min_SmoothTransition_BCDepRate0_' #Prescribed by user
Simu.NumericalMethod = 'FEM' #FEM, FVM or FDM
Simu.Type = 'Transient' #Steady or Transient
    
# Only for transient simulations:
Simu.Timestepping_Method = 'theta' #Only theta for now
Simu.Theta = 1.0#0<=theta<=1 -> in particular, 0: full explicit, 0.5: Cranck-Nicholson, 1: full implicit
Simu.Timestep_Size = 15*60 # Delta_t in seconds (?)
Simu.Total_TimeSteps = 4*24*20#4*24*20#4*24*20#60*38 #60*24 # Total number of timesteps
Simu.Output_Interval = 4*24#4*24#12*38 # One output every X timesteps

Simu.Save_csv = True ### Do we want to save csv files ? True or False
Simu.Output_Path = "./Output/Output_Fig6"  #path towards where the output will be stored

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
#For the heat and vapor solvers
Mat.Keff = 'MatLaw_keffCalonne_Phii' #Parameterization used for calculation of Keff (MatLaw Calonne/Hansen or USF). Can also be a constant (prescribe a float then).
Mat.Deff = 'MatLaw_DeffCalonne_Phii' #Parameterization used for calculation of Deff (MatLaw Calonne/Hansen or USF). Can also be a constant (prescribe a float then).
Mat.rhoCeff = 'MatLaw_rhoCeffWithAir_Phii' #Parameterization used for calculation of rhoCeff (MatLaw NoAir/WithAir or USF). Can also be a constant (prescribe a float then).
Mat.SWD = 'MatLaw_RhovSatLibbrecht_Temperature' #Parameterization used for calculation of rhovSat (MatLaw Clausius-Clapeyron/Libbrecht/Hansen or USF). Can also be a constant (prescribe a float then).
Mat.dSWDdT = 'MatLaw_dRhovSatdTLibbrecht_Temperature' #Parameterization used for calculation of drhovSatdT.Must be consistent with SWD
#For the settlement solver
Mat.viscosity = 'MatLaw_ViscosityVionnet_Phii_Temperature' #Parameterization of viscosity
Mat.IsCompaction = True# Boolean activating/deactivating compaction when updating Phii with settlement solver
Mat.IsDeposition = True # Boolean activating/deactivating deposition when updating Phii with settlement solver

####~~~~~~~~~~~~####       
#### Body Force ####
####~~~~~~~~~~~~####
BF = model_classes.BodyForce() #Create object BF of class BodyForce
    
####~~~~~~~~~~~~~~~~~~~####
#### Initial Condition ####
####~~~~~~~~~~~~~~~~~~~####
IC = model_classes.InitialCondition() #Create object IC of class InitialCondition
#For the heat solver
IC.Temperature = 263 #Constant or defined through USF
#IC.Temperature = 'USF_TinitSchurholtScenario2_Coordinate'#273 #Constant or defined through USF

#For the vapor solver
IC.Rhov = 'MatLaw_RhovSatLibbrecht_Temperature' #Constant or defined through USF

#For the settlement solver
IC.Phii =  'USF_PhiiinitSimson_Coordinate' #'USF_Phiiinit_Coordinate' #Constant or defined through USF
####~~~~~~~~~~~~~~~~~~~~~####
#### Boundary Conditions ####
####~~~~~~~~~~~~~~~~~~~~~####
## BC TOP ##
BC1 = model_classes.BoundaryCondition() #Create object BC_Top of class BoundaryCondition
BC1.Tag = "top" #To which boundary does BC1 apply ?
#For the heat solver
BC1.BCType_Heat = 'Given Temperature' #Can be Heat Flux (Neumann), 'Given Temperature' (Dirichlet), 'Adiabatic' (Natural BC) or 'Energy Balance' (Robin, not implemented yet)
BC1.Temperature = 253 #'USF_BCtopSchurholtscenario1_Time' #If BCType_Heat = 'Given Temperature' : Constant or defined through USF_Ttop_Variables
#BC1.HeatFlux = -0.9 #'USF_HeatFluxtop_Variables' #If BCType_Heat = 'Heat Flux' : Constant or defined through USF_HeatFluxtop_Variables'

###For the diagnostic of deposition rate of the hansen model
BC1.BCType_DepositionRate = 'Given Deposition Rate'
BC1.DepositionRate = 0

# BC BOTTOM ##
BC2 = model_classes.BoundaryCondition() #Create object BC_Bot of class BoundaryCondition
BC2.Tag = "bot" #To which boundary does BC2 apply ?
#For the heat solver
BC2.BCType_Heat = 'Given Temperature' #Can be Heat Flux (Neumann), 'Given Temperature' (Dirichlet), 'Adiabatic' (Natural BC) or 'Energy Balance' (Robin, not implemented yet)
BC2.Temperature = 273 #If BCType_Heat = 'Given Temperature' : Constant or defined through USF_Tbot_Variables
#BC2.HeatFlux = 19.9 #'USF_HeatFluxbot_Variables' #If BCType_Heat = 'Heat Flux' : Constant or defined through USF_HeatFluxbot_Variables'

###For the diagnostic of deposition rate of the hansen model
BC2.BCType_DepositionRate = 'Given Deposition Rate'
BC2.DepositionRate = 0
###~~~~~~~~~####
### SOLVER  ####
###~~~~~~~~~####
###Solver 1
Solver1 = model_classes.Solver() #Create object Solver1 of class Solver
Solver1.ExecSolver = 1 # Execute this solver every N time steps
Solver1.File = "Solvers.HansenHeatVap_Solver_2DOFs_MixedForm" #In which file the solver is stored
Solver1.Name = "HansenHeatVap" #Solver Name
Solver1.DOFs = 2 ### For this solver DOF1= H, DOF2= T
Solver1.MassLumping = False ## Do we want to lump the mass matrix ?
Solver1.Nonlin_Max_It = 100 #Maximum number of non-linear iterations
Solver1.Nonlin_Thres = 1e-5 #convergence tolerance of non-linearization
Solver1.Nonlin_Newton_After_It = 101  #Linearization method switch from Picard to Newton after Iteration 101 (so never for now)
Solver1.Nonlin_Newton_After_Thres = 1e-6 #Linearization method switch from Picard to Newton after when gap between two iterations is smaller than 1e-6 (so never for now)
Solver1.Relaxation_Factor = 1.0 #Solution used for linearization at iteration i+1 is w*U_i+1 + (1-w)Ui. If 1 means no relaxation.
#
######Solver 2
Solver2 = model_classes.Solver() #Create object Solver2 of class Solver
Solver2.ExecSolver = 1 # Execute this solver every N time steps
Solver2.DOFs = 1
Solver2.File = "Solvers.DepositionRateHansen_Solver" #In which file the solver is stored
Solver2.Name = "DepositionRateHansen" #Solver Name
Solver2.MassLumping = True ## Do we want to lump the mass matrix ?

####Solver 3
Solver3 = model_classes.Solver() #Create object Solver2 of class Solver
Solver3.ExecSolver = 1 # Execute this solver every N time steps
Solver3.DOFs = 1
Solver3.File = "Solvers.Stress_Solver" #In which file the solver is stored
Solver3.Name = "Stress_Solver" #Solver Name

######Solver 4
Solver4 = model_classes.Solver() #Create object Solver2 of class Solver
Solver4.ExecSolver = 1 # Execute this solver every N time steps
Solver4.DOFs = 1
Solver4.File = "Solvers.Settlement_Solver" #In which file the solver is stored
Solver4.Name = "Settlement_Solver" #Solver Name

#####Solver 5
Solver5 = model_classes.Solver() #Create object Solver2 of class Solver
Solver5.ExecSolver = 1 # Execute this solver every N time steps
Solver5.DOFs = 1
Solver5.File = "Solvers.EnergyConservation_Solver" #In which file the solver is stored
Solver5.Name = "EnergyConservation" #Solver Name