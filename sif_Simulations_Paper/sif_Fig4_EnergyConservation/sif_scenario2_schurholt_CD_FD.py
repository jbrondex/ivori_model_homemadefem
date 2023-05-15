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
#Simu.Name = 'Simu_Schurholt_Scenario2_CalParCalEq_dt5mn_38h_Implicit_Kreactlumped_cdiag_' #Prescribed by user
Simu.Name = 'Simu_Papier_Fig4_CD_FD_ProperMass_cLumped_MassnotLumped_CalPar_Alpha5minus3_NoFluxBCs_SettlementOff_DepOn_201nodes_5d_Output12h_Implicit_dt1min_' #Prescribed by user
Simu.NumericalMethod = 'FEM' #FEM, FVM or FDM
Simu.Type = 'Transient' #Steady or Transient
    
# Only for transient simulations:
Simu.Timestepping_Method = 'theta' #Only theta for now
Simu.Theta = 1.0#0<=theta<=1 -> in particular, 0: full explicit, 0.5: Cranck-Nicholson, 1: full implicit
Simu.Timestep_Size = 60 # Delta_t in seconds (?)
Simu.Total_TimeSteps = 60*24#15*4*24*5+1#60*38 #60*24 # Total number of timesteps
Simu.Output_Interval = 60*2#15*4*12#12*38 # One output every X timesteps

Simu.Save_csv = True ### Do we want to save csv files ? True or False
Simu.Output_Path = "./Output/Output_Fig4"  #path towards where the output will be stored

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
Snowpack.Htot = 1 #0.1 #Total Snowpack height [m]
Snowpack.Nnodes = 201 #201 #Total number of nodes
# else if Snowpack_AutoMesh = False:
#Snowpack.MeshFile = 'mesh_1D.mesh' #Name of mesh file to download

####~~~~~~~~~~####
#### Material ####
####~~~~~~~~~~####
Mat = model_classes.Material() #Create object Mat of class Material
#For the heat and vapor solvers
Mat.Keff = 'MatLaw_keffCalonne_Phii' #Parameterization used for calculation of Keff (MatLaw Calonne/Hansen or USF). Can also be a constant (prescribe a float then).
Mat.Deff = 'MatLaw_DeffCalonne_Phii' #Parameterization used for calculation of Deff (MatLaw Calonne/Hansen or USF). Can also be a constant (prescribe a float then).
Mat.rhoCeff = 'MatLaw_rhoCeffNoAir_Phii' #Parameterization used for calculation of rhoCeff (MatLaw NoAir/WithAir or USF). Can also be a constant (prescribe a float then).
Mat.SWD = 'MatLaw_RhovSatLibbrecht_Temperature' #Parameterization used for calculation of rhovSat (MatLaw Clausius-Clapeyron/Libbrecht/Hansen or USF). Can also be a constant (prescribe a float then). 
Mat.dSWDdT = 'MatLaw_dRhovSatdTLibbrecht_Temperature' #Parameterization used for calculation of drhovSatdT.Must be consistent with SWD   
#For the settlement solver
Mat.viscosity = 'MatLaw_ViscosityVionnet_Phii_Temperature' #Parameterization of viscosity
Mat.IsCompaction = False# Boolean activating/deactivating compaction when updating Phii with settlement solver
Mat.IsDeposition = True # Boolean activating/deactivating deposition when updating Phii with settlement solver

####~~~~~~~~~~~~####       
#### Body Force ####
####~~~~~~~~~~~~####
BF = model_classes.BodyForce() #Create object BF of class BodyForce
#For the heat and vapor solvers
BF.s = 3770 #volumic SSA [m^2/m^3]
BF.alpha = 5e-3##10**-5 #sticking coefficient
BF.vkin = 'MatLaw_vkin_Temperature' #kinetic velocity. Constant or defined through "MatLaw_vkin_Temperature" or USF to write.

####~~~~~~~~~~~~~~~~~~~####
#### Initial Condition ####
####~~~~~~~~~~~~~~~~~~~####
IC = model_classes.InitialCondition() #Create object IC of class InitialCondition
#For the heat solver
#IC.Temperature = 263 #Constant or defined through USF
IC.Temperature = 'USF_TinitSchurholtScenario2_Coordinate'#273 #Constant or defined through USF

#For the vapor solver
IC.Rhov = 'MatLaw_RhovSatLibbrecht_Temperature' #Constant or defined through USF
        
#For the settlement solver
IC.Phii =  'USF_PhiiinitSchurholtScenario2_Coordinate' #'USF_Phiiinit_Coordinate' #Constant or defined through USF
    
####~~~~~~~~~~~~~~~~~~~~~####
#### Boundary Conditions ####
####~~~~~~~~~~~~~~~~~~~~~####
## BC TOP ##
BC1 = model_classes.BoundaryCondition() #Create object BC_Top of class BoundaryCondition
BC1.Tag = "top" #To which boundary does BC1 apply ?
#For the heat solver
BC1.BCType_Heat = 'Adiabatic'
# BC1.BCType_Heat = 'Given Temperature' #Can be Heat Flux (Neumann), 'Given Temperature' (Dirichlet), 'Adiabatic' (Natural BC) or 'Energy Balance' (Robin, not implemented yet)
# BC1.Temperature = 253 #'USF_BCtopSchurholtscenario1_Time' #If BCType_Heat = 'Given Temperature' : Constant or defined through USF_Ttop_Variables
# BC1.HeatFlux = -0.9 #'USF_HeatFluxtop_Variables' #If BCType_Heat = 'Heat Flux' : Constant or defined through USF_HeatFluxtop_Variables'
#For the vapor solver
BC1.BCType_Rhov = 'No Flux'
# BC1.BCType_Rhov = 'Given Rhov' #Can be 'Rhov Flux' (Neumann), 'Given Rhov' (Dirichlet), 'No Flux' (Natural BC) or 'Rhov Balance' (Robin, not implemented yet)
# BC1.Rhov = 'MatLaw_RhovSatLibbrecht_Temperature' #If BCType_Rhov = 'Given Rhov' : Constant or defined through USF_Rhovbot_Variables
#BC1.RhovFlux = -1.89e-8 #If BCType_Rhov = 'Rhov Flux' : Constant or defined through USF_RhovFluxbot_Variables'

# BC BOTTOM ##
BC2 = model_classes.BoundaryCondition() #Create object BC_Bot of class BoundaryCondition
BC2.Tag = "bot" #To which boundary does BC2 apply ?
#For the heat solver
BC2.BCType_Heat = 'Adiabatic'
# BC2.BCType_Heat = 'Given Temperature' #Can be Heat Flux (Neumann), 'Given Temperature' (Dirichlet), 'Adiabatic' (Natural BC) or 'Energy Balance' (Robin, not implemented yet)
# BC2.Temperature = 273 #If BCType_Heat = 'Given Temperature' : Constant or defined through USF_Tbot_Variables
#BC2.HeatFlux = 19.9 #'USF_HeatFluxbot_Variables' #If BCType_Heat = 'Heat Flux' : Constant or defined through USF_HeatFluxbot_Variables'
#For the vapor solver
BC2.BCType_Rhov = 'No Flux'
# BC2.BCType_Rhov = 'Given Rhov' #Can be 'Rhov Flux' (Neumann), 'Given Rhov' (Dirichlet), 'No Flux' (Natural BC) or 'Rhov Balance' (Robin, not implemented yet)
# BC2.Rhov = 'MatLaw_RhovSatLibbrecht_Temperature' #If BCType_Rhov = 'Given Rhov' : Constant or defined through USF_Rhovbot_Variables
#BC2.RhovFlux = 0.0# 'USF_RhovFluxbot_Variables' #If BCType_Rhov = 'Rhov Flux' : Constant or defined through USF_RhovFluxbot_Variables'

###~~~~~~~~~####
### SOLVER  ####
###~~~~~~~~~####
#####Solver 1
Solver1 = model_classes.Solver() #Create object Solver1 of class Solver
Solver1.ExecSolver = 1 # Execute this solver every N time steps
Solver1.File = "Solvers.CalonneHeat_Solver_1DOFs_KreactLumped" #In which file the solver is stored
Solver1.Name = "CalonneHeat" #Solver Name
Solver1.DOFs = 1 ### For this solver DOF1 = T
Solver1.MassLumping = True ## Do we want to lump the mass matrix ?
Solver1.Nonlin_Max_It = 100 #Maximum number of non-linear iterations
Solver1.Nonlin_Thres = 1e-5 #convergence tolerance of non-linearization
Solver1.Nonlin_Newton_After_It = 101  #Linearization method switch from Picard to Newton after Iteration 101 (so never for now)
Solver1.Nonlin_Newton_After_Thres = 1e-6 #Linearization method switch from Picard to Newton after when gap between two iterations is smaller than 1e-6 (so never for now)
Solver1.Relaxation_Factor = 1.0 #Solution used for linearization at iteration i+1 is w*U_i+1 + (1-w)Ui. If 1 means no relaxation. 

#######Solver 2
Solver2 = model_classes.Solver() #Create object Solver1 of class Solver
Solver2.ExecSolver = 1 # Execute this solver every N time steps
Solver2.File = "Solvers.CalonneVap_Solver_1DOFs_KreactLumped" #In which file the solver is stored
Solver2.Name = "CalonneVap" #Solver Name
Solver2.DOFs = 1 ### For this solver DOF1= Rhov
Solver2.MassLumping = True ## Do we want to lump the mass matrix ?
Solver2.Nonlin_Max_It = 100 #Maximum number of non-linear iterations
Solver2.Nonlin_Thres = 1e-5 #convergence tolerance of non-linearization
Solver2.Nonlin_Newton_After_It = 101  #Linearization method switch from Picard to Newton after Iteration 101 (so never for now)
Solver2.Nonlin_Newton_After_Thres = 1e-6 #Linearization method switch from Picard to Newton after when gap between two iterations is smaller than 1e-6 (so never for now)
Solver2.Relaxation_Factor = 1.0 #Solution used for linearization at iteration i+1 is w*U_i+1 + (1-w)Ui. If 1 means no relaxation.

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
