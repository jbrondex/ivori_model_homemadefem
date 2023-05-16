# A finite-element framework to explore the numerical solution of the coupled problem of heat conduction, water vapor diffusion and settlement in dry snow

**Authors:** Julien Brondex; Kevin Fourteau, Marie Dumont, Pascal Hagenmuller, Neige Calonne, Francois Tuzet, and Henning Lowe

**Contact author:** Julien Brondex - julien.brondex@meteo.fr 

*Univ. Grenoble Alpes, Université de Toulouse, Météo-France, CNRS, CNRM, Centre Etudes de la Neige (Snow Research Center)*

## Motivations

We aim at pinpointing the most suitable numerical treatments of the coupled problem of heat conduction, water vapor diffusion and settlement in dry snow,
based on criteria of time step requirements, conservation of energy, and conservation of mass. 
To this end, we unify previous model developments within a comprehensive, stand-alone, finite-element core written in python. 
Each step of the numerical formulation and solution is coded internally, except for the inversion of the linearized system of equations 
which relies on the standard numpy linear algebra library. This work is at the origin of an article submitted to GMD on 2023/05/18. 

## Quick install

Install the following packages, preferably in a conda or virtual environment:

python==3.8.10

numpy==1.17.4

matplotlib==3.1.2

pandas==0.25.3

vtk==9.1.0

six==1.14.0

## General Structure of the Model
![image](/images/Architecture_solvers.png)


## Structure github repository

```
├── LICENSE                             <- To be created
│
├── README.md                           <- File describing the github repository. You are currently reading this file.
│
├── Figures                            
│   ├── Fig1_Paper                      <- Contains First plot (i.e, Figure 3) of submitted version of the paper (Figs. 1 and 2 are images).
│   :    
│   ├── Fign_Paper                      <- Contains plot n (i.e, Figure n+2) of submitted version of the paper (Figs. 1 and 2 are images).               
│   :                                                
│   ├── PythonScripts                   <- Contains all the scripts to produce the Figs from the simulation outputs
│
├── images                              <- Contains images corresponding to Fig. 1 and 2 of the submitted version of the paper
│
├── Output                              <- Contains the output of simulations               
│   ├── Output_Fig1                  
│   :                                
│                                    
├── sif_Simulations_Paper               <- Contains all the sifs (simulation input files) corresponding to all the simulations described in the paper (+ some additionnal simulations). Description of a sif file is given below. 
│   ├── sif_Fig1_CC3DOFsVariousVersion  <- Contains all the sifs of simulations illustrated in Fig. 3 of the submitted version of the paper. 
│   :                           
│   :                           
│
├── Solvers                             <- Each of the file in this directory corresponds to a solver. Solvers are described carefully in the paper.
│   :                          
│
├── constants.py                        <- Where all the constants of the model are defined.
│
├── FEM_classes.py                      <- Where all the required classes are defined (don't touch).
│
├── Gauss_Integ.py                      <- Parameters regarding Gaussian points required for Gaussian quadratures (don't touch).
│
├── main.py                             <- The executable file. Basic users simply have to change the 'Directory_sif' and 'Name_sif' lines to specify which sif (i.e., simulation) they want to run.
│
├── MatLaw.py                           <- Where parameterizations of material properties (e.g., effective parameters, saturation water vapor density) are defined.                   
│
├── model_classes.py                    <- Definition of classes to manipulate sections of the sif (don't touch).
│
├── Ref_Elements.py                     <- Parameters regarding shape/test functions and their deivatives in the reference element (don't touch).
│
└── USF.py                              <- Where user-defined parameterizations (e.g., initial distribution of a field as a function of space, or evolution of a boundary condition as a function of time) are defined.
```
## Run the code

The two important files to run the code are the main.py and the sif_<...>.py.
The main.py is the executable. The sif_<...>.py is the configuration file in
which the user can prescribe all the parameters of the simulation 
(total duration, time step, material properties, initial and boundary conditions, solvers, ... ).
More details on how to write a sif are given in the section below.

In the main.py file, the only two lines that the user has to change are the "Directory_sif" and "Name_sif".
They correspond to respectively, the path to the directory where the sif_<...>.py file to run is located, and the full
name of this sif_<...>.py file. Once these lines have been filled, the main.py can be run. In some situations, the importlib.import_module might not work.
In that case, comment lines 32 to 34 of the main.py file and use instead line 26 to import the sif_<...>.py file you want to run.

## How to write a sif ?
The sif_<...>.py (solver input file) basically corresponds to the user interface. All the sif_<...>.py that have been used 
to run the simulations presented in the article are stored in the folders 'sif_Simulations_Paper'. 
The user can re-run these sif_<...>.py as such, or use them as examples to design its own numerical experiment and
write the corresponding sif_<...>.py. The sif_<...>.py files available in the folders 'sif_Simulations_Paper' are generally well-commented,
and they should be easy to adapt. Below, we detail the various sections of a sif_<...>.py file.

**Simulation**

This section enables to set the simulation parameters: name of simulation, value of theta (temporal discretization),
time step, total number of time step, output intervals, whether or not the outputs should be saved, and path to the folder where they should be saved.
For now, the only possibility for keywords, resp., Simu.NumericalMethod, Simu.Type, and Simu.Timestepping_Method, is, resp., 'FEM', 'Transient', and
'theta'.

**Constants of nature**

This section enables to set all the constants required in the model. These constants are defined in the constants.py file.

**Snowpack**

This section enables to prescribe the initial geometry. The model being 1D, it basically corresponds
to the initial total height of the snowpack and total number of nodes. This will automatically generate a regular mesh (Snowpack.Auto_Mesh = True). 
We are working on the possibility to prescribe an initial irregular mesh from a file <...>.mesh
(keyword Snowpack.MeshFile) but this option is not ready yet.

**Material**

This section enables to prescribe all material properties. The main material properties are the
effective heat conductivity, heat capacity, vapor diffusion coefficient, and viscosity, as well as the model used
to calculate the saturation water vapor density and its T-derivative from temperature. In this section, the feedback
of deposition and/or settlement on the ice phase can be activated or de-activated (keywords IsCompaction and IsDeposition). Material properties can either be prescribed as constants
or as functions of other field variables. In the latter case, they need to be defined in the MatLaw.py file. In the sif_<...>.py file, 
the keyword must then be a string with the following syntax: 'MatLaw_<_Name_> _ <_Variable1_>_<_Variable2_>'.

**Body Force**

This section enables to prescribe parameters related to body forces. This typically corresponds to source terms
of equations to be solved. These parameters can either be prescribed as constants or as functions of other
field variables. In the latter case, they need to be defined in the MatLaw.py file. In the sif_<...>.py file, 
the keyword must then be a string with the following syntax: 'MatLaw_<_Name_> _ <_Variable1_>_<_Variable2_>'.

**Initial Conditions**

This section enables the prescription of the initial conditions. They can either be prescribed as constants or
as functions of other variables. In the latter case, they can be defined in the USF.py file (e.g., for a space-dependent initial field),
or in the MatLaw.py file (e.g., when the vapor density is initially saturated). In the sif_<...>.py file, 
the keyword must then be a string with the following syntax: 'USF_<_Name_> _ <_Variable1_>_ <_Variable2_>' (resp. 'MatLaw_<_Name_> _ <_Variable1_>_<_Variable2_>')

**Boundary Conditions**

This section enables the prescription of the boundary conditions. BC1 corresponds to the top BC, and BC2 to the bottom BC. In the current stage of development, 
BCs are required on heat and vapor whenever the system of Calonne or Hansen is solved. In the absence of any prescription, the natural (i.e., No Flux)
BC applies. Otherwise, BCs can be of Dirichlet and/or Neumann types (Robin BCs are not implemented yet). Examples of BC implementation can be found in available
sif_<...>.py files. Again, BCs can either be prescribed as constants or as functions of other variables (or time). In the latter case, they can be defined in the USF.py file (e.g., a time-dependent BC),
or in the MatLaw.py file (e.g., vapor density is saturated at a boundary). In the sif_<...>.py file, the keyword must then be a string with the following syntax: 'USF_<_Name_> _ <_Variable1_>_ <_Variable2_>' (resp. 'MatLaw_<_Name_> _ <_Variable1_>_<_Variable2_>')

**Solvers**

It is the section in which all solvers required to solve the model equations are called. Solvers will be
executed in the same order as they appear in the sif_<...>.py. All solvers are stored in the "Solvers" folder.
The numerical methods on which they rely are thoroughly described in the submitted paper. Each solver 
comes with specific keywords, the meaning of which can be deduced from the comments in available sif_<...>.py file.  