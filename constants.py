# Standard constants of physics (also in sif for solvers)
rhoi = 917  # [kgm-3] ice density
rhoa = 1.335  # [kgm-3] dry air density
LCal = 2.6e9  # [Jm-3] latent heat of sublimation Calonne et al. (2014)
Lm = LCal / rhoi  # [Jkg-1] LCal/rhoi
mH2O = 2.991507e-26  # [kg] mass h2o molecule
kB = 1.38e-23  # [JK-1] Boltzmann constant
D0 = 2.036e-5  # [m2s-1] Diffusion coefficient in air
ka = 0.024  # [Wm-1K-1] thermal conductivity air
ki = 2.3  # [Wm-1K-1] thermal conductivity ice
Ca = 1005  # [Jkg-1K-1] heat capacity of dry air
Ci = 2000  # [Jkg-1K-1] heat capacity of ice
Tfus = 273  # [K] Melting temperature of water
g = 9.80665  # [ms-2] gravitational constant

# model geometry
Z_field = 0.5  # [m] height
Z_lab = 0.02  # [m] height

### Libbrecht (1999) Saturation Water Vapor Density
T_ref_L = 6150  # [K] reference temperature
a0 = 3.6636e12
a1 = -1.3086e8
a2 = -3.3793e6
f = 461.31

### Calonne et al. (2014) Saturation Water Vapor Density
rho_ref = 4.837e-3  # [kgm-3] reference density
T_ref_C = 273  # [K] reference temperature

#### Calonne et al. (2011) k_eff formula
ka0 = 0.024
ka1 = -1.23e-4
ka2 = 2.5e-6

### Hansen and Foslien (2015) Saturation Water Vapor Density
c1 = -2445.5646
c2 = 8.2312
c3 = -1.667006e-2
c4 = 1.20514e-5
c5 = -6.757169
c6 = 133.3224
R_v = 461.5  # (N-m)/(kg-K)

### snow viscosity formulation from Vionnet et al. (2012)
eta_0 = 7.62237e6  # [kgs-1-m-1]
c_eta = 250  # [kgm-3]
a_eta = 0.1  # [K-1]
b_eta = 0.023  # [m3kg-1]


### snow viscosity for non-linear Glen's flow law
D_rate_literature = 10e-6  # [s-1] Deformation rate for snow intermediate value as mentioned in Bartelt et al. (2002)

### power law to restrict ice volume growth
pl1 = 690
pl2 = 650
