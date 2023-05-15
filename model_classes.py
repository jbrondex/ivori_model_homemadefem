#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 15:51:05 2022

@author: brondexj
"""
### Here we define all classes required in the sif
### Later we could use def __init__(self) method to prescribe default value of some attributes

class Simulation:
    "Define the parameters of the simulation"
    
class Constants:
    "Set the constants of nature that are required in the model"

class Geometry:
    "Define the computational domain"
    
class Material:
    "Define material properties"       
    Material_list = []
   
    def __init__(self):
        self.__class__.Material_list.append(self)

class BodyForce:
    "Define body forces"
    BodyForce_list = []
   
    def __init__(self):
        self.__class__.BodyForce_list.append(self)

class InitialCondition:
    "Define initial conditions of solution fields"
    InitialCondition_list = []
   
    def __init__(self):
        self.__class__.InitialCondition_list.append(self)
        
class BoundaryCondition:
    "Define boundary conditions of solution fields"
    BC_list = []
   
    def __init__(self):
        self.__class__.BC_list.append(self)

class Solver:
    "Define Equation to solve and numerical parameters to enable solving"
    Solver_list = []
   
    def __init__(self):
        self.__class__.Solver_list.append(self)


        
    

