import numpy as np

def evaluate_basis_ref(xi):
    return np.array([(1-xi)/2., (xi+1)/2.])
   
def evaluate_dbasis_ref(xi):
    return np.array([-1/2.,1/2.])
   
def evaluate_ddbasis_ref(xi):
    return np.array([0,0])
   
