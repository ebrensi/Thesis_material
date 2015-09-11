from scipy import zeros, dot, random, mat, linalg, diag, sqrt, sum, hstack, ones
from scipy.linalg import norm, eig

def band_Arnoldi(nmax, R, MultA, tol={}, n0=0, result0={}):
    
    if not tol:  # if tol is empty...
        tol.defl_flag = 1;
        tol.defl_tol = sqrt(eps);
        tol.normA_flag = 1;
    