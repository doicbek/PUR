'''
Author: Dominic Beck (dobeck@stanford.edu)
'''
import healpy as hp
import numpy as np
from libcpp.vector cimport vector

cdef extern from "libmcm.cpp":
    cdef void compute_mcm_namaster(int lmax, int lmax_mask, int npcl, vector[double] pcl_masks, int s1, int s2, int pure_e1, int pure_b1, int pure_e2, int pure_b2, int do_teb, int l_toeplitz, int l_exact, int dl_band, vector[vector[vector[double]]] &xi_00, vector[vector[vector[vector[double]]]] &xi_0s, vector[vector[vector[vector[double]]]] &xi_pp, vector[vector[vector[vector[double]]]] &xi_mm)

def compute_mcm(clmask,lmax,pe1=False,pe2=False,pb1=False,pb2=False):
    """Compute the mode-coupling matrix for masked pseudo-Cl estimation.

    Computes the coupling matrix that describes how a sky mask mixes power
    between multipoles and between E/B polarization modes. Supports pure
    pseudo-Cl estimation which eliminates E/B leakage.

    Parameters
    ----------
    clmask : array-like, shape (lmax_mask+1,)
        Power spectrum of the mask, as returned by ``healpy.anafast(mask)``.
        This array is not modified.
    lmax : int
        Maximum multipole of the output coupling matrix.
        Must be <= len(clmask) - 1.
    pe1 : bool, optional
        Use pure E-mode estimator for field 1. Default: False.
    pe2 : bool, optional
        Use pure E-mode estimator for field 2. Default: False.
    pb1 : bool, optional
        Use pure B-mode estimator for field 1. Default: False.
    pb2 : bool, optional
        Use pure B-mode estimator for field 2. Default: False.

    Returns
    -------
    mcm : dict of dict of ndarray
        Nested dictionary of (lmax+1, lmax+1) coupling matrices.
        Keys: 'TT', 'TE', 'TB', 'EE', 'BB', 'EB', 'BE'.

    Raises
    ------
    ValueError
        If lmax exceeds the length of the mask power spectrum.
    """
    cdef int lmax_,lmax_mask
    cdef int s1, s2

    lmax_=lmax

    lmax_mask=len(clmask)-1

    if lmax_ > lmax_mask:
        raise ValueError(f"lmax ({lmax}) must be <= len(clmask)-1 ({lmax_mask}). "
                         "Provide a mask power spectrum with enough multipoles.")

    s1=2
    s2=2

    cdef vector[vector[vector[double]]] xi_00
    cdef vector[vector[vector[vector[double]]]] xi_0s
    cdef vector[vector[vector[vector[double]]]] xi_pp
    cdef vector[vector[vector[vector[double]]]] xi_mm

    cdef vector[double] clmask_c = vector[double](lmax_mask+1,0.)
    ll=np.arange(lmax_mask+1)
    clmask_norm = clmask*(2*ll+1)/(4*np.pi)
    for l in range(lmax_mask+1):
        clmask_c[l]=clmask_norm[l]

    compute_mcm_namaster(lmax_,lmax_mask,1,clmask_c,s1,s2,pe1,pb1,pe2,pb2,1,-1,-1,-1,xi_00,xi_0s,xi_pp,xi_mm)
    
    cm={}
    cm['TT']={}
    cm['TT']['TT']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['TE']={}
    cm['TE']['TE']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['TB']={}
    cm['TB']['TB']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['EE']={}
    cm['EE']['EE']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['EE']['BB']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['BB']={}
    cm['BB']['EE']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['BB']['BB']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['EB']={}
    cm['EB']['EB']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['EB']['BE']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['BE']={}
    cm['BE']['EB']=np.zeros((lmax+1,lmax+1), dtype=np.float64)
    cm['BE']['BE']=np.zeros((lmax+1,lmax+1), dtype=np.float64)

    sgn=1
    if s1+s2==0: sgn=-1
    
    for l2 in range(lmax+1):
        for l3 in range(lmax+1):
            fac=(2*l3+1.)*sgn
            cm['TT']['TT'][l2][l3]=fac*xi_00[0][l2][l3]
            cm['TE']['TE'][l2][l3]=fac*xi_0s[0][pe2][l2][l3]
            cm['TB']['TB'][l2][l3]=fac*xi_0s[0][pb2][l2][l3]
            cm['EE']['BB'][l2][l3]=fac*xi_mm[0][pe1+pe2][l2][l3]
            cm['EB']['BE'][l2][l3]=-fac*xi_mm[0][pe1+pb2][l2][l3]
            cm['BE']['EB'][l2][l3]=-fac*xi_mm[0][pb1+pe2][l2][l3]
            cm['BB']['EE'][l2][l3]=fac*xi_mm[0][pb1+pb2][l2][l3]
            cm['EE']['EE'][l2][l3]=fac*xi_pp[0][pe1+pe2][l2][l3]
            cm['EB']['EB'][l2][l3]=fac*xi_pp[0][pe1+pb2][l2][l3]
            cm['BE']['BE'][l2][l3]=fac*xi_pp[0][pb1+pe2][l2][l3]
            cm['BB']['BB'][l2][l3]=fac*xi_pp[0][pb1+pb2][l2][l3]

    return cm