'''
Author: Ari Cukierman (ajcukier@stanford.edu)
'''
import healpy as hp
import numpy as np

def get_spinned_windows(w,lmax=None,mmax=None):    
    nside = hp.npix2nside(len(w))
    if lmax is None:
        lmax = 3 * nside - 1
    if mmax is None:
        mmax = lmax

    wlm = hp.map2alm(w, lmax=lmax, mmax=mmax)
    ell = np.arange(lmax+1)
    filter_1 = -np.sqrt((ell+1.)*ell)
    filter_2 = -np.sqrt((ell-1.)*(ell+2.))
    
    filter_1[:1] = 0
    filter_2[:1] = 0
    
    wlm1_e = hp.almxfl(wlm, filter_1, mmax=mmax)
    wlm2_e = hp.almxfl(wlm, filter_2, mmax=mmax)
    wlm1_b = np.zeros_like(wlm1_e)
    wlm2_b = np.zeros_like(wlm2_e)

    w1_full = hp.alm2map_spin(np.array([wlm1_e,wlm1_b]), nside, 1, lmax=lmax, mmax=mmax)
    w2_full = hp.alm2map_spin(np.array([wlm2_e,wlm2_b]), nside, 2, lmax=lmax, mmax=mmax)

    w1 = []
    w2 = []

    binw=np.zeros_like(w)
    binw[w!=0]=1

    w1_plus = w1_full[0]*binw
    w1_minus = w1_full[1]*binw
    w2_plus = w2_full[0]*binw
    w2_minus = w2_full[1]*binw

    return w1_plus,w1_minus,w2_plus,w2_minus

def map2alm_pure(maps,mask,lmax=None,mmax=None):
    """Computes pure alm of a masked Healpix map. The input maps must all be
    in ring ordering. See K. Smith (2006, astro-ph/0608662) for details
    on the purfication scheme.
    Parameters
    ----------
    maps : array-like, shape (Npix,) or (n, Npix)
      The input map or a list of n input maps. Must be in ring ordering.
    mask : array-like, shape (Npix,) 
      The mask which will be multiplied to the input maps.
    lmax : int, scalar, optional
      Maximum l of the power spectrum. Default: 3*nside-1
    mmax : int, scalar, optional
      Maximum m of the alm. Default: lmax
    iter : int, scalar, optional
      Number of iteration (default: 3)
    Returns
    -------
    alms : array or tuple of array
      tuple of 3 alm (almT, almE, almB).
    """
    nside = hp.npix2nside(len(mask))

    if lmax is None:
        lmax = 3 * nside - 1
    if mmax is None:
        mmax = lmax

    w1_plus,w1_minus,w2_plus,w2_minus = get_spinned_windows(mask,lmax,mmax)

    p2 = np.array([mask*maps[1],mask*maps[2]])
    p1 = np.array([w1_plus*maps[1] + w1_minus*maps[2], w1_plus*maps[2] - w1_minus*maps[1]])
    p0 = np.array([w2_plus*maps[1] + w2_minus*maps[2], w2_plus*maps[2] - w2_minus*maps[1]])

    alm = hp.map2alm(maps[0]*mask,lmax=lmax,mmax=mmax)
    
    s2eblm = hp.map2alm_spin(p2,2,lmax=lmax,mmax=mmax)
    s1eblm = hp.map2alm_spin(p1,1,lmax=lmax,mmax=mmax)

    s0eblm = s1eblm.copy()
    s0eblm[0] = hp.map2alm(p0[0],lmax=lmax,mmax=mmax)
    s0eblm[1] = hp.map2alm(p0[1],lmax=lmax,mmax=mmax)

    ell = np.arange(lmax+1)
    filter_1 = np.zeros(lmax+1)
    filter_2 = np.zeros(lmax+1)

    filter_1[2:] = 2*np.sqrt(1./((ell[2:]+2.)*(ell[2:]-1.)))
    filter_2[2:] = np.sqrt(1./((ell[2:]+2.)*(ell[2:]+1.)*ell[2:]*(ell[2:]-1.)))
    
    for k in range(2):
        s1eblm[k] = hp.almxfl(s1eblm[k],filter_1,mmax=mmax)
        s0eblm[k] = hp.almxfl(s0eblm[k],filter_2,mmax=mmax)

    elm_p = s2eblm[0] + s1eblm[0] + s0eblm[0]
    blm_b = s2eblm[1] + s1eblm[1] + s0eblm[1]

    return np.array([alm,elm_p,blm_b])
