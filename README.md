# PUR
lightweight pure pseudo-Cl estimator for HEALPix 

***Use with caution - not tested!***

## Installation
```
pip install . [--user]
```

## Usage
Compute pure E- and B-mode spherical harmonic coefficients
```
import pur
alm=pur.map2alm_pure(map,mask)
cl=hp.alm2cl(alm)
```

Compute mode-coupling matrix corresponding to these pseudeo-power spectra
```
import healpy
mcl=healpy.anafast(mask)
mcm=pur.compute_mcm(mcl,lmax=500,pe1=False,pe2=False,pb1=True,pb2=True) # This assumes non-pure E-modes and pure B-modes in the mode-coupling matrix
```

That's it!

## Credits
The framework of E- and B-purification and mode-coupling-matrix computation is outlined in
- [K. Smith, Pseudo-Cℓ estimators which do not mix E and B modes (2006)](https://arxiv.org/abs/astro-ph/0511629 )
- [J. Grain et al, Polarized CMB power spectrum estimation using the pure pseudo-cross-spectrum approach (2009)](https://arxiv.org/abs/0903.2350)

Parts of this package evolved from code provided in
- NaMaster (https://github.com/LSSTDESC/NaMaster)
- pspy (https://github.com/simonsobs/pspy)
- healpy (https://github.com/healpy/healpy)

## Authors
- [Ari Cukierman](https://kipac.stanford.edu/people/ari-cukierman)
- [Dominic Beck](https://kipac.stanford.edu/people/dominic-beck)
