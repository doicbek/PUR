# PUR

Lightweight pure pseudo-Cl estimator for HEALPix.

PUR computes pure E- and B-mode spherical harmonic coefficients from masked CMB
(Cosmic Microwave Background) polarization data in HEALPix format, and the
corresponding mode-coupling matrices needed to obtain unbiased power spectrum
estimates.

## Installation

### Requirements

- Python >= 3.8
- NumPy
- HEALPy
- Cython
- A C++ compiler (e.g. g++)

### Install

```bash
pip install .
```

Or for development (editable install):

```bash
pip install -e .
```

## Quick Start

### Compute pure E/B-mode alm from masked maps

```python
import numpy as np
import healpy as hp
import pur

# Load or generate T, Q, U maps and a binary mask
nside = 256
maps = np.array([T_map, Q_map, U_map])  # shape (3, Npix)
mask = np.ones(hp.nside2npix(nside))
mask[bad_pixels] = 0

# Compute pure spherical harmonic coefficients
alm = pur.map2alm_pure(maps, mask, lmax=500)
# alm[0] = T_alm, alm[1] = pure E_alm, alm[2] = pure B_alm

# Compute power spectra
cl = hp.alm2cl(alm)
```

### Compute the mode-coupling matrix

The mode-coupling matrix (MCM) describes how the mask mixes power between
multipoles. It is needed to deconvolve the mask effect from pseudo-Cl estimates.

```python
# Compute the mask power spectrum
mcl = hp.anafast(mask)

# Compute MCM with pure B-modes for both fields
mcm = pur.compute_mcm(mcl, lmax=500, pe1=False, pe2=False, pb1=True, pb2=True)

# Access specific coupling blocks:
# mcm['TT']['TT'] - temperature auto-coupling (lmax+1, lmax+1)
# mcm['EE']['EE'] - E-mode auto-coupling
# mcm['BB']['BB'] - B-mode auto-coupling
# mcm['EE']['BB'] - E-to-B leakage (should be suppressed by purification)
```

## API Reference

### `pur.map2alm_pure(maps, mask, lmax=None, mmax=None)`

Computes pure spherical harmonic coefficients of masked CMB polarization maps.
Uses the purification scheme from K. Smith (2006) to avoid E/B mixing caused
by incomplete sky coverage.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `maps` | array, shape (3, Npix) | Input T, Q, U maps in ring ordering |
| `mask` | array, shape (Npix,) | Sky mask (0 = masked, nonzero = observed) |
| `lmax` | int, optional | Maximum multipole l. Default: 3*nside - 1 |
| `mmax` | int, optional | Maximum m of the alm. Default: lmax |

**Returns:**

| Return | Type | Description |
|--------|------|-------------|
| `alms` | array, shape (3, Nalm) | Pure alm arrays: (almT, almE, almB) |

**Example:**

```python
import pur
import numpy as np
import healpy as hp

nside = 128
npix = hp.nside2npix(nside)
maps = np.array([np.random.randn(npix) for _ in range(3)])
mask = np.ones(npix)
mask[:1000] = 0  # mask out some pixels

alm = pur.map2alm_pure(maps, mask, lmax=200)
cl_TT, cl_EE, cl_BB, cl_TE, cl_EB, cl_TB = hp.alm2cl(alm)
```

---

### `pur.get_spinned_windows(w, lmax=None, mmax=None)`

Computes the spin-1 and spin-2 window functions from a mask, used internally
by `map2alm_pure`. These are the differential-operator-weighted windows
needed for the purification scheme.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `w` | array, shape (Npix,) | The mask in HEALPix ring ordering |
| `lmax` | int, optional | Maximum multipole l. Default: 3*nside - 1 |
| `mmax` | int, optional | Maximum m. Default: lmax |

**Returns:**

| Return | Type | Description |
|--------|------|-------------|
| `w1_plus` | array, shape (Npix,) | Spin-1 positive-parity window |
| `w1_minus` | array, shape (Npix,) | Spin-1 negative-parity window |
| `w2_plus` | array, shape (Npix,) | Spin-2 positive-parity window |
| `w2_minus` | array, shape (Npix,) | Spin-2 negative-parity window |

---

### `pur.compute_mcm(clmask, lmax, pe1=False, pe2=False, pb1=False, pb2=False)`

Computes the mode-coupling matrix (MCM) for cross-power spectrum estimation
between two fields observed through the same mask. The MCM accounts for
how the mask mixes power between different multipoles and between E and B modes.

When purification flags are enabled, the MCM is computed for the corresponding
pure pseudo-Cl estimator, which eliminates ambiguous modes at the boundary
of the observed region.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `clmask` | array, shape (lmax_mask+1,) | Power spectrum of the mask (from `hp.anafast(mask)`) |
| `lmax` | int | Maximum multipole for the coupling matrix. Must be <= len(clmask)-1 |
| `pe1` | bool | Enable pure E-mode for field 1. Default: False |
| `pe2` | bool | Enable pure E-mode for field 2. Default: False |
| `pb1` | bool | Enable pure B-mode for field 1. Default: False |
| `pb2` | bool | Enable pure B-mode for field 2. Default: False |

**Returns:**

| Return | Type | Description |
|--------|------|-------------|
| `mcm` | dict of dict | Nested dictionary of coupling matrices |

The returned dictionary has the following structure. Each value is an
`(lmax+1, lmax+1)` numpy array:

```
mcm['TT']['TT']  - T-T coupling
mcm['TE']['TE']  - T-E coupling
mcm['TB']['TB']  - T-B coupling
mcm['EE']['EE']  - E-E same-parity coupling
mcm['EE']['BB']  - B-to-E leakage (opposite-parity)
mcm['BB']['BB']  - B-B same-parity coupling
mcm['BB']['EE']  - E-to-B leakage (opposite-parity)
mcm['EB']['EB']  - E-B same-parity coupling
mcm['EB']['BE']  - B-E opposite-parity coupling
mcm['BE']['EB']  - E-B opposite-parity coupling
mcm['BE']['BE']  - B-E same-parity coupling
```

**Example:**

```python
import pur
import healpy as hp

mask = hp.read_map("mask.fits")
mcl = hp.anafast(mask)

# Auto-spectrum MCM with pure B-modes
mcm = pur.compute_mcm(mcl, lmax=500, pb1=True, pb2=True)

# Cross-spectrum MCM with different purification settings
mcm_cross = pur.compute_mcm(mcl, lmax=500, pe1=True, pe2=False, pb1=True, pb2=True)
```

---

### Convenience re-exports

PUR also re-exports the following from HEALPy for convenience:

- `pur.map2alm` - Standard (non-pure) spherical harmonic transform
- `pur.alm2cl` - Compute power spectra from alm coefficients

## Background

### The Problem

CMB polarization can be decomposed into E-modes (gradient-like, even parity) and
B-modes (curl-like, odd parity). On a full sky, this decomposition is clean.
However, real observations only cover part of the sky. When a mask is applied,
the standard pseudo-Cl estimator mixes E and B modes, contaminating the
faint B-mode signal with leaked E-mode power.

### The Solution

The **pure pseudo-Cl** method (Smith 2006) constructs modified estimators that
are orthogonal to the ambiguous modes at the mask boundary. This eliminates
E-to-B leakage at the cost of slightly higher variance. PUR implements this
approach by:

1. Computing spin-weighted window functions from the mask derivatives
2. Applying purification filters in harmonic space
3. Computing the corresponding mode-coupling matrices that account for the
   purification

### When to Use Pure Estimation

- **Pure B-modes** (`pb1=True, pb2=True`): Recommended for B-mode science
  (e.g., primordial gravitational waves, CMB lensing B-modes). Eliminates E-to-B
  leakage which would otherwise dominate the B-mode signal.
- **Pure E-modes** (`pe1=True, pe2=True`): Rarely needed since E-modes are
  much stronger than B-modes, so B-to-E leakage is usually negligible.
- **No purification**: Appropriate for temperature-only analysis or when E/B
  mixing is handled by other means.

## Credits

The framework of E- and B-purification and mode-coupling-matrix computation is outlined in:

- [K. Smith, Pseudo-Cl estimators which do not mix E and B modes (2006)](https://arxiv.org/abs/astro-ph/0511629)
- [J. Grain et al, Polarized CMB power spectrum estimation using the pure pseudo-cross-spectrum approach (2009)](https://arxiv.org/abs/0903.2350)

Parts of this package evolved from code provided in:

- [NaMaster](https://github.com/LSSTDESC/NaMaster) (LSST DESC)
- [pspy](https://github.com/simonsobs/pspy) (Simons Observatory)
- [healpy](https://github.com/healpy/healpy)

## Authors

- [Ari Cukierman](https://kipac.stanford.edu/people/ari-cukierman)
- [Dominic Beck](https://kipac.stanford.edu/people/dominic-beck)

## License

GPLv2
