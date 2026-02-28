# Bugs Found and Fixed

This document describes the bugs identified in the PUR codebase and the fixes applied.

## Critical Bugs

### 1. Mode-coupling matrix ignores field-1 purification settings (`mcm.pyx`)

**Symptom:** Setting `pe1=True` or `pb1=True` in `compute_mcm()` had no effect on the output. Only `pe2` and `pb2` influenced the coupling matrices.

**Root cause:** The Cython wrapper used `pe2+pe2` to index into the coupling matrix arrays instead of `pe1+pe2`. For example:
```python
# BUG: pe2+pe2 ignores pe1 entirely
cm['EE']['EE'][l2][l3] = fac * xi_pp[0][pe2+pe2][l2][l3]

# FIX: pe1+pe2 correctly combines both fields' purification
cm['EE']['EE'][l2][l3] = fac * xi_pp[0][pe1+pe2][l2][l3]
```

This affected all spin-spin coupling blocks (EE, BB, EB, BE).

**Impact:** Cross-spectrum estimates between a purified and non-purified field would use the wrong coupling matrix, producing incorrect debiased power spectra.

### 2. Process crash on invalid `lmax` in `compute_mcm()` (`libmcm.cpp`)

**Symptom:** Calling `compute_mcm(clmask, lmax=N)` where `N > len(clmask)-1` caused an immediate process abort (`terminate called after throwing an instance of 'int'`).

**Root cause:** The C++ `throw` statement used C's comma operator:
```cpp
// BUG: comma operator evaluates both sides, throws the int `nfin`
throw("Output array is too small %d\n", nfin);
```
This threw an `int` value instead of a string exception, which could not be caught by Python.

**Fix:**
1. Added input validation in `mcm.pyx` to raise a Python `ValueError` before reaching C++.
2. Replaced C++ `throw` with `throw std::runtime_error(...)`.

### 3. Input array mutation in `compute_mcm()` (`mcm.pyx`)

**Symptom:** After calling `compute_mcm(clmask, ...)`, the `clmask` array was permanently modified. A second call with the same array would produce different (wrong) results.

**Root cause:** The normalization `clmask *= (2*ll+1) / (4*np.pi)` modified the input array in-place.

**Fix:** Create a new array for the normalized values:
```python
clmask_norm = clmask * (2*ll+1) / (4*np.pi)
```

## Medium Bugs

### 4. `map2alm_pure()` silently accepts invalid input shapes (`purealm.py`)

**Symptom:** Passing a 1D array to `map2alm_pure()` produced garbage results instead of raising an error. `maps[1]` on a 1D array returns a single float (the pixel at index 1), which then gets broadcast in arithmetic operations, producing a valid-looking but completely wrong result.

**Fix:** Added input validation:
```python
maps = np.asarray(maps)
if maps.ndim != 2 or maps.shape[0] != 3:
    raise ValueError("maps must have shape (3, Npix) containing T, Q, U maps")
```

### 5. Missing C++ `#include` directives (`libmcm.cpp`)

**Symptom:** The C++ file relied on Cython's generated code to transitively provide necessary headers (`<vector>`, `<cmath>`, etc.). This is fragile and could break with different Cython versions or compiler configurations.

**Fix:** Added explicit includes:
```cpp
#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
```

### 6. Uninitialized variables in `drc3jj()` (`libmcm.cpp`)

**Symptom:** Compiler warnings about `denom`, `c1old`, `sumfor`, `sumbac`, and other variables potentially being used uninitialized. While the code logic ensures they are set before use in practice, the compiler cannot prove this statically.

**Fix:** Initialized all local variables at declaration.

## Low-Severity Bugs

### 7. `numpy.sqrt` RuntimeWarning in `get_spinned_windows()` (`purealm.py`)

**Symptom:** Every call to `get_spinned_windows()` or `map2alm_pure()` emitted:
```
RuntimeWarning: invalid value encountered in sqrt
```

**Root cause:** `filter_2 = -np.sqrt((ell-1.)*(ell+2.))` computes `sqrt(-2)` for `ell=0`, producing NaN (which was then overwritten with 0).

**Fix:** Compute filters only for valid ell ranges:
```python
filter_1 = np.zeros(lmax+1)
filter_2 = np.zeros(lmax+1)
filter_1[1:] = -np.sqrt((ell[1:]+1.)*ell[1:])
filter_2[2:] = -np.sqrt((ell[2:]-1.)*(ell[2:]+2.))
```

### 8. Deprecated `numpy.distutils` in `setup.py`

**Symptom:** Build produced a deprecation warning. `numpy.distutils` was removed in NumPy 2.0+ for Python >= 3.12.

**Fix:** Replaced with standard `setuptools.Extension` and added the `pur/` directory to `include_dirs` so the C++ `#include "libmcm.cpp"` directive works reliably.
