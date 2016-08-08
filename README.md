A python package to calculate cosmological distances in a FRW metric.

Description
-----------

A package to calculate cosmological distances.  The workhorse is a class Cosmo
that wraps and the underlying C code.  The geometry can be non-flat.

In order to take advantage of the speed of the underlying C code, you should,
if at all possible, send numpy arrays to the methods rather than repeatedly
call the python methods.

Installation
------------

```bash
    python setup.py install --prefix=/some/path
```

Examples
--------

```python
import cosmology
c=cosmology.Cosmo()

# comoving distance to z=0.5
c.Dc(0.0, 0.5) 

# angular diameter distance between z=0.5 and z=0.9
c.Da(0.5, 0.9)

# luminosity distance between z=0.2 and a sequence of redshifts
c.Dl(0.2, [0.3, 0.4, 0.5])

# new cosmology
c=cosmology.Cosmo(H0=70.0, omega_m=0.25)

# inverse critical density for lensing, lens at 0.2 and
# source at 0.3
c.sigmacritinv(0.2, 0.3)
```

Tests
----------

```python
import cosmology
cosmology.test()
```

TODO
----
 - Add equation of state
 - Add *evolving* equation of state.
