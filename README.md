# PySurf96

## Modelling Surface Wave Dispersion Curves

This a slim wrapper around `surf96` from _Computer programs in seismology_ by R. Hermann (http://www.eas.slu.edu/eqc/eqccps.html)

In this realisation we wrap the Fortran77 code through `f2py`, which is approximately **5x faster** than calling a subprocess.

## Installation

This packages is for Python2 and Python3
Prerequisits are numpy and a Fortran77 compiler, like GNU GCC.

```
sudo python setup.py install
```

Or through pip:

```
pip install git+https://github.com/miili/pysurf96
```

## Example

```
import numpy as np
from pysurf96 import surf96

# Thickness in km
thicknesses = np.array([5., 23., 8., 0])

vs = np.array([2, 3.6, 3.8, 3.3])
vp = vs * 1.73
rho = vp * .32 + .77

# Periods we are interested in
periods = np.linspace(1., 20., 20)

dispersion_velocities = surf96(thickness, vp, vs, rho, periods)
```

## Citation

    Herrmann, R. B. (2013) Computer programs in seismology: An evolving tool for instruction and research, Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096 
