# PySurf96

_Modelling Surface Wave Dispersion Curves_

This is a slim wrapper around the program `surf96` from _Computer programs in seismology_ by R. Hermann (http://www.eas.slu.edu/eqc/eqccps.html) for forward modelling of Rayleigh and Love wave dispersion curves.

In this realisation the Fortran77 code is wrapped by `f2py`, which makes the forward computation approximately **5x faster** compared over calling a Python subprocess.

For more useful software in seismology, visit https://pyrocko.org.

## Installation

This package is for Python 2 and Python 3.

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

## Citations and Acknowledgments

> Herrmann, R. B. (2013) Computer programs in seismology: An evolving tool for instruction and research, Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096

Thanks to Hongjian Fang for creating the Fortran subroutine (https://github.com/caiweicaiwei/SurfTomo)

