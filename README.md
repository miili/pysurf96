# PySurf96

[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org/)

_Modelling Surface Wave Dispersion Curves_

This is a slim wrapper around the program `surf96` from _Computer programs in seismology_ by R. Hermann (<http://www.eas.slu.edu/eqc/eqccps.html>) for forward modelling of Rayleigh and Love wave dispersion curves.

In this implementation the Fortran77 code is wrapped by `f2py`, which makes the forward computation approximately **8x faster** compared over calling a Python subprocess.

More useful software for seismology at <https://pyrocko.org>.

## Installation

This package is for Python 3.

Prerequisits is a Fortran77 compiler, like GNU GCC.

```
pip install .
```

Or through pip:

```
pip install git+https://github.com/miili/pysurf96
```

## Documentation

Essentially this is a single function, `surf96`. Here is the docstring:

```
Calculate synthetic surface wave dispersion curves for a given earth model, wave type and periods.

This is a slim Fortran wrapper around surf96 from Computer Programs in Seismology from R. Hermann (2013)

Args:
    thickness (np.ndarray): Layer thickness in kilometers.
    vp (np.ndarray): Layer Vp velocity.
    vs (np.ndarray): Layer Vs velocity.
    rho (np.ndarray): Layer density in g/m^3.
    periods (np.ndarray): The periods in seconds, where wave velocity is calculated
    wave (WaveType, optional): The wave type, "love" or "rayleigh". Defaults to "love".
    mode (int, optional): Mode of the wave, 1: fundamental, 2: second-mode, etc... Defaults to 1.
    velocity (Velocity, optional): "group" or "phase" velocity. Defaults to "group".
    flat_earth (bool, optional): Assume a flat earth. Defaults to True.
Raises:
    ValueError: Raised when input values are unexpected.
    Surf96Error: If surf96 fortran code raises an error,
        this may be due to low velocity zone.
Returns:
    np.ndarray: The surface wave velocities at defined periods.
```

## Example

```python
import numpy as np
from pysurf96 import surf96

# Define the velocity model in km and km/s
thickness = np.array([5.0, 23.0, 8.0, 0])
vs = np.array([2, 3.6, 3.8, 3.3])
vp = vs * 1.73
rho = vp * 0.32 + 0.77

# Periods we are interested in
periods = np.linspace(1.0, 20.0, 20)

velocities = surf96(
    thickness,
    vp,
    vs,
    rho,
    periods,
    wave="love",
    mode=1,
    velocity="group",
    flat_earth=True,

```

## Citations and Acknowledgments

> Herrmann, R. B. (2013) Computer programs in seismology: An evolving tool for instruction and research, Seism. Res. Lettr. 84, 1081-1088, doi:10.1785/0220110096

Thanks to Hongjian Fang for creating the Fortran subroutine (<https://github.com/caiweicaiwei/SurfTomo>)
