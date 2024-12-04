from __future__ import annotations

from typing import Literal

import numpy as np

from .surfdisp96_ext import surfdisp96  # noqa

MAXLAYER = 100
MAXPERIODS = 60


WaveType = Literal["love", "rayleigh"]
Velocity = Literal["group", "phase"]


class Surf96Error(Exception):
    pass


def surf96(
    thickness: np.ndarray,
    vp: np.ndarray,
    vs: np.ndarray,
    rho: np.ndarray,
    periods: np.ndarray,
    wave: WaveType = "love",
    mode: int = 1,
    velocity: Velocity = "group",
    flat_earth: bool = True,
) -> np.ndarray:
    """Calculate synthetic surface wave dispersion curves.

    Calculate synthetic surface wave dispersion curves for a given earth model, wave
    type and periods.

    This is a slim Fortran wrapper around surf96 from Computer Programs in Seismology
    from R. Hermann (2013)

    Args:
        thickness (np.ndarray): Layer thickness in [km].
        vp (np.ndarray): Layer Vp velocity in [km/s].
        vs (np.ndarray): Layer Vs velocity in [km/s].
        rho (np.ndarray): Layer density in [g/m^3].
        periods (np.ndarray): The periods in seconds, where wave velocity is calculated
        wave (WaveType, optional): The wave type, "love" or "rayleigh".
            Defaults to "love".
        mode (int, optional): Mode of the wave, 1: fundamental, 2: second-mode, etc...
            Minimum is fundamental mode (1). Defaults to 1.
        velocity (Velocity, optional): "group" or "phase" velocity. Defaults to "group".
        flat_earth (bool, optional): Assume a flat earth. Defaults to True.

    Raises:
        ValueError: Raised when input values are unexpected.
        Surf96Error: If surf96 fortran code raises an error,
            this may be due to low velocity zone.

    Returns:
        np.ndarray: The surface wave velocities at defined periods.
    """
    if not (thickness.size == vp.size == vs.size == rho.size):
        raise ValueError("Thickness, vp/vs velocities and rho have different sizes.")
    if not (thickness.ndim == vp.ndim == vs.ndim == rho.ndim == 1):
        "Thickness, vp/vs velocities or rho have more than one dimension"
    if thickness.size > MAXLAYER:
        raise ValueError(f"Maximum number of layers is {MAXLAYER}")
    if periods.size > MAXPERIODS:
        raise ValueError(f"Maximum number of periods is {MAXPERIODS}")
    if wave not in ("love", "rayleigh"):
        raise ValueError("Wave type has to be either love or rayleigh")
    if velocity not in ("group", "phase"):
        raise ValueError("Velocity has to be group or phase")
    if mode <= 0:
        raise ValueError("Mode has to be at least 1 (fundamental mode)")

    nlayers = thickness.size
    kmax = periods.size

    _thk = np.empty(MAXLAYER)
    _vp = np.empty(MAXLAYER)
    _vs = np.empty(MAXLAYER)
    _rho = np.empty(MAXLAYER)

    _thk[:nlayers] = thickness
    _vp[:nlayers] = vp
    _vs[:nlayers] = vs
    _rho[:nlayers] = rho

    iflsph = 0 if flat_earth else 1
    iwave = 1 if wave == "love" else 2
    igr = 0 if velocity == "phase" else 1
    mode = int(mode)

    t = np.empty(MAXPERIODS)
    t[:kmax] = periods

    result = np.zeros(MAXPERIODS)

    error = surfdisp96(
        _thk, _vp, _vs, _rho, nlayers, iflsph, iwave, mode, igr, kmax, t, result
    )
    if error:
        raise Surf96Error(
            "surf96 threw an error! "
            "This may be due to low velocity zone causing"
            " reverse phase velocity dispersion,"
            " and mode jumping. Due to looking for Love waves in a halfspace"
            " which is OK if there are Rayleigh data."
        )

    return result[:kmax]


__all__ = ["surf96", "Surf96Error"]
