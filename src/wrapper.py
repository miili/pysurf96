import numpy as num
from .surfdisp96_ext import surfdisp96  # noqa


MAXLAYER = 100
MAXPERIODS = 60


class Surf96Error(Exception):
    pass


def surf96(thickness, vp, vs, rho, periods,
           wave='love', mode=1, velocity='group', flat_earth=True):
    '''Calculate synthetic surface wave dispersion curves

    A slim Fortran wrapper around surf96 from Computer Programs in Seismology
    from R. Hermann (2013).

    Parameters
    ----------
    thickness : numpy.array
        Layer thickness in kilometers
    vp : numpy.array
        Layer Vp velocity
    vs : numpy.array
        Layer Vs velocity
    rho : numpy.array
        Layer density in g/m^3
    periods : numpy.array
        The periods in seconds, where wave velocity is calculated
    wave : str
        The wave type, "love" or "rayleigh"
    mode : int
        Mode of the wave, 1: fundamental, 2: second-mode, etc...
    velocity : str
        "group" or "phase" velocity
    flat_earth : bool
        Assume a flat earth

    Returns
    -------
    numpy.array
        The surface wave velocities at defined periods.

    Raises
    ------
    Surf96Error
        If surf96 fortran code raises an error,
        this may be due to low velocity zone.
    '''
    assert thickness.size == vp.size == vs.size == rho.size, \
        'thickness, vp/vs velocities and rho have different sizes.'
    assert thickness.ndim == vp.ndim == vs.ndim == rho.ndim == 1, \
        'thickness, vp/vs velocities or rho have more than one dimension'
    assert thickness.size <= MAXLAYER, 'maximum number of layers is 100'
    assert periods.size <= MAXPERIODS, 'maximum number of periods is 60'
    assert wave in ['love', 'rayleigh'], 'wave has to be love or rayleigh'
    assert velocity in ['group', 'phase'], 'velocity has to be group or phase'
    assert mode > 0, 'mode has to be at least 1 (fundamental mode)'

    nlayers = thickness.size
    kmax = periods.size

    _thk = num.empty(MAXLAYER)
    _vp = num.empty(MAXLAYER)
    _vs = num.empty(MAXLAYER)
    _rho = num.empty(MAXLAYER)

    _thk[:nlayers] = thickness
    _vp[:nlayers] = vp
    _vs[:nlayers] = vs
    _rho[:nlayers] = rho

    iflsph = 1 if flat_earth else 0
    iwave = 1 if wave == 'love' else 2
    igr = 0 if velocity == 'phase' else 1
    mode = int(mode)

    t = num.empty(MAXPERIODS)
    t[:kmax] = periods

    result = num.zeros(MAXPERIODS)

    error = surfdisp96(_thk, _vp, _vs, _rho, nlayers, iflsph, iwave,
                       mode, igr, kmax, t, result)
    if error > 0:
        raise Surf96Error(
            'surf96 threw an error! '
            'This may be due to low velocity zone causing'
            ' reverse phase velocity dispersion,'
            ' and mode jumping. Due to looking for Love waves in a halfspace'
            ' which is OK if there are Rayleigh data.')

    return result[:kmax]
