import numpy as num
import pysurf96


km = 1e3


def test_surf():
    '''
    c----- parameters
    c     thkm, vpm, vsm, rhom: model for dispersion calculation
    c     nlayer - I4: number of layers in the model
    c     iflsph - I4: 0 flat earth model, 1 spherical earth model
    c     iwave - I4: 1 Love wave, 2 Rayleigh wave
    c     mode - I4: ith mode of surface wave, 1 fundamental, 2 first higher, ....
    c     igr - I4: 0 phase velocity, > 0 group velocity
    c     kmax - I4: number of periods (t) for dispersion calculation
    c     t - period vector (t(NP))
    c     cg - output phase or group velocities (vector,cg(NP))
    c----- 
    '''
    th = num.zeros(100)
    th[:4] = num.array([5., 23., 8., 0])

    vs = num.zeros(100)
    vs[:4] = num.array([2.7, 3.6, 3.8, 4.4])

    vp = vs * 1.73
    rho = vp * .32 + .77

    nlayer = 4
    iflsph = 0
    iwave = 2
    mode = 1
    igr = 0
    kmax = 21

    t = num.zeros(60)
    t[:21] = num.linspace(1, 40, 21)
    result = num.zeros(60)

    pysurf96.surf96(th, vp, vs, rho, nlayer, iflsph, iwave,
                    mode, igr, kmax, t, result)

    return t, result


if __name__ == '__main__':
    t, cg = test_surf()
    for p, c in zip(t, cg):
        print(p, c)
