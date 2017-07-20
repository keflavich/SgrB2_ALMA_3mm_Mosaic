import numpy as np
from astropy import units as u

def gas_depletion_law(sigma_gas0, t, c=0.3, alpha=2., k=1e-3*u.pc**2/u.Myr/u.M_sun):
    #if alpha != 2:
    #    k = k.to(u.pc**2/u.Myr/u.M_sun).value
    #    sigma_gas0 = sigma_gas0.to(u.M_sun/u.pc**2).value
    #    t = u.Quantity(t, u.Myr).value
    if alpha == 1:
        return c*sigma_gas0 * (1-np.exp(-k*t))
    beta = 1/(1.-alpha)
    t0 = k**-1 / (alpha-1) * sigma_gas0**(1-alpha)
    #print("t0 = {0}".format(t0))
    sigma_star = c * sigma_gas0 * (1 - (t/t0+1)**beta)
    return u.Quantity(sigma_star, u.M_sun/u.pc**2)

def sigma_gas_of_t(sigma_gas0, t, c=0.3, alpha=2., k=1e-3*u.pc**2/u.Myr/u.M_sun):
    #if alpha != 2:
    #    k = k.to(u.pc**2/u.Myr/u.M_sun).value
    #    sigma_gas0 = sigma_gas0.to(u.M_sun/u.pc**2).value
    #    t = u.Quantity(t, u.Myr).value
    if alpha == 1:
        return sigma_gas0 * np.exp(-k*t)
    t0 = k**-1 / (alpha-1) * sigma_gas0**(1-alpha)
    beta = 1/(1.-alpha)
    return sigma_gas0 * (t/t0+1)**beta

if __name__ == "__main__":

    import pylab as pl

    sigma_gas0 = np.logspace(1,5) * u.M_sun / u.pc**2

    pl.figure(1)
    pl.clf()
    pl.title("alpha=2")
    for time in (0.01, 0.1, 0.5, 1, 2, 4, 10)*u.Myr:
        pl.loglog(sigma_gas_of_t(sigma_gas0, time),
                  gas_depletion_law(sigma_gas0, time), label=time,
                  color='k')

    time = np.logspace(-2,1)*u.Myr
    for sigma_gas_ii in (10,100,1000,1e4)*u.Msun/u.pc**2:
        pl.loglog(sigma_gas_of_t(sigma_gas_ii, time),
                  gas_depletion_law(sigma_gas_ii, time), label=sigma_gas_ii,
                  color='k')

    #pl.legend(loc='best')
    pl.xlabel("$\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]")
    pl.ylabel("$\Sigma_{*}$ [M$_\odot$ pc$^{-2}$]")

    pl.figure(2)
    pl.clf()
    pl.title("alpha=1")
    for time in (0.01, 0.1, 0.5, 1, 2, 4, 10)*u.Myr:
        pl.loglog(sigma_gas_of_t(sigma_gas0, time, alpha=1, k=0.1/u.Myr),
                  gas_depletion_law(sigma_gas0, time, alpha=1, k=0.1/u.Myr), label=time,
                  color='b')

    time = np.logspace(-2,1)*u.Myr
    for sigma_gas_ii in (10,100,1000,1e4)*u.Msun/u.pc**2:
        pl.loglog(sigma_gas_of_t(sigma_gas_ii, time, alpha=1, k=0.1/u.Myr),
                  gas_depletion_law(sigma_gas_ii, time, alpha=1, k=0.1/u.Myr), label=sigma_gas_ii,
                  color='b')

    #pl.legend(loc='best')
    pl.xlabel("$\Sigma_{gas}$ [M$_\odot$ pc$^{-2}$]")
    pl.ylabel("$\Sigma_{*}$ [M$_\odot$ pc$^{-2}$]")
