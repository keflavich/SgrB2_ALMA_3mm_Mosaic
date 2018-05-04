import numpy as np
import scipy.optimize
from astropy import units as u 
from astropy import constants

epsilon_ff = 0.005 # quoted anywhere from 0.004-0.1

# estimate sigma_vdisp based on equation 17:
mcloud = 1e6 * u.M_sun
rcloud = 10 * u.pc
sigma_vdisp = ((constants.G*mcloud/(3*rcloud))**0.5).to(u.km/u.s)

# eqn21
def Sigma_sfr_eqn21(Sigma_gas, epsilon_ff=0.005):
    return (2.31 * epsilon_ff * constants.G * Sigma_gas**2 / sigma_vdisp).to(u.M_sun/u.pc**2/u.yr)

def Pstar(rho, Tstar, m_0, r_0):
    """
    Probability of observing a star as a function of density

    Parameters
    ----------
    rho : cm^-3
        gas density
    Tstar : yr
        Observable timescale of stellar state
    tff_0 : s
        Initial free-fall timescale (e.g., tff for midplane density)
    """
    rho_0 = m_0 / (4/3 * np.pi * r_0**3)
    tff_0 = tff(rho_0)
    return Tstar / (Tstar + tff_0 - t_of_rho(rho, r_0, rho_0, m_0))

def tff(density):
    """
    """
    return ((32*constants.G*density/(3*np.pi))**-0.5).to(u.Myr)

def density_of_collapse(time, r_init, rho_0, m_0):

    rad, beta = radius_of_collapse(time, r_init, rho_0)

    return (m_0 / (4/3. * np.pi * rad**3)).to(u.Da/u.cm**3)

def radius_of_collapse(time, r_init, rho_0):

    time = u.Quantity(time, u.Myr)
    tff_0 = tff(rho_0).to(u.Myr)
    
    def func(beta):
        return (beta + 0.5*np.sin(2*beta) - np.pi/2 * time / tff_0).value

    beta = scipy.optimize.newton(func, 5)
    #print(beta)

    rad = r_init * np.cos(beta)**2

    return rad, beta

# def t_of_rho(rho, r_0, rho_0, m_0):
#     def func(time):
#         if isinstance(time, u.Quantity) and time.unit != u.Myr:
#             time = time.value
#         tt = u.Quantity(time, u.Myr)
#         return (density_of_collapse(tt, r_0, rho_0, m_0) - rho)
# 
#     tff_0 = tff(rho_0).to(u.Myr).value
#     return scipy.optimize.newton(func, tff_0/2., maxiter=2000) * u.Myr

def t_of_rho(rho, r_0, rho_0, m_0):
    """
    Analytic solution
    """

    tff_0 = tff(rho_0).to(u.Myr)

    r_eff = (4 * np.pi * m_0 / (3 * rho))**(1/3.)

    term1 = (np.arccos((r_eff/r_0)**0.5)).decompose().value
    term2 = ((1-r_eff/r_0)**0.5 * (r_eff/r_0)**0.5).decompose().value

    return ((term1 + term2) * 2 * tff_0 / np.pi).to(u.yr)


if __name__ == "__main__":

    r0 = 20*u.pc
    m0 = 1e6*u.M_sun
    #r0 = 1*u.pc
    #m0 = 1e4*u.M_sun
    rho_0 = m0/(4/3. * np.pi * r0**3)

    assert density_of_collapse(0, r0, rho_0, m0) == rho_0
    assert density_of_collapse(tff(rho_0), r0, rho_0, m0) > 1e20 * u.Da/u.cm**3

    import pylab as pl
    pl.figure(1).clf()
    pl.figure(2).clf()
    pl.figure(3).clf()
    pl.figure(4).clf()

    for Tstar in (3e4, 1e5, 3e5)*u.yr:

        times = np.logspace(-1, np.log10(tff(rho_0).value)*0.995) * u.Myr
        times = np.linspace(0.001, (tff(rho_0).value)*0.995)*u.Myr
        densities = np.logspace(1, 9, 100)*(2.8*u.Da/u.cm**3)

        #radii,beta = zip(*[radius_of_collapse(t, r0, rho_0) for t in times])
        #radii = u.Quantity(radii)
        #densities = [density_of_collapse(t, r0, rho_0, m0) for t in times]
        t_dens = u.Quantity([t_of_rho(rho, r0, rho_0, m0) for rho in densities], u.Myr)
        radii = ((m0 * (4*np.pi/3) / densities)**(1/3)).to(u.pc)

        probs = [Pstar(rho, Tstar, m0, r0) for rho in densities]


        pl.figure(1)
        pl.subplot(2,1,1)
        pl.semilogy(u.Quantity(t_dens), u.Quantity(densities))
        pl.xlabel("Time (Myr)")
        pl.ylabel("$\\rho$ cm$^{-3}$")
        pl.subplot(2,1,2)
        pl.semilogy(u.Quantity(t_dens), radii)
        pl.xlabel("Time (Myr)")
        pl.ylabel("$r$")

        pl.figure(2)
        pl.semilogx(u.Quantity(densities), probs, label='$T_*={0:0.2f}$ Myr'.format(Tstar.to(u.Myr).value))
        pl.ylabel("P($\\rho$)")
        pl.xlabel("$\\rho$ cm$^{-3}$")
        pl.legend(loc='best')

        #pl.figure(3)
        #pl.plot(times, t_dens)
        #pl.xlabel("Input time")
        #pl.ylabel("Time of density")

        pl.figure(4)
        pl.semilogx(u.Quantity(densities), t_dens)
        pl.ylabel("Time of density")
        pl.xlabel("$\\rho$ cm$^{-3}$")
