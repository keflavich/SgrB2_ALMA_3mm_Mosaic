import numpy as np
import scipy.optimize
from astropy import units as u 
from astropy import constants

epsilon_ff = 0.005 # quoted anywhere from 0.004-0.1

# estimate sigma_vdisp based on equatio 17:
mcloud = 1e6 * u.M_sun
rcloud = 10 * u.pc
sigma_vdisp = ((constants.G*mcloud/(3*rcloud))**0.5).to(u.km/u.s)

# eqn21
def Sigma_sfr_eqn21(Sigma_gas):
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

    time = u.Quantity(time, u.Myr)
    tff_0 = tff(rho_0).to(u.Myr)
    
    def func(beta):
        return (beta + 0.5*np.sin(2*beta) - np.pi/2 * time / tff_0).value

    beta = scipy.optimize.newton(func, 5)
    #print(beta)

    rad = r_init * np.cos(beta)**2

    return (m_0 / (4/3. * np.pi * rad**3)).to(u.Da/u.cm**3)

def t_of_rho(rho, r_0, rho_0, m_0):
    def func(time):
        if isinstance(time, u.Quantity) and time.unit != u.Myr:
            time = time.value
        tt = u.Quantity(time, u.Myr)
        return (density_of_collapse(tt, r_0, rho_0, m_0) - rho)

    tff_0 = tff(rho_0).to(u.Myr).value
    return scipy.optimize.newton(func, tff_0/2., maxiter=2000) * u.Myr
