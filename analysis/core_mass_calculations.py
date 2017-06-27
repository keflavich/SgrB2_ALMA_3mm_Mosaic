import numpy as np
from dust_emissivity import dust
from astropy import units as u
from astropy import constants
import radio_beam
from constants import frequency, distance


def snu_dust_of_density(density=1e4*u.cm**-3, Td=40*u.K, radius=4000*u.au,
                        distance=distance, cfreq=frequency):
    # from the HII region notebook; not used here...
    mass = (density * 2.8 * u.Da * 4/3. * radius**3).to(u.M_sun)
    beam = radio_beam.Beam((radius/distance).to(u.arcsec,u.dimensionless_angles()))
    flux = dust.snuofmass(nu=cfreq, mass=mass, beamomega=beam, temperature=Td, distance=distance)
    return flux


detection_limit = 0.5*u.mJy
volume_unresolved = ((0.5*u.arcsec/2.35 * distance)**3 * 4/3. * np.pi).to(u.cm**3, u.dimensionless_angles())
particle_mass = 2.8*u.Da

mthresh_20K_2pt5mjy = dust.massofsnu(nu=frequency, snu=2.5*u.mJy,
                                     distance=distance, temperature=20*u.K)
mthresh_20K = dust.massofsnu(nu=frequency, snu=detection_limit,
                             distance=distance, temperature=20*u.K)
mthresh_40K = dust.massofsnu(nu=frequency, snu=detection_limit,
                             distance=distance, temperature=40*u.K)
dens_20K = (mthresh_20K/particle_mass/volume_unresolved).to(u.cm**-3)
dens_20K_2pt5mjy = (mthresh_20K_2pt5mjy/particle_mass/volume_unresolved).to(u.cm**-3)
dens_40K = (mthresh_40K/particle_mass/volume_unresolved).to(u.cm**-3)
tff_20K = ((constants.G*mthresh_20K/volume_unresolved)**-0.5).to(u.yr)
tff_20K_2pt5mjy = ((constants.G*mthresh_20K_2pt5mjy/volume_unresolved)**-0.5).to(u.yr)
tff_40K = ((constants.G*mthresh_40K/volume_unresolved)**-0.5).to(u.yr)


mthresh_80K = dust.massofsnu(nu=frequency, snu=detection_limit,
                             distance=distance, temperature=80*u.K)
dens_80K = (mthresh_80K/particle_mass/volume_unresolved).to(u.cm**-3)
tff_80K = ((constants.G*mthresh_80K/volume_unresolved)**-0.5).to(u.yr)

print("Mass of detection limit (T=20 K): {0}"
      .format(mthresh_20K))
print("Density of 20K threshold mass: {0}, 10^{1:0.2f}"
      .format(dens_20K, np.log10(dens_20K.value)))
print("T_ff(20K): {0}".format(tff_20K))

print()
print("Mass of detection limit (T=40 K): {0}"
      .format(mthresh_40K))
print("Density of 40K threshold mass: {0}, 10^{1:0.2f}"
      .format(dens_40K, np.log10(dens_40K.value)))
print("T_ff(40K): {0}".format(tff_40K))

print()
print("Mass of 2.5 mJy detection limit (T=20 K): {0}"
      .format(mthresh_20K_2pt5mjy))
print("Density of 20K 2.5mjythreshold mass: {0}, 10^{1:0.2f}"
      .format(dens_20K_2pt5mjy, np.log10(dens_20K_2pt5mjy.value)))
print("T_ff(20K) 2.5mjy: {0}".format(tff_20K_2pt5mjy))

print()
print("Mass of detection limit (T=80 K): {0}"
      .format(mthresh_80K))
print("Density of 80K threshold mass: {0}, 10^{1:0.2f}"
      .format(dens_80K, np.log10(dens_80K.value)))
print("T_ff(80K): {0}".format(tff_80K))

