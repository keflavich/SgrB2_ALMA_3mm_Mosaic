import numpy as np
from astropy import constants
from sedfitter.sed import SEDCube
from astropy.table import Table
import paths
from astropy.table import Table, Column
from astropy import units as u
from astropy import coordinates
import masscalc
from constants import distance

cont_tbl = core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"),
                                      format='ascii.ipac')

highconf = core_phot_tbl['color']=='green'
lowconf = core_phot_tbl['color']=='orange'
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'

models = ['s-pbhmi', 's-pbsmi', 'sp--h-i', 's-p-hmi', 'sp--hmi', 'sp--s-i',
          's-p-smi', 'sp--smi', 'spubhmi', 'spubsmi', 'spu-hmi', 'spu-smi',
          's---s-i', 's---smi', 's-ubhmi', 's-ubsmi', 's-u-hmi', 's-u-smi']

minimum_error = 0.1

allpars = {}

for modeltype in models:
    seds = SEDCube.read('{modeltype}/flux.fits'.format(**locals()))
    pars = Table.read('{modeltype}/parameters.fits'.format(**locals()))

    luminosities = (constants.sigma_sb * (pars['star.temperature']*u.K)**4 *
                    (4*np.pi*(pars['star.radius']*u.R_sun)**2)).to(u.L_sun)

    flux_sgrb2 = seds.val[:,4,6]*(1/8.4)**2
    
    allpars[modeltype] = {}

    for row in cont_tbl:

        flux = row['peak']
        minflux,maxflux = (flux*(1-minimum_error), flux*(1+minimum_error))*u.Jy

        valid = (flux_sgrb2 > minflux) & (flux_sgrb2 < maxflux)

        vpars = pars[valid]

        allpars[modeltype][row['name']] = vpars
