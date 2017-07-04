from spectral_cube import SpectralCube
from astropy import coordinates
import paths

from astropy import units as u

sgb2m = coordinates.SkyCoord.from_name('Sgr B2 (M)')

cube = (SpectralCube.read(paths.Fpath('merge/lines/HCN.r2_TP_7m_12m_feather.fits'))
        .with_spectral_unit(u.km/u.s, velocity_convention='radio'))

xp,yp = map(int, cube.wcs.celestial.wcs_world2pix(sgb2m.ra.deg, sgb2m.dec.deg, 0))

s19cube = SpectralCube.read(paths.Fpath('tp/tp_concat.spw19.image.fits'))
s19cube = s19cube.to(u.K, s19cube.beam.jtok_equiv(s19cube.spectral_axis))
hcncube = s19cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=cube.wcs.wcs.restfrq*u.Hz)
hcnspec = hcncube.spectral_slab(-165*u.km/u.s, 130*u.km/u.s).mean(axis=(1,2))

import pylab as pl
pl.clf()
(5*hcnspec).quicklook(color='k')
hcnsgb2spec = cube[:,yp,xp]
(hcnsgb2spec-4*hcnsgb2spec.unit).quicklook(color='k')
ax = pl.gca()
ax.set_xlim(hcnspec.spectral_axis.min().value, hcnspec.spectral_axis.max().value)
pl.savefig(paths.fpath("absorption_emission_spectrum_sgrb2.png"), bbox_inches='tight')
