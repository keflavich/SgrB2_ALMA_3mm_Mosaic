from spectral_cube import SpectralCube
from astropy import units as u

for line, restfreq in (('HC3N','90979.02MHz'),
                       ('HNC','90.663574GHz'),
                       ('H41a','92034.43MHz'),
                       #('SiO','86.8469850GHz'),
                       ('CH3CN','91971.465MHz'),
                       ('HCOp','89.18853GHz'),
                       ('HCN','88.631847GHz'),
                       ('C3H2440-515','91.51601GHz'),
                       ('H2CS303-202','103040.548MHz'),
                       ('H2CO615-616','101332.993MHz'),
                       ('H2CS322-221','103039.93MHz'),
                       ('H2CS321-220','103051.867MHz'),
                       ('H2CS313212','101477.885MHz'),
                       ('CFp','102587.476MHz'),
                      ):
    cube = SpectralCube.read('SgrB2_b3_12M_TE.{0}.image.pbcor.fits'.format(line))
    cutout = cube[:, 1800:2700, 1800:2200]
    cutout.with_spectral_unit(u.km/u.s, velocity_convention='radio').write('{0}_cutout.fits'.format(line), overwrite=True)
    med = cutout.median(axis=0)
    medsub = cutout-med
    medsub.with_spectral_unit(u.km/u.s, velocity_convention='radio').write('{0}_cutout_medsub.fits'.format(line), overwrite=True)

"""
# run on shamash:
for fn in glob.glob("*12M_TE*fits"):
    print(fn)
    cube = SpectralCube.read(fn)[:,2151:2399,1851:2079]
    cube.write('cutouts/{0}'.format(fn.replace(".fits","_M_cutout.fits")), overwrite=True)
for fn in glob.glob("*12M_TE*fits"):
    print(fn)
    cube = SpectralCube.read(fn)[:,2589:2652,1975:2035]
    cube.write('cutouts/{0}'.format(fn.replace(".fits","_N_cutout.fits")))
"""
