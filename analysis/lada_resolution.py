# resolution of the extinction measurements used in Lada, Lombardi, Alves 2010
from astropy import units as u
from constants import distance as d_sgrb2

lada_resolution = {'orionA':(371*u.pc, 0.05*u.deg),
                   'orionB':(398*u.pc, 0.05*u.deg),
                   'california':(450*u.pc, 0.0222*u.deg),
                   'perseus':(240*u.pc, 0.0416*u.deg),
                   'taurus':(153*u.pc, 0.0418*u.deg),
                   'ophiucus':(119*u.pc, 0.0332*u.deg),
                   'RCrA': (148*u.pc, 0.05*u.deg),
                   'Lupus 3': (230*u.pc, 0.0222*u.deg),
                   'Lupus 4': (162*u.pc, 0.0222*u.deg),
                   'Lupus 1': (144*u.pc, 0.0222*u.deg),
                   'pipe':(130*u.pc, 0.0332*u.deg),
}


for region,(distance,angres) in lada_resolution.items():
    print("{0:12s}: {1:0.2f}, or {2:0.2f} at {3}"
          .format(region, (distance*angres).to(u.pc, u.dimensionless_angles()),
                  (distance*angres/d_sgrb2).to(u.arcsec)), d_sgrb2)
