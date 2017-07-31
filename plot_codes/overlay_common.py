import numpy as np
import paths
from astropy import coordinates
from astropy.table import Table

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"),
                           format='ascii.ipac')
hiis = np.array(['HII' in row['Classification'] for row in core_phot_tbl], dtype='bool')
strong = np.array(['S' in row['Classification'] for row in core_phot_tbl], dtype='bool')
weak = np.array(['W' == row['Classification'][0] for row in core_phot_tbl], dtype='bool')
xray = np.array(['X' in row['Classification'] for row in core_phot_tbl], dtype='bool')
methanolmaser = np.array(['M' in row['Classification'] for row in core_phot_tbl], dtype='bool')
watermaser = np.array(['W' == row['Classification'][5] for row in core_phot_tbl], dtype='bool')
maser = methanolmaser | watermaser
print("Found {0} HII regions, {1} strong sources, {2} weak sources, {3} X-ray sources, "
      "{4} methanol masers, and {5} water masers (total {6} masers)"
      .format(hiis.sum(), strong.sum(), weak.sum(), xray.sum(), methanolmaser.sum(),
              watermaser.sum(), maser.sum()))
assert all(weak == ~strong)
not_measured = core_phot_tbl['peak'] == 0.0
measured = ~not_measured
cores = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                             frame='fk5')
print("Found {0} HII regions, {1} strong sources, {2} weak sources, {3} X-ray sources, "
      "{4} methanol masers, and {5} water masers (total {6} masers) out of the measured data"
      .format((measured&hiis).sum(), (measured&strong).sum(),
              (measured&weak).sum(), (measured&xray).sum(),
              (measured&methanolmaser).sum(), (measured&watermaser).sum(),
              (measured&maser).sum()))

def plotcores(ax, alpha=0.5, show_unmeasured=False,
              markersize=None, markersize_override=False, **kwargs):

    all_coredots = []
    for (mask, color, marker, markersize_) in [
                                  (measured & weak & ~hiis, 'orange', 's', None),
                                  (measured & strong & ~hiis, 'red', 'o', None),
                                  (measured & hiis, 'cyan', 'o', None),
                                  #(hiis & not_measured, 'cyan', 's'),
                                  (measured & xray, 'green', 'x', None),
                                  (measured & methanolmaser, 'magenta', '+', None),
                                  (measured & watermaser, 'blue', '+', None),
    ]:

        if (not markersize_override) and (markersize_ is not None):
            markersize = markersize_

        coredots, = ax.plot(cores.ra[mask], cores.dec[mask], linestyle='none',
                            markeredgecolor=color,
                            markersize=markersize,
                            alpha=alpha, marker=marker, color=color, **kwargs)
        all_coredots.append(coredots)
    return all_coredots
