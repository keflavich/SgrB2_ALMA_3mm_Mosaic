import glob
import pyspeckit
from astroquery.splatalogue import Splatalogue
from astropy import units as u
import paths
import pylab as pl


plot_kwargs = {'color':'r', 'linestyle':'--'}
annotate_kwargs = {'color': 'r'}


spectra_to_species = {'core_173_k1': [('hcn', 'Hydrogen Cyanide', 5500),
                                      ('hnc', 'Hydrogen Isocyanide', 5500),
                                      ('cyanoacetylene', 'HC3N', 1500),
                                      #('formamide', 'NH2CHO', 500),
                                     ],
                        }
snu_min = {'core_173_k1': 0.01,
          }
velo = {'core_173_k1': 55*u.km/u.s,
       }

pl.figure(1).clf()

for target,species_list in spectra_to_species.items():
    spectra = pyspeckit.Spectra(glob.glob(paths.fspath("*{0}*fits".format(target))))

    for species_tuple in species_list:
        species_name, chemid, tmax = species_tuple

        for ii in range(4):
            cat = Splatalogue.query_lines(spectra[ii].xarr.min(),
                                          spectra[ii].xarr.max(),
                                          chemical_name=chemid,
                                          energy_max=tmax,
                                          energy_type='eu_k', noHFS=True,
                                          line_lists=['SLAIM'])
            spectra[ii].plotter(figure=pl.figure(1))
            spectra[ii].plotter.axis.set_ylim(snu_min[target],
                                              spectra[ii].plotter.axis.get_ylim()[1])
            spectra[ii].plotter.line_ids(["{0}_{1}".format(a,b)
                                          for a,b in zip(cat['Species'],
                                                         cat['Resolved QNs'])],
                                         cat['Freq-GHz']*u.GHz,
                                         velocity_offset=velo[target],
                                         plot_kwargs=plot_kwargs,
                                         annotate_kwargs=annotate_kwargs)
            spectra[ii].plotter.savefig(paths.fpath('line_id_spectra/{target}_{species_name}_{chemid}_spw{0}.pdf'
                                                    .format(ii, target=target,
                                                            chemid=chemid,
                                                            species_name=species_name,)),
                                        bbox_extra_artists=[])
