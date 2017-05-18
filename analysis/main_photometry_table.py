import paths
import numpy as np
from astropy.table import Table,Column
from astropy import coordinates
from astropy import units as u
from latex_info import (latexdict, format_float, round_to_n,
                        strip_trailing_zeros, exp_to_tex)

cont_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"), format='ascii.ipac')

coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'])
cont_tbl.remove_column('RA')
cont_tbl.remove_column('Dec')
cont_tbl.add_column(Column(name="Coordinates", data=coords))

for row in cont_tbl:
    row['name'] = row['name'][5:].replace("_"," ")

cont_tbl.rename_column('name','ID')

cont_tbl.add_column(Column(name='$S_{nu,max}$', data=cont_tbl['peak'].to(u.mJy/u.beam)))
cont_tbl.remove_column('peak')
cont_tbl.add_column(Column(name='$T_{B,max}$', data=cont_tbl['peak_K'].to(u.K)))
cont_tbl.remove_column('peak_K')
cont_tbl.add_column(Column(name='$S_{nu,tot}$',
                           data=(cont_tbl['sum']*cont_tbl['beam_area'].to(u.sr)*u.beam).to(u.mJy,
                                                                                           u.dimensionless_angles())))
cont_tbl.remove_column('sum')
cont_tbl.add_column(Column(name='$\sigma_{bg}$',
                           data=u.Quantity(cont_tbl['bgmad'],
                                           u.Jy/u.beam).to(u.mJy/u.beam)))
cont_tbl.remove_column('bgmad')
cont_tbl.add_column(Column(name='$\\alpha$', data=cont_tbl['alpha']))
cont_tbl.remove_column('alpha')
cont_tbl.add_column(Column(name='$E[\\alpha]$', data=cont_tbl['alphaerror']))
cont_tbl.remove_column('alphaerror')
cont_tbl.add_column(Column(name='$M(20 K)$', data=cont_tbl['peak_mass_20K']))
cont_tbl.remove_column('peak_mass_20K')
cont_tbl.add_column(Column(name='$N(\hh,20 K)$', data=cont_tbl['peak_col_20K']))
cont_tbl.remove_column('peak_col_20K')

cont_tbl.add_column(Column(name='Classification',
                           data=['S' if row['color'] == 'green' else 'W' for row in cont_tbl]))
cont_tbl.remove_column('color')

cont_tbl.remove_column('npix')
cont_tbl.remove_column('beam_area')
cont_tbl.remove_column('peak_90GHz')
cont_tbl.remove_column('sum_90GHz')
cont_tbl.remove_column('bgmad_90GHz')
cont_tbl.remove_column('peak_100GHz')
cont_tbl.remove_column('sum_100GHz')
cont_tbl.remove_column('bgmad_100GHz')

formats = {'Coordinates': lambda x: x.to_string('hmsdms'),
           '$S_{nu,max}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$S_{nu,tot}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$\sigma_{bg}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$\\alpha$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$E[\\alpha]$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$N(\hh,20 K)$': format_float,
           '$M(20 K)$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
          }

cont_tbl.write(paths.texpath("continuum_photometry.tex"), formats=formats, overwrite=True,
               latexdict=latexdict)
