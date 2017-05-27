import paths
import numpy as np
from astropy.table import Table,Column
from astropy import coordinates
from astropy import units as u
from latex_info import (latexdict, format_float, round_to_n,
                        strip_trailing_zeros, exp_to_tex)
from astropy.io import fits
from astropy import wcs
from files import contfilename

latexdict = latexdict.copy()

sgrb2contfile = fits.open(contfilename)

cont_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"), format='ascii.ipac')

coords = coordinates.SkyCoord(cont_tbl['RA'], cont_tbl['Dec'])
cont_tbl.remove_column('RA')
cont_tbl.remove_column('Dec')
cont_tbl.add_column(Column(name="Coordinates", data=coords))

for row in cont_tbl:
    row['name'] = row['name'][5:].replace("_"," ")

cont_tbl.rename_column('name','ID')

cont_tbl.add_column(Column(name='$S_{\\nu,max}$', data=cont_tbl['peak'].to(u.mJy/u.beam)))
cont_tbl.remove_column('peak')
cont_tbl.add_column(Column(name='$T_{B,max}$', data=cont_tbl['peak_K'].to(u.K)))
cont_tbl.remove_column('peak_K')
mywcs = wcs.WCS(sgrb2contfile[0].header)
pix_area = wcs.utils.proj_plane_pixel_area(mywcs)*u.deg**2
cont_tbl.add_column(Column(name='$S_{\\nu,tot}$',
                           data=(u.Quantity(cont_tbl['sum'], u.Jy/u.beam) /
                                 (u.Quantity(cont_tbl['beam_area'], u.sr) /
                                 pix_area) * u.beam).to(u.mJy,)))
cont_tbl.remove_column('sum')
cont_tbl.add_column(Column(name='$\sigma_{bg}$',
                           data=u.Quantity(cont_tbl['bgmad'],
                                           u.Jy/u.beam).to(u.mJy/u.beam)))
cont_tbl.remove_column('bgmad')
cont_tbl.add_column(Column(name='$\\alpha$', data=cont_tbl['alpha']))
cont_tbl.remove_column('alpha')
cont_tbl.add_column(Column(name='$E(\\alpha)$', data=cont_tbl['alphaerror']))
cont_tbl.remove_column('alphaerror')
cont_tbl.add_column(Column(name='$M_{40K}$', data=cont_tbl['peak_mass_40K']))
cont_tbl.remove_column('peak_mass_40K')
cont_tbl.add_column(Column(name='$N(\hh)_{40 K}$', data=cont_tbl['peak_col_40K']))
cont_tbl.remove_column('peak_col_40K')

cont_tbl.remove_column('color')

for colname in ['SIMBAD_ID', 'SIMBAD_OTYPE', 'Caswell_Name', 'npix',
                'beam_area', 'peak_90GHz', 'sum_90GHz', 'bgmad_90GHz',
                'peak_100GHz', 'sum_100GHz', 'bgmad_100GHz', 'Caswell_V_CH3OH',
                'Caswell_matchdistance', 'Muno_xray_ID', 'Muno_xray_Counts',
                'peak_mass_20K', 'peak_col_20K',
                'Muno_xray_matchdistance',]:
    cont_tbl.remove_column(colname)

formats = {'Coordinates': lambda x: x.to_string('hmsdms', sep=":"),
           '$S_{\\nu,max}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$S_{\\nu,tot}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$\sigma_{bg}$': lambda x: strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$\\alpha$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           '$E(\\alpha)$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.3f}'.format(round_to_n(x,2))),
           '$N(\hh)_{40 K}$': format_float,
           '$M_{40K}$': lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
           "$T_{B,max}$": lambda x: '-' if np.isnan(x) else strip_trailing_zeros('{0:0.2f}'.format(round_to_n(x,2))),
          }

# shorter-form units
cont_tbl['$S_{\\nu,max}$'].unit = 'mJy bm$^{-1}$'
cont_tbl['$\sigma_{bg}$'].unit = 'mJy bm$^{-1}$'

# caption needs to be *before* preamble.
#latexdict['caption'] = 'Continuum Source IDs and photometry'
latexdict['header_start'] = '\label{tab:photometry}'#\n\\footnotesize'
latexdict['preamble'] = '\caption{Continuum Source IDs and photometry}\n\\resizebox{\\textwidth}{!}{'
latexdict['col_align'] = 'l'*len(cont_tbl.columns)
latexdict['tabletype'] = 'table'
latexdict['tablefoot'] = ("}\par\n"
                          "The Classification column consists of three letter codes "
                          "as described in Section \\ref{sec:classification}.  "
                          "In column 1, "
                          "\\texttt{S} indicates a strong source, "
                          "\\texttt{W} indicates weak or low-confidence source. "
                          "In column 2, an \\texttt{X} indicates a match with the "
                          "\citet{Muno2009a} Chandra X-ray source catalog, while an"
                          "underscore indicates there was no match."
                          "In column 3, \\texttt{M} indicates a match with the, "
                          "\citet{Caswell2010a} Methanol Multibeam Survey "
                          "\\methanol maser catalog, while an underscore indicates "
                          "there was no match.  Finally, we include the SIMBAD "
                          "\\citep{Wenger2000a} "
                          "source object type classification if one was found."
                         )

cont_tbl.sort('$S_{\\nu,max}$')
cont_tbl[:-35:-1].write(paths.texpath("continuum_photometry.tex"),
                        formats=formats, overwrite=True, latexdict=latexdict)
