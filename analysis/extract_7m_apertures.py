import numpy as np
from astropy import units as u
import pyspeckit
import spectral_cube
import pyregion
import os
import pylab as pl
from restfreqs import restfreqs

files = [
'SgrB2_a_03_7M.lowres.spw0.image.pbcor.fits',
'SgrB2_a_03_7M.lowres.spw1.image.pbcor.fits',
'SgrB2_a_03_7M.lowres.spw2.image.pbcor.fits',
'SgrB2_a_03_7M.lowres.spw3.image.pbcor.fits',
]

basepath = '/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic/'
regs = pyregion.open(os.path.join(basepath, "regions/spectral_apertures.reg"))

for spw,fn in enumerate(files):
    cube = spectral_cube.SpectralCube.read(os.path.join(basepath, 'FITS', fn))

    xarr = pyspeckit.units.SpectroscopicAxis(cube.spectral_axis.value,
                                             unit=str(cube.spectral_axis.unit),
                                             refX=cube.wcs.wcs.restfrq,
                                             refX_units='Hz')

    spectra = {}
    for region_number,reg in enumerate(regs):
        name = reg.attr[1]['text']
        if name not in spectra:
            #sp = cube.get_apspec(reg.coord_list,coordsys='galactic',wunit='degree')
            shape = pyregion.ShapeList([reg])
            #mask = shape.get_mask(header=noisehdr, shape=noise.shape)
            scube = cube.subcube_from_ds9region(shape)
            data = scube.apply_numpy_function(np.nanmean, axis=(1,2))
            #error = ((noise[mask & noiseokmask]**2).sum()**0.5/np.count_nonzero(mask))
            sp = pyspeckit.Spectrum(data=data,
                                    error=np.ones(data.size),#*error,
                                    xarr=xarr, header=cube.wcs.to_header())
            sp.xarr.convert_to_unit('GHz')
            #sp.header['ERROR'] = error
            #sp.error[:] = sp.stats((218.5e9,218.65e9))['std']
            sp.specname = reg.attr[1]['text']
            # Error is already computed above; this is an old hack
            #sp.error[:] = sp.stats((218e9,218.1e9))['std']
            spectra[name] = sp
            sp.unit = "$T_{A}$ [K]"
        else:
            sp = spectra[name]

        sp.plotter.figure = pl.figure(1)
        sp.plotter()
        #linesel = (lfreq > sp.xarr.as_unit('GHz').min()) & (lfreq < sp.xarr.as_unit('GHz').max())
        try:
            #sp.plotter.line_ids(names[linesel], lfreq[linesel], xval_units='GHz')
            sp.plotter.line_ids(restfreqs.keys(), [x.to(u.GHz).value for x in restfreqs.values()], xval_units='GHz')
        except RuntimeError as ex:
            print ex

        spname = sp.specname.replace(" ","_")
        sp.plotter.savefig(os.path.join(basepath, 'figures', "{0}_{1}.png".format(spname, 'spw'+str(spw))))
