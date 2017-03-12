import numpy as np
from astropy import wcs
import aplpy
import paths
from astropy.io import fits
import reproject
from matplotlib.colors import Normalize,LogNorm
from matplotlib.colors import rgb_to_hsv,hsv_to_rgb
import PIL
from PIL import ImageEnhance
import pyavm
from astropy import log

if 'fncont' not in locals():
    log.info("Loading and reprojecting FITS files")
    fncont = paths.Fpath('merge/SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase.image.pbcor.fits')
    fnhc3n = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HC3N.image.pbcor_max_medsub.fits')
    fnhcn = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HCN.image.pbcor_max_medsub.fits')
    fnhnc = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HNC.image.pbcor_max_medsub.fits')
    fnhcop = paths.Fpath('merge/max/SgrB2_b3_7M_12M.HCOp.image.pbcor_max_medsub.fits')

    outhdr = fits.getheader(fncont)
    cont = fits.getdata(fncont)
    hcn = reproject.reproject_interp(fnhcn, outhdr, order=1)[0]
    hc3n = reproject.reproject_interp(fnhc3n, outhdr, order=1)[0]
    hnc = reproject.reproject_interp(fnhnc, outhdr, order=1)[0]
    hcop = reproject.reproject_interp(fnhcop, outhdr, order=1)[0]

    log.info("Done loading and reprojecting FITS files")


# BEGIN IMAGE MAKING HERE
red,green,blue,alpha = 0,1,2,3

# FIRST LAYER: hc3n = purpley
rgb_im = np.zeros(hc3n.shape + (4,))
rgb_im[:,:,alpha]=1.0
monochrome_hc3n = Normalize(vmin=np.nanpercentile(hc3n,10),
                            vmax=np.nanpercentile(hc3n,99.0), clip=True)(hc3n)
rgb_hc3n = np.zeros(hc3n.shape + (3,))
rgb_hc3n[:,:,blue] = np.nan_to_num(monochrome_hc3n)
hsv_hc3n = rgb_to_hsv(rgb_hc3n)
hue, saturation, value = 0,1,2
hsv_hc3n[:,:,hue] = 330/360.
rgb_hc3n = hsv_to_rgb(hsv_hc3n)
rgb_im[:,:,:alpha] += np.nan_to_num(rgb_hc3n)

# SECOND LAYER: hcn = orangish?  ...
monochrome_hcn = Normalize(vmin=np.nanpercentile(hcn,10),
                           vmax=np.nanpercentile(hcn,99.995), clip=True)(hcn)
rgb_hcn = np.zeros(hc3n.shape + (3,))
rgb_hcn[:,:,blue] = np.nan_to_num(monochrome_hcn)
hsv_hcn = rgb_to_hsv(rgb_hcn)
hue, saturation, value = 0,1,2
hsv_hcn[:,:,hue] = 65/360.
rgb_hcn = hsv_to_rgb(hsv_hcn)
#rgb_im[:,:,:alpha] += rgb_hcn


monochrome_hnc = Normalize(vmin=np.nanpercentile(hnc,10),
                           vmax=np.nanpercentile(hnc,99.9995), clip=True)(hnc)
rgb_hnc = np.zeros(hc3n.shape + (3,))
#rgb_im[:,:,red] += np.nan_to_num(monochrome_hnc)
rgb_hnc[:,:,blue] = np.nan_to_num(monochrome_hnc)
hsv_hnc = rgb_to_hsv(rgb_hnc)
hsv_hnc[:,:,hue] = 25/360.
rgb_hnc = hsv_to_rgb(hsv_hnc)
rgb_im[:,:,:alpha] += rgb_hnc

monochrome_hcop = Normalize(vmin=np.nanpercentile(hcop,10),
                            vmax=np.nanpercentile(hcop,99.9995),
                            clip=True)(hcop)
rgb_hcop = np.zeros(hc3n.shape + (3,))
rgb_hcop[:,:,blue] = np.nan_to_num(monochrome_hcop)
hsv_hcop = rgb_to_hsv(rgb_hcop)
hsv_hcop[:,:,hue] = 45/360.
rgb_hcop = hsv_to_rgb(hsv_hcop)
rgb_im[:,:,:alpha] += rgb_hcop


vmin = 0.0002
monochrome_cont = LogNorm(vmin=vmin, vmax=np.nanpercentile(cont,99.9995),
                          clip=True)(cont*(cont > vmin))
monochrome_cont = monochrome_cont.filled(np.nan)
# vmin = 0.0001
# monochrome_cont = Normalize(vmin=vmin, vmax=np.nanpercentile(cont,99.9),
#                             clip=True)(cont*(cont > vmin))
#monochrome_cont[monochrome_cont.mask] = 0.0
#monochrome_cont.mask[:] = False
rgb_cont = monochrome_cont[:,:,None]
rgb_cont = np.zeros(hc3n.shape + (3,))
rgb_cont[:,:,blue] = np.nan_to_num(monochrome_cont)
hsv_cont = rgb_to_hsv(rgb_cont)
hue, saturation, value = 0,1,2
hsv_cont[:,:,hue] = 150/360.
rgb_cont = hsv_to_rgb(hsv_cont)
rgb_im[:,:,:alpha] += rgb_cont

#rgb_im[rgb_im>1] = 1
#rgb_im /= rgb_im.max()

avm = pyavm.AVM.from_header(outhdr)

im = PIL.Image.fromarray((rgb_im[:,:,:alpha]/rgb_im[:,:,:alpha].max()*255).astype('uint8')[::-1,:])
outname = paths.fpath("rgb_overview_default.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Contrast(im).enhance(1.5)
outname = paths.fpath("rgb_overview_contrast.png")
im.save(outname)
avm.embed(outname, outname)
im = ImageEnhance.Brightness(im).enhance(1.5)
outname = paths.fpath("rgb_overview_brightness.png")
im.save(outname)
avm.embed(outname, outname)
#im = ImageEnhance.Brightness(im).enhance(1.5)
#outname = paths.fpath("rgb_overview_brightness_1.5.png")
#im.save(outname)
#avm.embed(outname, outname)


FF = aplpy.FITSFigure(outname)
FF.show_rgb(outname)
#FF.show_regions(paths.rpath('overview_labels.reg'), layer='labels')
FF.save(paths.fpath("rgb_overview_aplpy_withlabels.png"), dpi=150)

