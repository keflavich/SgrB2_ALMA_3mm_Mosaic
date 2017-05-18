import numpy as np
import paths
from astropy.io import fits
import pylab as pl
from astropy import wcs
import astropy.visualization
from visualization import make_scalebar, hide_labels
from astropy.nddata import Cutout2D
from astropy import coordinates
from astropy import units as u
import warnings

warnings.filterwarnings('ignore', category=wcs.FITSFixedWarning, append=True)

ra1m,dec1m = (266.84197, -28.391225)
ra2m,dec2m = (266.82652, -28.378154)

ra1s,dec1s = (266.83741, -28.398368)
ra2s,dec2s = (266.83265, -28.393495)

init = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_init.image.pbcor.fits'))
sc1 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal1.image.pbcor.fits'))
sc2 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal2.image.pbcor.fits'))
sc3 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal3.image.pbcor.fits'))
sc4 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase.image.pbcor.fits'))
sc5 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale.image.tt0.pbcor.fits'))
sc5tt = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))
sc6tt = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal6_ampphase_taylorterms_multiscale_deeper_mask1.5mJy.image.tt0.pbcor.fits'))

rinit = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_init.residual.fits'))
rsc1 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal1.residual.fits'))
rsc2 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal2.residual.fits'))
rsc3 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal3.residual.fits'))
rsc4 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase.residual.fits'))
rsc5 = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale.residual.tt0.fits'))
rsc5tt = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.residual.tt0.fits'))
rsc6tt = fits.open(paths.mergepath('continuum/SgrB2_selfcal_full_TCTE7m_selfcal6_ampphase_taylorterms_multiscale_deeper_mask1.5mJy.residual.tt0.fits'))

ims = [init,sc1,sc2,sc3,sc4,sc5,sc5tt,sc6tt]
resids = [rinit,rsc1,rsc2,rsc3,rsc4,rsc5,rsc5tt,rsc6tt]
assert len(ims) == len(resids)
nims = len(ims)

for name, ((ra1,dec1),(ra2,dec2)),(vmin,vmax) in [
         ('SgrB2M', ((ra1m,dec1m),(ra2m,dec2m)), [-0.001, 0.1]),
         ('SgrB2S', ((ra1s,dec1s),(ra2s,dec2s)), [-0.001, 0.05]),]:

    fig = pl.figure(1, figsize=(20,6), dpi=75)
    fig.clf()

    for ii, (fh,rfh) in enumerate(zip(ims, resids)
                                 ):

        mywcs = wcs.WCS(fh[0].header)
        center = coordinates.SkyCoord((ra1+ra2)/2, (dec1+dec2)/2, frame='fk5',
                                      unit=(u.deg, u.deg))
        size = max([np.abs(ra2-center.ra.deg), np.abs(dec2-center.dec.deg)]) * 2.1 * u.deg
        cutout_im = Cutout2D(fh[0].data, position=center, size=size, wcs=mywcs)
        cutout_res = Cutout2D(rfh[0].data, position=center, size=size,
                              wcs=mywcs)

        ax = fig.add_subplot(2,nims,ii+1, projection=cutout_im.wcs)
        im = ax.imshow(cutout_im.data*1e3, cmap='gray',
                       norm=astropy.visualization.simple_norm(fh[0].data,
                                                              stretch='asinh',
                                                              min_cut=vmin*1e3,
                                                              max_cut=vmax*1e3,
                                                              asinh_a=0.001),
                       transform=ax.get_transform(cutout_im.wcs),
                       origin='lower',)
        ax2 = fig.add_subplot(2,nims,ii+nims+1, projection=cutout_im.wcs)
        im2 = ax2.imshow(cutout_res.data*1e3, cmap='gray',
                         norm=astropy.visualization.simple_norm(fh[0].data,
                                                                stretch='asinh',
                                                                min_cut=vmin*1e3,
                                                                max_cut=vmax*1e3,
                                                                asinh_a=0.001),
                         transform=ax.get_transform(cutout_res.wcs),
                         origin='lower',)

        #(x1,y1),(x2,y2) = mywcs.wcs_world2pix([[ra1,dec1]],0)[0], mywcs.wcs_world2pix([[ra2,dec2]],0)[0]
        (x1,y1),(x2,y2) = cutout_im.wcs.wcs_world2pix([[ra1,dec1]],0)[0], cutout_im.wcs.wcs_world2pix([[ra2,dec2]],0)[0]
        ax.axis([x1,x2,y1,y2])
        ax2.axis([x1,x2,y1,y2])
        #tr_fk5 = ax.get_transform("fk5")
        hide_labels(ax)
        hide_labels(ax2)
        #ax.plot([ra1,ra2], [dec1,dec2], transform=tr_fk5, marker='o', color='r', zorder=50)

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.canvas.draw()

    x1 = ax.bbox.x1 / (fig.bbox.x1-fig.bbox.x0)
    # y0 = ax.bbox.y0 / (fig.bbox.y1-fig.bbox.y0)
    y0 = ax2.bbox.y0 / (fig.bbox.y1-fig.bbox.y0)
    # single-ax version height = (ax.bbox.y1-ax.bbox.y0) / (fig.bbox.y1-fig.bbox.y0)
    height = (ax.bbox.y1-ax2.bbox.y0) / (fig.bbox.y1-fig.bbox.y0)
    cax_bbox = [x1 + 0.02, y0, 0.02, height]
    cb = pl.colorbar(mappable=im, cax=fig.add_axes(cax_bbox))
    cb.set_label('mJy/beam')


    fig.savefig(paths.fpath("selfcal_progression_TCTE7m_{0}.pdf".format(name)),
                bbox_inches='tight', dpi=300)
