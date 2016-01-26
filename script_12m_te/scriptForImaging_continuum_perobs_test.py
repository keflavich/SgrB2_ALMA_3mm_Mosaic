myimsize = [4096, 4096]
phasecenter = 'J2000 17:47:19.438000 -28.23.29.78000'

contvis='SgrB2_TE_contsplit.ms'

rootname = 'SgrB2_TE_continuum_spw{0}'

for spw in range(16):
    os.system('rm -rf '+rootname.format(spw)+'.cont_robust0.5_tclean.*')
    myimagebase = rootname.format(spw)+'.cont_robust0.5_tclean'
    tclean(vis=contvis,
         imagename=myimagebase,
         field='SgrB2',
         gridder='mosaic',
         phasecenter=phasecenter,
         spw="{0}".format(spw),
         specmode="mfs",
         niter=10000,
         threshold="0.60mJy",
         deconvolver="clark",
         interactive=False,
         imsize=myimsize,
         cell="0.125arcsec",
         outframe='LSRK',
         weighting="briggs",
         robust = 0.5,
         savemodel='none')
    impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

