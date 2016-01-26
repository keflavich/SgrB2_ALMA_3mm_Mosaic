import glob

contvis='SgrB2_TE_contsplit.ms'
centralpointingvis='SgrB2_TE_center.ms'
phasecenter='J2000 17:47:20.129 -28.23.04.07'


ptg_numbers = [3, 53, 54, 55, 56, 57, 65, 66, 67, 68, 69, 76, 77, 78, 79, 87,
               88, 89, 90, 91, 98, 99, 100, 101, 102, 110, 111, 112, 113]

rmtables(centralpointingvis)
split(vis=contvis, outputvis=centralpointingvis,
      field=",".join([str(x-3) for x in ptg_numbers]),
      datacolumn='data')

outname = 'SgrB2_center_TE'

rmtables(glob.glob(outname+".*"))
myimagebase = outname
tclean(vis=centralpointingvis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="5.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=[1024,1024],
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5, 
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image

gaincal(vis=centralpointingvis, caltable='phase_0.cal', solint='int', gaintype='G', calmode='p')
selfcal1vis = 'selfcal_SgrB2_TE_centralpointing_iter1.ms'
rmtables([selfcal1vis])
split(vis=centralpointingvis, outputvis=selfcal1vis, datacolumn='data')
applycal(vis=selfcal1vis, field="", gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_center_TE_selfcal1'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal1vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="1.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=[1024,1024],
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image


gaincal(vis=selfcal1vis, caltable='phase_1.cal', solint='int', gaintype='G', calmode='p')
selfcal2vis = 'selfcal_SgrB2_TE_centralpointing_iter2.ms'
rmtables([selfcal2vis])
split(vis=centralpointingvis, outputvis=selfcal2vis, datacolumn='data')
applycal(vis=selfcal2vis, field="", gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_center_TE_selfcal2'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal2vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="1.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=[1024,1024],
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image



gaincal(vis=selfcal2vis, caltable='phase_2.cal', solint='int', gaintype='G', calmode='p')
selfcal3vis = 'selfcal_SgrB2_TE_centralpointing_iter3.ms'
rmtables([selfcal3vis])
split(vis=centralpointingvis, outputvis=selfcal3vis, datacolumn='data')
applycal(vis=selfcal3vis, field="", gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_center_TE_selfcal3'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal3vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="1.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=[1024,1024],
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image
