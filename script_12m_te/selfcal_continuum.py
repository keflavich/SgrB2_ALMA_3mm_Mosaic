"""
Interesting plotcal commands:

    plotcal(caltable='ampphase_3.cal', xaxis='time', yaxis='phase', iteration='field', subplot=331)
    plotcal(caltable='ampphase_3.cal', xaxis='time', yaxis='amp', iteration='field', subplot=331)
    plotcal(caltable='ampphase_3.cal', xaxis='phase', yaxis='amp', iteration='field', subplot=331)
    plotcal(caltable='ampphase_3.cal', xaxis='snr', yaxis='amp', iteration='field', subplot=331)
    plotcal(caltable='ampphase_3.cal', xaxis='time', yaxis='amp', antenna='DV18')

ANTENNA DV18 MAY BE BAD: it has a bunch of weird gain amplitudes

plotms(vis='SgrB2_TE_contsplit.ms', iteraxis='antenna', gridrows=3, gridcols=3, xaxis='uvdist', yaxis='amp')
"""
import glob
import selfcal_heuristics

contvis='SgrB2_TE_contsplit.ms'
selfcal0vis = 'selfcal_SgrB2_TE_full_selfcal_iter0.ms'
phasecenter='J2000 17:47:19.242 -28.23.33.22'
imsize=[4096,4096]

rmtables(selfcal0vis)
split(vis=contvis, outputvis=selfcal0vis,
      datacolumn='data')

outname = 'SgrB2_selfcal_full_TE'

rmtables(glob.glob(outname+".*"))
myimagebase = outname
tclean(vis=selfcal0vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="10.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=imsize,
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

rmtables(['phase_0.cal'])
gaincal(vis=selfcal0vis, caltable='phase_0.cal', solint='int', gaintype='G',
        calmode='p')

okfields, not_ok_fields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

# plotcal(caltable='phase_0.cal', xaxis='time', yaxis='phase', iteration='antenna', field=okfields_str, subplot=331)

print("Iteration 0: Self-calibrating on field IDs {0}".format(okfields_str))
# flag out the fields that do not have solutions
# (this should be overkill/redundant with not applying bad cals in the applycal step)
flagdata(vis='phase_0.cal', field=",".join(map(str, not_ok_fields)), mode='manual')

selfcal1vis = 'selfcal_SgrB2_TE_full_selfcal_iter1.ms'
rmtables([selfcal1vis])
applycal(vis=selfcal0vis, field=okfields_str, gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal0vis, outputvis=selfcal1vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TE_selfcal1'
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
       threshold="5.0mJy",
       deconvolver="clark",
       interactive=False,
       imsize=imsize,
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


rmtables(['phase_1.cal'])
gaincal(vis=selfcal1vis, caltable='phase_1.cal', solint='int', gaintype='G', calmode='p',
        field=okfields_str)

okfields2, not_ok_fields2 = selfcal_heuristics.goodenough_field_solutions('phase_1.cal')
if not okfields2:
    raise ValueError("No fields have valid solutions.")
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

#print("Iteration 1: Self-calibrating on field IDs {0}".format(okfields_str))
#flagdata(vis='phase_1.cal', field=",".join(map(str, not_ok_fields)), mode='manual')

selfcal2vis = 'selfcal_SgrB2_TE_full_selfcal_iter2.ms'
rmtables([selfcal2vis])
applycal(vis=selfcal1vis, field=okfields_str, gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal1vis, outputvis=selfcal2vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TE_selfcal2'
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
       imsize=imsize,
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



rmtables(['phase_2.cal'])
gaincal(vis=selfcal2vis, caltable='phase_2.cal', solint='int', gaintype='G', calmode='p')

#okfields, not_ok_fields = selfcal_heuristics.goodenough_field_solutions('phase_2.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

print("Iteration 2: Self-calibrating on field IDs {0}".format(okfields_str))
flagdata(vis='phase_2.cal', field=",".join(map(str, not_ok_fields)), mode='manual')

selfcal3vis = 'selfcal_SgrB2_TE_full_selfcal_iter3.ms'
rmtables([selfcal3vis])
applycal(vis=selfcal2vis, field=okfields_str, gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal2vis, outputvis=selfcal3vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TE_selfcal3'
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
       threshold="0.5mJy",
       deconvolver="clark",
       interactive=False,
       imsize=imsize,
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



rmtables(['ampphase_3.cal'])
gaincal(vis=selfcal3vis, caltable='ampphase_3.cal', solint='int', gaintype='G',
        calmode='ap', minsnr=7)
selfcal_heuristics.flag_extreme_amplitudes('ampphase_3.cal')

# use the fields that were OK before
#okfields, not_ok_fields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal4vis = 'selfcal_SgrB2_TE_full_selfcal_iter4_ampphase.ms'
rmtables([selfcal4vis])
applycal(vis=selfcal3vis, field=okfields_str, gaintable=["ampphase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal3vis, outputvis=selfcal4vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TE_selfcal4_ampphase'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal4vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="0.5mJy",
       deconvolver="clark",
       interactive=False,
       imsize=imsize,
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


rmtables(['ampphase_4.cal'])
gaincal(vis=selfcal3vis, caltable='ampphase_4.cal', solint='int', gaintype='G',
        calmode='ap')

selfcal_heuristics.flag_extreme_amplitudes('ampphase_4.cal')
#okfields, not_ok_fields = selfcal_heuristics.goodenough_field_solutions('phase_3.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal5vis = 'selfcal_SgrB2_TE_full_selfcal_iter5_ampphase.ms'
rmtables([selfcal5vis])
applycal(vis=selfcal4vis, field=okfields_str, gaintable=["ampphase_4.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal4vis, outputvis=selfcal5vis, datacolumn='corrected')
# could try selfcal5vis....


outname = 'SgrB2_selfcal_full_TE_selfcal5_ampphase'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal5vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="0.5mJy",
       deconvolver="clark",
       interactive=False,
       imsize=imsize,
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


# test to see if multiscale is any good
outname = 'SgrB2_selfcal_full_TE_selfcal4_ampphase_multiscale'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal4vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="0.5mJy",
       deconvolver="multiscale",
       scales=[0,4,12],
       interactive=False,
       imsize=imsize,
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

