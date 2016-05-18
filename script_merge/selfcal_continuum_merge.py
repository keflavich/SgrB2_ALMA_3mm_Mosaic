import os
import glob
import selfcal_heuristics

vis_tc='SgrB2_TC_contsplit.ms'
#vis_te = 'TE/selfcal_SgrB2_TCTE_full_selfcal_iter4_ampphase.ms'
vis_te='SgrB2_TE_contsplit.ms'

contvis = cont_merge_ms = "SgrB2_TETC_cont.ms"
if not os.path.exists(cont_merge_ms):
    concat(vis=[vis_tc, vis_te], concatvis=cont_merge_ms)

phasecenter='J2000 17:47:19.242 -28.23.33.22'
imsize=[4096,4096]

outname = 'SgrB2_nocalTE_nocalTC_merge_continuum'
os.system('rm -rf ' + outname + "*")
myimagebase = outname


selfcal0vis = 'selfcal_SgrB2_TCTE_fullfield_iter0.ms'

rmtables(selfcal0vis)
split(vis=contvis, outputvis=selfcal0vis,
      datacolumn='data')

outname = 'SgrB2_selfcal_full_TCTE_init'

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

okfields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])
print("Self-calibration fields: {0}".format(okfields_str))

selfcal1vis = 'selfcal_SgrB2_TCTE_full_selfcal_iter1.ms'
rmtables([selfcal1vis])
applycal(vis=selfcal0vis, field=okfields_str, gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal0vis, outputvis=selfcal1vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE_selfcal1'
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
gaincal(vis=selfcal1vis, caltable='phase_1.cal', solint='int', gaintype='G', calmode='p')

okfields = selfcal_heuristics.goodenough_field_solutions('phase_1.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal2vis = 'selfcal_SgrB2_TCTE_full_selfcal_iter2.ms'
rmtables([selfcal2vis])
applycal(vis=selfcal1vis, field=okfields_str, gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal1vis, outputvis=selfcal2vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE_selfcal2'
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

okfields = selfcal_heuristics.goodenough_field_solutions('phase_2.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal3vis = 'selfcal_SgrB2_TCTE_full_selfcal_iter3.ms'
rmtables([selfcal3vis])
applycal(vis=selfcal2vis, field=okfields_str, gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal2vis, outputvis=selfcal3vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE_selfcal3'
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
#okfields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal4vis = 'selfcal_SgrB2_TCTE_full_selfcal_iter4_ampphase.ms'
rmtables([selfcal4vis])
applycal(vis=selfcal3vis, field=okfields_str, gaintable=["ampphase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal3vis, outputvis=selfcal4vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE_selfcal4_ampphase'
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
#okfields = selfcal_heuristics.goodenough_field_solutions('phase_3.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal5vis = 'selfcal_SgrB2_TCTE_full_selfcal_iter5_ampphase.ms'
rmtables([selfcal5vis])
applycal(vis=selfcal4vis, field=okfields_str, gaintable=["ampphase_4.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal4vis, outputvis=selfcal5vis, datacolumn='corrected')
# could try selfcal5vis....


outname = 'SgrB2_selfcal_full_TCTE_selfcal5_ampphase'
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
# (it's not)
outname = 'SgrB2_selfcal_full_TCTE_selfcal4_ampphase_multiscale'
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
