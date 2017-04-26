import os
import glob
import selfcal_heuristics

vis_tc='SgrB2_TC_contsplit.ms'
#vis_te = 'TE/selfcal_SgrB2_TCTE_full_selfcal_iter4_ampphase.ms'
vis_te='SgrB2_TE_contsplit.ms'
vis_7m='SgrB2_ACA_contsplit.ms'

contvis = cont_merge_ms = combvis = 'SgrB2_ACA_TE_TC_contmerge.ms'
if not os.path.exists(cont_merge_ms):
    concat(vis=[vis_tc, vis_te, vis_7m], concatvis=cont_merge_ms)

phasecenter='J2000 17:47:19.242 -28.23.33.22'
imsize=[4096,4096]

outname = 'SgrB2_nocalTE_nocalTC_nocal7m_merge_continuum'
os.system('rm -rf ' + outname + "*")
myimagebase = outname


selfcal0vis = 'selfcal_SgrB2_TCTE7m_fullfield_iter0.ms'

rmtables(selfcal0vis)
split(vis=contvis, outputvis=selfcal0vis,
      datacolumn='data')

outname = 'SgrB2_selfcal_full_TCTE7m_init'

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

selfcal1vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter1.ms'
rmtables([selfcal1vis])
applycal(vis=selfcal0vis, field=okfields_str, gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal0vis, outputvis=selfcal1vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal1'
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

okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_1.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal2vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter2.ms'
rmtables([selfcal2vis])
applycal(vis=selfcal1vis, field=okfields_str, gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal1vis, outputvis=selfcal2vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal2'
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

okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_2.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal3vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter3.ms'
rmtables([selfcal3vis])
applycal(vis=selfcal2vis, field=okfields_str, gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal2vis, outputvis=selfcal3vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal3'
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
#okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal4vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter4_ampphase.ms'
rmtables([selfcal4vis])
applycal(vis=selfcal3vis, field=okfields_str, gaintable=["ampphase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal3vis, outputvis=selfcal4vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase'
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


# rmtables(['ampphase_4.cal'])
# gaincal(vis=selfcal3vis, caltable='ampphase_4.cal', solint='int', gaintype='G',
#         calmode='ap')
# 
# selfcal_heuristics.flag_extreme_amplitudes('ampphase_4.cal')
# #okfields = selfcal_heuristics.goodenough_field_solutions('phase_3.cal')
# #okfields_str = ",".join(["{0}".format(x) for x in okfields])
# 
# selfcal5vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter5_ampphase.ms'
# rmtables([selfcal5vis])
# applycal(vis=selfcal4vis, field=okfields_str, gaintable=["ampphase_4.cal"],
#          interp="linear", applymode='calonly', calwt=False)
# split(vis=selfcal4vis, outputvis=selfcal5vis, datacolumn='corrected')
# # could try selfcal5vis....


# test to see if multiscale is any good
# (it's not)
outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase_multiscale'
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
       scales=[0,4,12,36],
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

# taylor terms
outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase_taylorterms'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal4vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver="mtmfs",
       niter=100000,
       threshold="0.5mJy",
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image


# Do some tapering
selfcal4vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter4_ampphase.ms'

outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase_taper1.5as_r0.5'
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
       interactive=False,
       imsize=[1800,1800],
       deconvolver="multiscale",
       scales=[0,4,12,36],
       cell="0.3arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust=0.5,
       uvtaper='1.5arcsec',
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image', pbimage=myimagebase+'.pb', outfile=myimagebase+'.image.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.pbcor', fitsimage=myimagebase+'.image.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb', fitsimage=myimagebase+'.pb.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model', fitsimage=myimagebase+'.model.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual', fitsimage=myimagebase+'.residual.fits', dropdeg=True, overwrite=True) # export the PB image


#
# test to see if multiscale is any good
outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase_multiscale'
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
       scales=[0,4,12,36],
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


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal4_ampphase_taylorterms_multiscale'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal4vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver="mtmfs",
       niter=100000,
       threshold="5.0mJy",
       scales=[0,4,12,36],
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image



rmtables(['ampphase_4.cal'])
gaincal(vis=selfcal4vis, caltable='ampphase_4.cal', solint='int', gaintype='G',
        calmode='ap')

selfcal_heuristics.flag_extreme_amplitudes('ampphase_4.cal')
okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('ampphase_4.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])

selfcal5vis = 'selfcal_SgrB2_TCTE7m_full_selfcal_iter5_ampphase.ms'
rmtables([selfcal5vis])
applycal(vis=selfcal4vis, field=okfields_str, gaintable=["ampphase_4.cal"],
         interp="linear", applymode='calonly', calwt=False)
split(vis=selfcal4vis, outputvis=selfcal5vis, datacolumn='corrected')


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal5vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver='mtmfs',
       niter=100000,
       threshold="5.0mJy",
       scales=[0,4,12],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image




cleanimage = 'SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase.image'
ia.open(cleanimage)
ia.calcmask(mask=cleanimage+" > 0.0025", name='clean_mask_2.5mJy')
ia.close()
makemask(mode='copy', inpimage=cleanimage,
         inpmask=cleanimage+":clean_mask_2.5mJy", output='clean_2.5mJy.mask',
         overwrite=True)
exportfits('clean_2.5mJy.mask', 'clean_2.5mJy.mask.fits', dropdeg=True, overwrite=True)


outname = 'SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal5vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver='mtmfs',
       niter=100000,
       scales=[0,4,12,36],
       threshold="5.0mJy",
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image

outname = 'SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal5vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver='mtmfs',
       niter=100000,
       scales=[0,4,12,36],
       threshold="0.1mJy",
       nterms=2,
       interactive=False,
       imsize=imsize,
       mask='clean_2.5mJy.mask',
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image


for robust in (2,1,0,-1,-2):
    outname = 'SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_robust{0}'.format(robust)
    os.system('rm -rf ' + outname + "*")
    myimagebase = outname
    tclean(vis=selfcal5vis,
           imagename=myimagebase,
           field='SgrB2',
           gridder='mosaic',
           spw="",
           phasecenter=phasecenter,
           specmode="mfs",
           deconvolver='mtmfs',
           nterms=2,
           niter=100000,
           threshold="0.5mJy",
           scales=[0,4,12,36],
           interactive=False,
           imsize=imsize,
           cell="0.125arcsec",
           outframe='LSRK',
           weighting="briggs",
           robust = robust,
           savemodel='modelcolumn')
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
