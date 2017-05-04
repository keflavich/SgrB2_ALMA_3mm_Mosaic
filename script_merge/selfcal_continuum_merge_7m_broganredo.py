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


def makefits(myimagebase):
    impbcor(imagename=myimagebase+'.image.tt0', pbimage=myimagebase+'.pb.tt0', outfile=myimagebase+'.image.tt0.pbcor', overwrite=True) # perform PBcorr
    exportfits(imagename=myimagebase+'.image.tt0.pbcor', fitsimage=myimagebase+'.image.tt0.pbcor.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.image.tt1', fitsimage=myimagebase+'.image.tt1.fits', dropdeg=True, overwrite=True) # export the corrected image
    exportfits(imagename=myimagebase+'.pb.tt0', fitsimage=myimagebase+'.pb.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt0', fitsimage=myimagebase+'.model.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.model.tt1', fitsimage=myimagebase+'.model.tt1.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.residual.tt0', fitsimage=myimagebase+'.residual.tt0.fits', dropdeg=True, overwrite=True) # export the PB image
    exportfits(imagename=myimagebase+'.alpha', fitsimage=myimagebase+'.alpha.fits', dropdeg=True, overwrite=True)
    exportfits(imagename=myimagebase+'.alpha.error', fitsimage=myimagebase+'.alpha.error.fits', dropdeg=True, overwrite=True)



selfcal_vis = 'selfcal_SgrB2_TCTE7m_try2_fullfield.ms'

rmtables(selfcal_vis)
split(vis=contvis, outputvis=selfcal_vis,
      datacolumn='data')


outname = 'SgrB2_selfcal_full_TCTE7m_try2_init'

rmtables(glob.glob(outname+".*"))
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="10.0mJy",
       deconvolver='mtmfs',
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)

rmtables(['phase_0.cal'])
gaincal(vis=selfcal_vis, caltable='phase_0.cal', solint='int', gaintype='G',
        calmode='p')

okfields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])
print("Self-calibration fields: {0}".format(okfields_str))
assert "[" not in okfields_str

applycal(vis=selfcal_vis, field=okfields_str, gaintable=["phase_0.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal1'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="5.0mJy",
       deconvolver='mtmfs',
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)


rmtables(['phase_1.cal'])
gaincal(vis=selfcal_vis, caltable='phase_1.cal', solint='int', gaintype='G',
        calmode='p')

okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_1.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])
assert "[" not in okfields_str

applycal(vis=selfcal_vis, field=okfields_str, gaintable=["phase_1.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal2'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="4.0mJy",
       deconvolver='mtmfs',
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)



rmtables(['phase_2.cal'])
gaincal(vis=selfcal_vis, caltable='phase_2.cal', solint='int', gaintype='G', calmode='p')

okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_2.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])
assert "[" not in okfields_str

applycal(vis=selfcal_vis, field=okfields_str, gaintable=["phase_2.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal3'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="0.5mJy",
       deconvolver='mtmfs',
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)



rmtables(['ampphase_3.cal'])
gaincal(vis=selfcal_vis, caltable='ampphase_3.cal', solint='int', gaintype='G',
        calmode='ap', minsnr=7)
selfcal_heuristics.flag_extreme_amplitudes('ampphase_3.cal')

# use the fields that were OK before
#okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('phase_0.cal')
#okfields_str = ",".join(["{0}".format(x) for x in okfields])

applycal(vis=selfcal_vis, field=okfields_str, gaintable=["ampphase_3.cal"],
         interp="linear", applymode='calonly', calwt=False)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal4_ampphase'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       niter=100000,
       threshold="2.5mJy",
       deconvolver='mtmfs',
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)



rmtables(['ampphase_4.cal'])
gaincal(vis=selfcal_vis, caltable='ampphase_4.cal', solint='int', gaintype='G',
        calmode='ap')

selfcal_heuristics.flag_extreme_amplitudes('ampphase_4.cal')
okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('ampphase_4.cal')
okfields_str = ",".join(["{0}".format(x) for x in okfields])
assert "[" not in okfields_str

# for this iteration, try flagging too
applycal(vis=selfcal_vis, field=okfields_str, gaintable=["ampphase_4.cal"],
         interp="linear", applymode='calflag', calwt=False)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal5_ampphase'
os.system('rm -rf ' + outname + "*")
myimagebase = outname
tclean(vis=selfcal_vis,
       imagename=myimagebase,
       field='SgrB2',
       gridder='mosaic',
       spw="",
       phasecenter=phasecenter,
       specmode="mfs",
       deconvolver='mtmfs',
       niter=100000,
       threshold="5.0mJy",
       scales=[0,4,12,36],
       nterms=2,
       interactive=False,
       imsize=imsize,
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)


cleanimage = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal5_ampphase.image.tt0'
ia.open(cleanimage)
ia.calcmask(mask=cleanimage+" > 0.0025", name='clean_mask_2.5mJy')
ia.close()
makemask(mode='copy', inpimage=cleanimage,
         inpmask=cleanimage+":clean_mask_2.5mJy", output='clean_2.5mJy.mask',
         overwrite=True)
exportfits('clean_2.5mJy.mask', 'clean_2.5mJy.mask.fits', dropdeg=True, overwrite=True)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal5_ampphase_deeper_mask2.5mJy'
myimagebase = outname
os.system('rm -rf ' + outname + "*")
tclean(vis=selfcal_vis,
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
makefits(myimagebase)



rmtables(['ampphase_5.cal'])
gaincal(vis=selfcal_vis, caltable='ampphase_5.cal', solint='int', gaintype='G',
        calmode='ap')

selfcal_heuristics.flag_extreme_amplitudes('ampphase_5.cal')
okfields,notokfields = selfcal_heuristics.goodenough_field_solutions('ampphase_5.cal',
                                                                     minsnr=3)
okfields_str = ",".join(["{0}".format(x) for x in okfields])
print("OK fields for round 5->6: {0}".format(okfields_str))
assert "[" not in okfields_str

applycal(vis=selfcal_vis, field=okfields_str, gaintable=["ampphase_5.cal"],
         interp="linear", applymode='calflag', calwt=False)

cleanimage = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal5_ampphase_deeper_mask2.5mJy.image.tt0'
ia.open(cleanimage)
ia.calcmask(mask=cleanimage+" > 0.0015", name='clean_mask_1.5mJy')
ia.close()
makemask(mode='copy', inpimage=cleanimage,
         inpmask=cleanimage+":clean_mask_1.5mJy", output='clean_1.5mJy.mask',
         overwrite=True)
exportfits('clean_1.5mJy.mask', 'clean_1.5mJy.mask.fits', dropdeg=True, overwrite=True)


outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_deeper_mask1.5mJy'
myimagebase = outname
os.system('rm -rf ' + outname + "*")
tclean(vis=selfcal_vis,
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
       mask='clean_1.5mJy.mask',
       cell="0.125arcsec",
       outframe='LSRK',
       weighting="briggs",
       robust = 0.5,
       savemodel='modelcolumn')
makefits(myimagebase)






for robust in (2,1,0,):
    outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_robust{0}'.format(robust)
    os.system('rm -rf ' + outname + "*")
    myimagebase = outname
    tclean(vis=selfcal_vis,
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
           savemodel='none')
    makefits(myimagebase)

    # Do some tapering to see if we can recover any larger angular scales
    outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_taper1.5as_r{0}'.format(robust)
    os.system('rm -rf ' + outname + "*")
    myimagebase = outname
    tclean(vis=selfcal_vis,
           imagename=myimagebase,
           field='SgrB2',
           gridder='mosaic',
           spw="",
           phasecenter=phasecenter,
           specmode="mfs",
           niter=100000,
           threshold="4.0mJy",
           interactive=False,
           imsize=[1800,1800],
           deconvolver="mtmfs",
           scales=[0,4,12,36],
           nterms=2,
           cell="0.3arcsec",
           outframe='LSRK',
           weighting="briggs",
           robust=robust,
           uvtaper='1.5arcsec',
           savemodel='none')
    makefits(myimagebase)



for robust in (-2,-1):
    outname = 'SgrB2_selfcal_full_TCTE7m_try2_selfcal6_ampphase_robust{0}'.format(robust)
    myimagebase = outname
    os.system('rm -rf ' + outname + "*")
    tclean(vis=selfcal_vis,
           imagename=myimagebase,
           field='SgrB2',
           gridder='mosaic',
           spw="",
           phasecenter=phasecenter,
           specmode="mfs",
           deconvolver='mtmfs',
           nterms=2,
           niter=100000,
           threshold="2.5mJy",
           scales=[0,4,12],
           interactive=False,
           imsize=[5120,5120],
           cell="0.1arcsec",
           outframe='LSRK',
           weighting="briggs",
           robust = robust,
           savemodel='none')
    makefits(myimagebase)
