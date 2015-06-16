from recipes.almahelpers import fixsyscaltimes # SACM/JAO - Fixes
__rethrow_casa_exceptions = True
h_init()
try:
    hifa_importdata(vis=['uid___A002_X95e355_X1f13', 'uid___A002_X95e355_X220a', 'uid___A002_X9cffbd_Xefe', 'uid___A002_X9d13e3_Xd4f'], session=['session_1', 'session_2', 'session_3', 'session_4'])
    fixsyscaltimes(vis = 'uid___A002_X95e355_X1f13.ms')# SACM/JAO - Fixes
    fixplanets(vis = 'uid___A002_X95e355_X1f13.ms', field = '1', fixuvw = T)# SACM/JAO - Fixes
    fixsyscaltimes(vis = 'uid___A002_X95e355_X220a.ms')# SACM/JAO - Fixes
    fixplanets(vis = 'uid___A002_X95e355_X220a.ms', field = '1', fixuvw = T)# SACM/JAO - Fixes
    fixsyscaltimes(vis = 'uid___A002_X9cffbd_Xefe.ms')# SACM/JAO - Fixes
    fixsyscaltimes(vis = 'uid___A002_X9d13e3_Xd4f.ms')# SACM/JAO - Fixes
    h_save() # SACM/JAO - Finish weblog after fixes
    h_init() # SACM/JAO - Restart weblog after fixes
    hifa_importdata(vis=['uid___A002_X95e355_X1f13', 'uid___A002_X95e355_X220a', 'uid___A002_X9cffbd_Xefe', 'uid___A002_X9d13e3_Xd4f'], session=['session_1', 'session_2', 'session_3', 'session_4'])
    hifa_flagdata(pipelinemode="automatic")
    hifa_fluxcalflag(pipelinemode="automatic")
    hif_refant(pipelinemode="automatic")
    hifa_tsyscal(pipelinemode="automatic")
    hifa_tsysflag(pipelinemode="automatic")
    hifa_wvrgcalflag(pipelinemode="automatic")
    hif_lowgainflag(pipelinemode="automatic")
    hif_setjy(pipelinemode="automatic")
    hif_bandpass(pipelinemode="automatic")
    hif_bpflagchans(pipelinemode="automatic")
    hifa_gfluxscale(pipelinemode="automatic")
    hifa_timegaincal(pipelinemode="automatic")
    hif_applycal(pipelinemode="automatic")
    hif_makecleanlist(intent='PHASE,BANDPASS,CHECK')
    hif_cleanlist(pipelinemode="automatic")
finally:
    h_save()
