
# central pointing includes these fields:
field="3,66,67,68,77,78,79,87,88,89,99,100,101,102,111,112"

clean(vis='calibrated.ms', imagename='sourceM_spw1_nocontrejection',
      spw='1,5,9,13', field=field, mode='mfs', niter=5000, imsize=[512,512],
      cell='0.5arcsec', weighting='briggs', robust=0.5, pbcor=True,
      usescratch=True)

split(vis='calibrated.ms', outputvis='sourceM_spw1.ms', field=field, spw='1,5,9,13')
flagmanager(vis='sourceM_spw1.ms',
            mode='save',
            versionname='before_cont_flags')
flagmanager(vis='sourceM_spw1.ms', mode='restore', versionname='before_cont_flags')
flagchannels=('0:88.5723813165~88.5799495711GHz,0:88.5965509037~88.6702803517GHz,0:88.6861492726~88.7137367813GHz,0:88.7259436435~88.7271643297GHz,'
              '0:88.7291174277~88.7444980741GHz,0:88.7510897797~88.7525546032GHz,0:88.819204071~88.8260399139GHz,0:88.831166796~88.8394674624GHz,'
              '0:88.8428853838~88.8531391481GHz,0:88.8655901476~88.8697404807GHz,0:88.8704728925~88.879994245GHz,0:88.9066052047~88.9078258909GHz,'
              '0:88.9153941455~88.9239389491GHz,0:88.9651981434~88.974719496GHz,0:88.9759401822~88.9776491429GHz,0:88.9783815547~88.9803346526GHz,'
              '0:88.9852173975~88.9947387501GHz,0:89.0008421812~89.0040159654GHz,0:89.0079221613~89.009631122GHz,0:89.0413689638~89.0430779246GHz,'
              '0:89.045519297~89.047472395GHz,0:89.0601675317~89.0748157664GHz,0:89.091417099~89.118272196GHz,0:89.1187604705~89.1297466465GHz,'
              '0:89.1316997444~89.1395121363GHz,0:89.1492776261~89.1507424495GHz,0:89.1566017434~89.1590431159GHz,0:89.1605079393~89.2310636031GHz,'
              '0:89.2364346225~89.2527918179GHz,0:89.2728110719~89.2889241301GHz,0:89.2920979143~89.3355543439GHz,0:89.3455639709~89.3548411862GHz,'
              '0:89.3926824592~89.4031803607GHz,0:89.4100162036~89.4175844581GHz,0:89.4217347913~89.4239320265GHz,0:89.4346740653~89.4402892219GHz,'
              '0:89.4578671035~89.459331927GHz,0:89.4647029464~89.471294652GHz,0:89.472271201~89.4810601418GHz,0:89.4827691025~89.4927787296GHz,'
              '0:89.4959525137~89.4988821607GHz,0:89.5196338265~89.5215869244GHz,0:89.5276903556~89.5313524142GHz,0:89.532084826~89.5342820612GHz,'
              '0:89.5355027474~89.5818888239GHz,0:89.5965370586~89.5984901566GHz,0:89.6026404897~89.6524444877GHz,0:89.6551299974~89.6878443882GHz,'
              '0:89.6885767999~89.6902857606GHz,0:89.6978540152~89.7037133091GHz,0:89.710549152~89.7576676402GHz,0:89.7586441892~89.7620621106GHz,'
              '0:89.7686538162~89.7759779336GHz,0:89.7786634433~89.7862316979GHz,0:89.7925792662~89.7981944228GHz,0:89.8050302657~89.8074716382GHz,'
              '0:89.8086923244~89.8113778341GHz,0:89.825293657~89.8306646764GHz,0:89.8316412254~89.8345708723GHz,0:89.8394536172~89.8450687739GHz,'
              '0:89.8484866953~89.8497073815GHz,0:89.9019527519~89.9068354968GHz,0:89.9083003203~89.9104975555GHz,0:89.9234368295~89.9251457902GHz,'
              '0:89.9258782019~89.9385733386GHz,0:89.9393057504~89.9410147111GHz,0:89.9466298677~89.9534657106GHz,0:89.9683580825~89.9773911606GHz,'
              '0:89.9847152779~89.9898421601GHz,0:90.0113262376~90.0240213743GHz,0:90.0364723738~90.0379371973GHz,0:90.0564916279~90.0640598825GHz,'
              '0:90.0660129804~90.0701633136GHz,0:90.096530136~90.0999480575GHz,0:90.1241176447~90.1287562524GHz,0:90.1351038207~90.1387658794GHz,'
              '0:90.1434044871~90.1648885646GHz,0:90.1668416626~90.1707478585GHz,0:90.1790485248~90.1863726421GHz,0:90.2066360335~90.2122511901GHz,'
              '0:90.2266552875~90.2278759738GHz,0:90.242035934~90.2454538554GHz,0:90.2571724432~90.2588814039GHz,0:90.2747503248~90.2762151483GHz,'
              '0:90.2840275401~90.2915957947GHz,0:90.3172302054~90.3233336365GHz,0:90.3313901656~90.3365170477GHz,0:90.3443294396~90.3528742431GHz,'
              '0:90.3631280074~90.3680107523GHz,0:90.3951099865~90.404631339GHz,1:88.5723813165~88.5799495711GHz,1:88.5965509037~88.6702803517GHz,'
              '1:88.6861492726~88.7137367813GHz,1:88.7259436435~88.7271643297GHz,1:88.7291174277~88.7444980741GHz,1:88.7510897797~88.7525546032GHz,'
              '1:88.819204071~88.8260399139GHz,1:88.831166796~88.8394674624GHz,1:88.8428853838~88.8531391481GHz,1:88.8655901476~88.8697404807GHz,'
              '1:88.8704728925~88.879994245GHz,1:88.9066052047~88.9078258909GHz,1:88.9153941455~88.9239389491GHz,1:88.9651981434~88.974719496GHz,'
              '1:88.9759401822~88.9776491429GHz,1:88.9783815547~88.9803346526GHz,1:88.9852173975~88.9947387501GHz,1:89.0008421812~89.0040159654GHz,'
              '1:89.0079221613~89.009631122GHz,1:89.0413689638~89.0430779246GHz,1:89.045519297~89.047472395GHz,1:89.0601675317~89.0748157664GHz,'
              '1:89.091417099~89.118272196GHz,1:89.1187604705~89.1297466465GHz,1:89.1316997444~89.1395121363GHz,1:89.1492776261~89.1507424495GHz,'
              '1:89.1566017434~89.1590431159GHz,1:89.1605079393~89.2310636031GHz,1:89.2364346225~89.2527918179GHz,1:89.2728110719~89.2889241301GHz,'
              '1:89.2920979143~89.3355543439GHz,1:89.3455639709~89.3548411862GHz,1:89.3926824592~89.4031803607GHz,1:89.4100162036~89.4175844581GHz,'
              '1:89.4217347913~89.4239320265GHz,1:89.4346740653~89.4402892219GHz,1:89.4578671035~89.459331927GHz,1:89.4647029464~89.471294652GHz,'
              '1:89.472271201~89.4810601418GHz,1:89.4827691025~89.4927787296GHz,1:89.4959525137~89.4988821607GHz,1:89.5196338265~89.5215869244GHz,'
              '1:89.5276903556~89.5313524142GHz,1:89.532084826~89.5342820612GHz,1:89.5355027474~89.5818888239GHz,1:89.5965370586~89.5984901566GHz,'
              '1:89.6026404897~89.6524444877GHz,1:89.6551299974~89.6878443882GHz,1:89.6885767999~89.6902857606GHz,1:89.6978540152~89.7037133091GHz,'
              '1:89.710549152~89.7576676402GHz,1:89.7586441892~89.7620621106GHz,1:89.7686538162~89.7759779336GHz,1:89.7786634433~89.7862316979GHz,'
              '1:89.7925792662~89.7981944228GHz,1:89.8050302657~89.8074716382GHz,1:89.8086923244~89.8113778341GHz,1:89.825293657~89.8306646764GHz,'
              '1:89.8316412254~89.8345708723GHz,1:89.8394536172~89.8450687739GHz,1:89.8484866953~89.8497073815GHz,1:89.9019527519~89.9068354968GHz,'
              '1:89.9083003203~89.9104975555GHz,1:89.9234368295~89.9251457902GHz,1:89.9258782019~89.9385733386GHz,1:89.9393057504~89.9410147111GHz,'
              '1:89.9466298677~89.9534657106GHz,1:89.9683580825~89.9773911606GHz,1:89.9847152779~89.9898421601GHz,1:90.0113262376~90.0240213743GHz,'
              '1:90.0364723738~90.0379371973GHz,1:90.0564916279~90.0640598825GHz,1:90.0660129804~90.0701633136GHz,1:90.096530136~90.0999480575GHz,'
              '1:90.1241176447~90.1287562524GHz,1:90.1351038207~90.1387658794GHz,1:90.1434044871~90.1648885646GHz,1:90.1668416626~90.1707478585GHz,'
              '1:90.1790485248~90.1863726421GHz,1:90.2066360335~90.2122511901GHz,1:90.2266552875~90.2278759738GHz,1:90.242035934~90.2454538554GHz,'
              '1:90.2571724432~90.2588814039GHz,1:90.2747503248~90.2762151483GHz,1:90.2840275401~90.2915957947GHz,1:90.3172302054~90.3233336365GHz,'
              '1:90.3313901656~90.3365170477GHz,1:90.3443294396~90.3528742431GHz,1:90.3631280074~90.3680107523GHz,1:90.3951099865~90.404631339GHz,'
              '2:88.5723813165~88.5799495711GHz,2:88.5965509037~88.6702803517GHz,2:88.6861492726~88.7137367813GHz,2:88.7259436435~88.7271643297GHz,'
              '2:88.7291174277~88.7444980741GHz,2:88.7510897797~88.7525546032GHz,2:88.819204071~88.8260399139GHz,2:88.831166796~88.8394674624GHz,'
              '2:88.8428853838~88.8531391481GHz,2:88.8655901476~88.8697404807GHz,2:88.8704728925~88.879994245GHz,2:88.9066052047~88.9078258909GHz,'
              '2:88.9153941455~88.9239389491GHz,2:88.9651981434~88.974719496GHz,2:88.9759401822~88.9776491429GHz,2:88.9783815547~88.9803346526GHz,'
              '2:88.9852173975~88.9947387501GHz,2:89.0008421812~89.0040159654GHz,2:89.0079221613~89.009631122GHz,2:89.0413689638~89.0430779246GHz,'
              '2:89.045519297~89.047472395GHz,2:89.0601675317~89.0748157664GHz,2:89.091417099~89.118272196GHz,2:89.1187604705~89.1297466465GHz,'
              '2:89.1316997444~89.1395121363GHz,2:89.1492776261~89.1507424495GHz,2:89.1566017434~89.1590431159GHz,2:89.1605079393~89.2310636031GHz,'
              '2:89.2364346225~89.2527918179GHz,2:89.2728110719~89.2889241301GHz,2:89.2920979143~89.3355543439GHz,2:89.3455639709~89.3548411862GHz,'
              '2:89.3926824592~89.4031803607GHz,2:89.4100162036~89.4175844581GHz,2:89.4217347913~89.4239320265GHz,2:89.4346740653~89.4402892219GHz,'
              '2:89.4578671035~89.459331927GHz,2:89.4647029464~89.471294652GHz,2:89.472271201~89.4810601418GHz,2:89.4827691025~89.4927787296GHz,'
              '2:89.4959525137~89.4988821607GHz,2:89.5196338265~89.5215869244GHz,2:89.5276903556~89.5313524142GHz,2:89.532084826~89.5342820612GHz,'
              '2:89.5355027474~89.5818888239GHz,2:89.5965370586~89.5984901566GHz,2:89.6026404897~89.6524444877GHz,2:89.6551299974~89.6878443882GHz,'
              '2:89.6885767999~89.6902857606GHz,2:89.6978540152~89.7037133091GHz,2:89.710549152~89.7576676402GHz,2:89.7586441892~89.7620621106GHz,'
              '2:89.7686538162~89.7759779336GHz,2:89.7786634433~89.7862316979GHz,2:89.7925792662~89.7981944228GHz,2:89.8050302657~89.8074716382GHz,'
              '2:89.8086923244~89.8113778341GHz,2:89.825293657~89.8306646764GHz,2:89.8316412254~89.8345708723GHz,2:89.8394536172~89.8450687739GHz,'
              '2:89.8484866953~89.8497073815GHz,2:89.9019527519~89.9068354968GHz,2:89.9083003203~89.9104975555GHz,2:89.9234368295~89.9251457902GHz,'
              '2:89.9258782019~89.9385733386GHz,2:89.9393057504~89.9410147111GHz,2:89.9466298677~89.9534657106GHz,2:89.9683580825~89.9773911606GHz,'
              '2:89.9847152779~89.9898421601GHz,2:90.0113262376~90.0240213743GHz,2:90.0364723738~90.0379371973GHz,2:90.0564916279~90.0640598825GHz,'
              '2:90.0660129804~90.0701633136GHz,2:90.096530136~90.0999480575GHz,2:90.1241176447~90.1287562524GHz,2:90.1351038207~90.1387658794GHz,'
              '2:90.1434044871~90.1648885646GHz,2:90.1668416626~90.1707478585GHz,2:90.1790485248~90.1863726421GHz,2:90.2066360335~90.2122511901GHz,'
              '2:90.2266552875~90.2278759738GHz,2:90.242035934~90.2454538554GHz,2:90.2571724432~90.2588814039GHz,2:90.2747503248~90.2762151483GHz,'
              '2:90.2840275401~90.2915957947GHz,2:90.3172302054~90.3233336365GHz,2:90.3313901656~90.3365170477GHz,2:90.3443294396~90.3528742431GHz,'
              '2:90.3631280074~90.3680107523GHz,2:90.3951099865~90.404631339GHz,3:88.5723813165~88.5799495711GHz,3:88.5965509037~88.6702803517GHz,'
              '3:88.6861492726~88.7137367813GHz,3:88.7259436435~88.7271643297GHz,3:88.7291174277~88.7444980741GHz,3:88.7510897797~88.7525546032GHz,'
              '3:88.819204071~88.8260399139GHz,3:88.831166796~88.8394674624GHz,3:88.8428853838~88.8531391481GHz,3:88.8655901476~88.8697404807GHz,'
              '3:88.8704728925~88.879994245GHz,3:88.9066052047~88.9078258909GHz,3:88.9153941455~88.9239389491GHz,3:88.9651981434~88.974719496GHz,'
              '3:88.9759401822~88.9776491429GHz,3:88.9783815547~88.9803346526GHz,3:88.9852173975~88.9947387501GHz,3:89.0008421812~89.0040159654GHz,'
              '3:89.0079221613~89.009631122GHz,3:89.0413689638~89.0430779246GHz,3:89.045519297~89.047472395GHz,3:89.0601675317~89.0748157664GHz,'
              '3:89.091417099~89.118272196GHz,3:89.1187604705~89.1297466465GHz,3:89.1316997444~89.1395121363GHz,3:89.1492776261~89.1507424495GHz,'
              '3:89.1566017434~89.1590431159GHz,3:89.1605079393~89.2310636031GHz,3:89.2364346225~89.2527918179GHz,3:89.2728110719~89.2889241301GHz,'
              '3:89.2920979143~89.3355543439GHz,3:89.3455639709~89.3548411862GHz,3:89.3926824592~89.4031803607GHz,3:89.4100162036~89.4175844581GHz,'
              '3:89.4217347913~89.4239320265GHz,3:89.4346740653~89.4402892219GHz,3:89.4578671035~89.459331927GHz,3:89.4647029464~89.471294652GHz,'
              '3:89.472271201~89.4810601418GHz,3:89.4827691025~89.4927787296GHz,3:89.4959525137~89.4988821607GHz,3:89.5196338265~89.5215869244GHz,'
              '3:89.5276903556~89.5313524142GHz,3:89.532084826~89.5342820612GHz,3:89.5355027474~89.5818888239GHz,3:89.5965370586~89.5984901566GHz,'
              '3:89.6026404897~89.6524444877GHz,3:89.6551299974~89.6878443882GHz,3:89.6885767999~89.6902857606GHz,3:89.6978540152~89.7037133091GHz,'
              '3:89.710549152~89.7576676402GHz,3:89.7586441892~89.7620621106GHz,3:89.7686538162~89.7759779336GHz,3:89.7786634433~89.7862316979GHz,'
              '3:89.7925792662~89.7981944228GHz,3:89.8050302657~89.8074716382GHz,3:89.8086923244~89.8113778341GHz,3:89.825293657~89.8306646764GHz,'
              '3:89.8316412254~89.8345708723GHz,3:89.8394536172~89.8450687739GHz,3:89.8484866953~89.8497073815GHz,3:89.9019527519~89.9068354968GHz,'
              '3:89.9083003203~89.9104975555GHz,3:89.9234368295~89.9251457902GHz,3:89.9258782019~89.9385733386GHz,3:89.9393057504~89.9410147111GHz,'
              '3:89.9466298677~89.9534657106GHz,3:89.9683580825~89.9773911606GHz,3:89.9847152779~89.9898421601GHz,3:90.0113262376~90.0240213743GHz,'
              '3:90.0364723738~90.0379371973GHz,3:90.0564916279~90.0640598825GHz,3:90.0660129804~90.0701633136GHz,3:90.096530136~90.0999480575GHz,'
              '3:90.1241176447~90.1287562524GHz,3:90.1351038207~90.1387658794GHz,3:90.1434044871~90.1648885646GHz,3:90.1668416626~90.1707478585GHz,'
              '3:90.1790485248~90.1863726421GHz,3:90.2066360335~90.2122511901GHz,3:90.2266552875~90.2278759738GHz,3:90.242035934~90.2454538554GHz,'
              '3:90.2571724432~90.2588814039GHz,3:90.2747503248~90.2762151483GHz,3:90.2840275401~90.2915957947GHz,3:90.3172302054~90.3233336365GHz,'
              '3:90.3313901656~90.3365170477GHz,3:90.3443294396~90.3528742431GHz,3:90.3631280074~90.3680107523GHz,3:90.3951099865~90.404631339GHz')
flagdata(vis='sourceM_spw1.ms', mode='manual',
         spw=flagchannels, flagbackup=False)

os.system('rm -rf phase.cal')
os.system('rm -rf sourceM_spw1_*')

split(vis='sourceM_spw1.ms',
      spw='',
      outputvis='sourceM_spw1_cont.ms',
      # number of channels to average together. change to appropriate value for
      # each spectral window in contspws (use listobs to find) and make sure to
      # use the native number of channels per SPW (that is, not the number of
      # channels left after flagging any lines)
      width=[80,],
      datacolumn='data')

clean(vis='sourceM_spw1_cont.ms', imagename='sourceM_spw1_cont_clean', field='',
      spw='', mode='mfs', niter=500, imsize=[512,512], cell='0.5arcsec',
      weighting='briggs', robust=0.5, pbcor=True, usescratch=True)

os.system('rm -rf phase.cal')
gaincal(vis='sourceM_spw1_cont.ms', caltable="phase.cal", solint="30s",
        calmode="p", gaintype='G')
plotcal(caltable="phase.cal", xaxis="time", yaxis="phase", subplot=331,
        iteration="antenna", plotrange=[0,0,-30,30], markersize=5, fontsize=10.0,
        figfile="sourceM_spw1_selfcal_phase_scan.png")
applycal(vis="sourceM_spw1_cont.ms", gaintable=["phase.cal"], interp="linear")
os.system('rm -rf sourceM_spw1_cont_selfcal*')
split(vis="sourceM_spw1_cont.ms", outputvis="sourceM_spw1_cont_selfcal.ms", datacolumn='corrected')
clean(vis='sourceM_spw1_cont_selfcal.ms', imagename='sourceM_spw1_cont_selfcal_clean', field='',
      spw='', mode='mfs', niter=5000, imsize=[512,512], cell='0.5arcsec',
      weighting='briggs', robust=0.5, pbcor=True, usescratch=True)


applycal(vis='sourceM_spw1.ms', gaintable=['phase.cal'], interp='linear')
clean(vis='sourceM_spw1.ms', imagename='sourceM_spw1_selfcal_clean', field='',
      #datacolumn='corrected',
      spw='', mode='channel', niter=5000, imsize=[512,512], cell='0.5arcsec',
      weighting='briggs', robust=0.5, pbcor=True, usescratch=True)
clearcal(vis='sourceM_spw1.ms')
clean(vis='sourceM_spw1.ms', imagename='sourceM_spw1_clean_noselfcal', field='',
      #datacolumn='data',
      spw='', mode='channel', niter=5000, imsize=[512,512], cell='0.5arcsec',
      weighting='briggs', robust=0.5, pbcor=True, usescratch=True)
exportfits('sourceM_spw1_clean_noselfcal.image', 'sourceM_spw1_clean_noselfcal.image.fits', dropdeg=True, overwrite=True)
exportfits('sourceM_spw1_selfcal_clean.image', 'sourceM_spw1_selfcal_clean.image.fits', dropdeg=True, overwrite=True)
