Notes April 2016
----------------
Why are there two TC files, one of which I thought was the TE file?  Why are they huge?
There is one file that is 4x larger, yet it contains only 25% more data.  This doesn't make any sense.

TC2, the large file:
634G    ./SgrB2_a_03_TC2.calibrated.ms
-rw-r--r-- 1 ginsburg users 211G Oct 14  2015 table.f25_TSM1
-rw-r--r-- 1 ginsburg users 211G Jun 17  2015 table.f24_TSM1
-rw-r--r-- 1 ginsburg users 211G Jun 17  2015 table.f17_TSM1
OBSID0: 06-Dec-2014/14:12:48.8
OBSID1: 06-Dec-2014/15:40:44.6
OBSID2: 01-Apr-2015/07:05:28.6
OBSID3: 02-Apr-2015/08:18:04.1
  7   Telescope Observation Date    Observer       Project
  8   ALMA      [                   4.92459e+09, 4.9246e+09]keflavich      uid://A001/X10e/X119
  9   ALMA      [                   4.9246e+09, 4.9246e+09]keflavich      uid://A001/X10e/X119
 10   ALMA      [                   4.93459e+09, 4.93459e+09]keflavich      uid://A001/X10e/X119
 11   ALMA      [                   4.93468e+09, 4.93468e+09]keflavich      uid://A001/X10e/X119
 12 Data records: 3670524       Total elapsed time = 1.00907e+07 seconds
 13    Observed from   06-Dec-2014/14:12:48.8   to   02-Apr-2015/09:11:29.3 (UTC)
 813 Spectral Windows:  (16 unique spectral windows and 1 unique polarization setups)
814   SpwID  Name                           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs
815   0      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92220.885      -244.141   1875000.0  91283.5068        1  XX
816   1      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90417.113      -244.141   1875000.0  89479.7356        2  XX
817   2      ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100429.129       244.141   1875000.0 101366.5068        3  XX
818   3      ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102289.244       244.141   1875000.0 103226.6223        4  XX
819   4      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92220.858      -244.141   1875000.0  91283.4796        1  XX
820   5      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90417.088      -244.141   1875000.0  89479.7103        2  XX
821   6      ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100429.102       244.141   1875000.0 101366.4796        3  XX
822   7      ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102289.215       244.141   1875000.0 103226.5932        4  XX
823   8      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92232.282      -244.141   1875000.0  91294.9038        1  XX
824   9      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90427.699      -244.141   1875000.0  89490.3210        2  XX
825   10     ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100440.526       244.141   1875000.0 101377.9038        3  XX
826   11     ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102301.456       244.141   1875000.0 103238.8340        4  XX
827   12     ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92232.204      -244.141   1875000.0  91294.8263        1  XX
828   13     ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90427.627      -244.141   1875000.0  89490.2490        2  XX
829   14     ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100440.448       244.141   1875000.0 101377.8263        3  XX
830   15     ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102301.373       244.141   1875000.0 103238.7509        4  XX

TC, small file:
 158G    ./SgrB2_a_03_TC.calibrated.ms
-rw-r--r-- 1 ginsburg users 2.5G Jan 16 10:18 table.f21_TSM1
-rw-r--r-- 1 ginsburg users 156G Jan 16 10:18 table.f17_TSM1
OBSID0: 06-Dec-2014/14:12:48.8
OBSID1: 01-Apr-2015/07:05:28.6
OBSID2: 02-Apr-2015/08:18:04.1
  7   Telescope Observation Date    Observer       Project
  8   ALMA      [                   4.92459e+09, 4.9246e+09]keflavich      uid://A001/X10e/X119
  9   ALMA      [                   4.93459e+09, 4.93459e+09]keflavich      uid://A001/X10e/X119
 10   ALMA      [                   4.93468e+09, 4.93468e+09]keflavich      uid://A001/X10e/X119
 11 Data records: 2711384       Total elapsed time = 1.00907e+07 seconds
 12    Observed from   06-Dec-2014/14:12:48.8   to   02-Apr-2015/09:11:29.3 (UTC)
 651 Spectral Windows:  (12 unique spectral windows and 1 unique polarization setups)
652   SpwID  Name                           #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs
653   0      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92220.885      -244.141   1875000.0  91283.5068        1  XX
654   1      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90417.113      -244.141   1875000.0  89479.7356        2  XX
655   2      ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100429.129       244.141   1875000.0 101366.5068        3  XX
656   3      ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102289.244       244.141   1875000.0 103226.6223        4  XX
657   4      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92232.282      -244.141   1875000.0  91294.9038        1  XX
658   5      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90427.699      -244.141   1875000.0  89490.3210        2  XX
659   6      ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100440.526       244.141   1875000.0 101377.9038        3  XX
660   7      ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102301.456       244.141   1875000.0 103238.8340        4  XX
661   8      ALMA_RB_03#BB_1#SW-01#FULL_RES   7680   TOPO   92232.204      -244.141   1875000.0  91294.8263        1  XX
662   9      ALMA_RB_03#BB_2#SW-01#FULL_RES   7680   TOPO   90427.627      -244.141   1875000.0  89490.2490        2  XX
663   10     ALMA_RB_03#BB_3#SW-01#FULL_RES   7680   TOPO  100440.448       244.141   1875000.0 101377.8263        3  XX
664   11     ALMA_RB_03#BB_4#SW-01#FULL_RES   7680   TOPO  102301.373       244.141   1875000.0 103238.7509        4  XX
