__rethrow_casa_exceptions = True
h_init()
try:
    hif_restoredata (vis=['uid___A002_X95e355_X1f13', 'uid___A002_X95e355_X220a', 'uid___A002_X9cffbd_Xefe', 'uid___A002_X9d13e3_Xd4f'], session=['session_1', 'session_1', 'session_1', 'session_1'])
finally:
    h_save()
