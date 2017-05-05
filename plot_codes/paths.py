import os
root = '/Users/adam/work/sgrb2/SgrB2_ALMA_3mm_Mosaic'
data = 'data/'
spectra = os.path.join(data, 'spectra/')
regions = 'regions/'
spectra_plots = os.path.join(spectra,'pngs/')
tablepath = os.path.join(root, 'tables/')
figurepath = os.path.join(root, 'figures/')
contpath = os.path.join(root, 'FITS/continuumdata',)

molpath = os.path.join(root, '../molecules/')


def spath(x):
    return os.path.join(root,spectra,x)

def rpath(x):
    return os.path.join(root,regions,x)

def dpath(x):
    return os.path.join(root,data,x)

def Fpath(x):
    return os.path.join(root,'FITS',x)

def tmpath(x):
    return os.path.join(root,'FITS/12m',x)

def mergepath(x):
    return os.path.join(root,'FITS/merge',x)

def acapath(x):
    return os.path.join(root,'FITS/7m',x)

def sppath(x):
    return os.path.join(root,spectra_plots,x)

def mpath(x):
    return os.path.join(molpath,x)

def tpath(x):
    return os.path.join(tablepath,x)

def fpath(x):
    return os.path.join(figurepath,x)

def cpath(x):
    return os.path.join(contpath, x)

def texpath(x):
    return os.path.join(root, 'tex_cores', x)
