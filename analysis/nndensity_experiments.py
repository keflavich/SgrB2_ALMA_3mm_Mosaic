"""
Some experiments to test whether there is a bias by selecting the highest-mass
stars when computing nearest-neighbor-based density
"""

import numpy as np
import imf
from scipy import spatial
import pylab as pl


# compute some basic fractions first
kroupa = imf.Kroupa()
mmax = 200
over10fraction = (kroupa.m_integrate(10, mmax)[0] /
                 kroupa.m_integrate(kroupa.mmin, mmax)[0])
over1fraction = (kroupa.m_integrate(1, mmax)[0] /
                 kroupa.m_integrate(kroupa.mmin, mmax)[0])
# we want the mean, not the median, which is not analytic (or at least it's
# easiest to compute numerically)
x = np.linspace(10,mmax,50000)
y = kroupa(x)
over10mean = (x*y).sum()/y.sum()
x = np.linspace(1,mmax,50000)
y = kroupa(x)
over1mean = (x*y).sum()/y.sum()
x = np.linspace(kroupa.mmin,mmax,50000)
y = kroupa(x)
meanmass = (x*y).sum()/y.sum() 


nplots = 6


KDTree = spatial.cKDTree

gsz = 50

for ii,(rand, name, gridlims) in enumerate([(np.random.rand, 'Uniform', (0,1)),
                                            (np.random.randn, 'Gaussian', (-4, 4)),]): 
    gridx, gridy = np.meshgrid(np.linspace(gridlims[0], gridlims[1], gsz), np.linspace(gridlims[0], gridlims[1], gsz))
    gridxy = np.array([gridx.ravel(), gridy.ravel()]).T

    nsources = 20000

    xpos = rand(nsources)
    ypos = rand(nsources)
    masses = imf.imf.make_cluster(nsources*5)[:nsources]
    xy = np.array([xpos, ypos])

    # per-source version
    kdt = KDTree(xy.T)
    distances, indices = kdt.query(xy.T, 21)
    nn11 = distances[:,10]

    selxy1 = xy[:, masses>1]
    kdtm1 = KDTree(selxy1.T)
    distancesm1, indicesm1 = kdtm1.query(selxy1.T, 21)
    nn11m1 = distancesm1[:,10]

    selxy10 = xy[:, masses>10]
    kdtm10 = KDTree(selxy10.T)
    distancesm10, indicesm10 = kdtm10.query(selxy10.T, 21)
    nn11m10 = distancesm10[:,10]

    print(name)
    print("nn11 ({2}): {0} +/- {1}".format(nn11.mean(), nn11.std(), distances.shape[0]))
    print("nn11m1 ({2}): {0} +/- {1}".format(nn11m1.mean(), nn11m1.std(), distancesm1.shape[0]))
    print("nn11m10 ({2}): {0} +/- {1}".format(nn11m10.mean(), nn11m10.std(), distancesm10.shape[0]))

    pl.figure(6+nplots*ii).clf()
    pl.title(name)
    pl.plot(distances.mean(axis=0), label="all")
    pl.plot(distancesm1.mean(axis=0), label="M>1")
    pl.plot(distancesm10.mean(axis=0), label="M>10")
    pl.legend(loc='best')


    pl.figure(1+nplots*ii)
    pl.clf()
    pl.title("{0} random distribution".format(name))
    pl.hist(nn11, label='Full sample', normed=True, bins=25)#np.linspace(gridlims[0], gridlims[1], 20))
    pl.hist(nn11m1, label='>1 Msun', alpha=0.5, normed=True, bins=25)#np.linspace(gridlims[0], gridlims[1], 20))
    pl.hist(nn11m10, label='>10 Msun', alpha=0.5, normed=True, bins=25)#np.linspace(gridlims[0], gridlims[1], 20))
    pl.legend(loc='best')

    nn = 11
    density = (nn-1) * meanmass / (np.pi*nn11**2)
    density_m1 = (nn-1) * over1mean/over1fraction / (np.pi*nn11m1**2)
    density_m10 = (nn-1) * over10mean/over10fraction / (np.pi*nn11m10**2)

    pl.figure(2+nplots*ii)
    pl.clf()
    pl.title("{0} random distribution".format(name))
    pl.hist(density, label='Full sample', normed=True, bins=50)
    pl.hist(density_m1, label='>1 Msun', alpha=0.5, normed=True, bins=50)
    pl.hist(density_m10, label='>10 Msun', alpha=0.5, normed=True, bins=50)
    pl.legend(loc='best')
    pl.xlabel("Mass density")

    distances_g, indices_g = kdt.query(gridxy, 11)
    nn11_g = distances_g[:,10]
    nn11_gridded = nn11_g.reshape([gsz, gsz])
    #nn11_gridded_density = np.pi*nn11_gridded**2

    kdt_m1 = KDTree(xy.T[masses>1, :])
    distances_g_m1, indices_g_m1 = kdt_m1.query(gridxy, 11)
    nn11_g_m1 = distances_g_m1[:,10]
    nn11_gridded_m1 = nn11_g_m1.reshape([gsz, gsz])

    kdt_m10 = KDTree(xy.T[masses>10, :])
    distances_g_m10, indices_g_m10 = kdt_m10.query(gridxy, 11)
    nn11_g_m10 = distances_g_m10[:,10]
    nn11_gridded_m10 = nn11_g_m10.reshape([gsz, gsz])

    vmax = max([nn11_g_m10.max(), nn11_g_m1.max(), nn11_g.max()])

    pl.figure(3+nplots*ii).clf()
    pl.suptitle("{0} random distribution".format(name))
    pl.subplot(1,3,1).imshow(nn11_gridded, extent=[-2, 2, -2, 2], origin='lower', interpolation='none', vmin=0, vmax=vmax)
    pl.subplot(1,3,2).imshow(nn11_gridded_m1, extent=[-2, 2, -2, 2], origin='lower', interpolation='none', vmin=0, vmax=vmax)
    pl.subplot(1,3,3).imshow(nn11_gridded_m10, extent=[-2, 2, -2, 2], origin='lower', interpolation='none', vmin=0, vmax=vmax)

    pl.figure(4+nplots*ii).clf()
    pl.scatter(xpos, ypos, s=masses, alpha=0.4, color='k')
    pl.scatter(selxy1[0,:], selxy1[1,:], s=masses[masses>1], alpha=0.4, color='b')
    pl.scatter(selxy10[0,:], selxy10[1,:], s=masses[masses>10], alpha=0.6, color='r')

    pl.figure(5+nplots*ii).clf()
    pl.subplot(3, 1, 1).plot(masses, density, '.')
    pl.subplot(3, 1, 2).plot(masses[masses>1], density_m1, '.')
    pl.subplot(3, 1, 3).plot(masses[masses>10], density_m10, '.')


pl.figure(nplots+nplots*ii+1).clf()
#simpler test
for nsources in (100, 1000, int(1e4)):
    for rand,name,ls in [(np.random.rand, 'Uniform', '-'), (np.random.randn, 'Gaussian', '--')]:
        xpos = rand(nsources)
        ypos = rand(nsources)
        xy = np.array([xpos, ypos])

        kdt = KDTree(xy.T)
        distances, indices = kdt.query(xy.T, 31)

        pl.plot(distances.mean(axis=0), label=str(nsources)+ " " + name, linestyle=ls)

pl.legend(loc='best')
