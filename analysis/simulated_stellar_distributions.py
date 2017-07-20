import numpy as np
from imf import imf
import pylab as pl
import paths

# make a 10^4 msun cluster
cluster = imf.make_cluster(2e4, massfunc='kroupa', mmax=200)
ncluster = len(cluster)

# create a Gaussian random spatial distribution
xx, yy = np.random.randn(ncluster), np.random.randn(ncluster)

# quicklook plot
pl.figure(1)
pl.clf()
pl.plot(xx, yy, 'k,', alpha=0.5)

massive = cluster > 8
pl.plot(xx[massive], yy[massive], 'r.', alpha=0.5)


from scipy.spatial.kdtree import KDTree

data = np.array([xx, yy])
nn_number = 11
surfdenses = {}
msurfdenses = {}


pl.figure(2)
pl.clf()
mass_thresholds = np.array((0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 24.0, 32.0))
for mass_threshold in mass_thresholds:
    mask = cluster>mass_threshold
    tree = KDTree(data.T[mask,:])
    dist, idx = tree.query(tree.data, k=nn_number)

    dist_nn = dist[:,nn_number-1]

    surfdens = (nn_number-1) / (np.pi*(dist_nn**2))

    surfdenses[mass_threshold] = surfdens

    # for computing mean mass
    mmax = 200
    mass_samples = np.linspace(mass_threshold, mmax, 50000)
    y = imf.kroupa(mass_samples)
    overthresholdmean = (mass_samples*y).sum()/y.sum()
    overthresholdfraction = (imf.kroupa.m_integrate(mass_threshold, mmax)[0] /
                             imf.kroupa.m_integrate(imf.kroupa.mmin, mmax)[0])

    msurfdens = surfdens * overthresholdmean / overthresholdfraction
    msurfdenses[mass_threshold] = msurfdens

    pl.hist(msurfdens, histtype='step', label="{0}".format(mass_threshold),
            normed=True, bins=np.logspace(0,4,21), linewidth=2, alpha=0.75)
    pl.gca().set_xscale('log')

pl.legend(loc='best')

pl.figure(3)
pl.clf()
#pl.errorbar(mass_thresholds,
#            [msurfdenses[ii].mean() for ii in msurfdenses],
#            yerr=[msurfdenses[ii].std() for ii in msurfdenses],
#            linestyle='none')
pl.boxplot([msurfdenses[ii] for ii in mass_thresholds],
           positions=mass_thresholds,
           widths=0.1*mass_thresholds,
          )
pl.xlabel("Mass threshold")
pl.ylabel("Mass surface density")
pl.gca().set_xscale('log')
pl.gca().set_xlim(0.3, 40)
pl.savefig(paths.fpath("SimulatedCluster_SurfaceDensityRecoveryVsThreshold_Gaussian.png"))

# create a uniform random spatial distribution
xx, yy = np.random.rand(ncluster), np.random.rand(ncluster)
data = np.array([xx, yy])

for mass_threshold in mass_thresholds:
    mask = cluster>mass_threshold
    tree = KDTree(data.T[mask,:])
    dist, idx = tree.query(tree.data, k=nn_number)

    dist_nn = dist[:,nn_number-1]

    surfdens = (nn_number-1) / (np.pi*(dist_nn**2))

    surfdenses[mass_threshold] = surfdens

    # for computing mean mass
    mmax = 200
    mass_samples = np.linspace(mass_threshold, mmax, 50000)
    y = imf.kroupa(mass_samples)
    overthresholdmean = (mass_samples*y).sum()/y.sum()
    overthresholdfraction = (imf.kroupa.m_integrate(mass_threshold, mmax)[0] /
                             imf.kroupa.m_integrate(imf.kroupa.mmin, mmax)[0])

    msurfdens = surfdens * overthresholdmean / overthresholdfraction
    msurfdenses[mass_threshold] = msurfdens

pl.figure(4)
pl.clf()
#pl.errorbar(mass_thresholds,
#            [msurfdenses[ii].mean() for ii in msurfdenses],
#            yerr=[msurfdenses[ii].std() for ii in msurfdenses],
#            linestyle='none')
pl.boxplot([msurfdenses[ii] for ii in mass_thresholds],
           positions=mass_thresholds,
           widths=0.1*mass_thresholds,
          )
pl.xlabel("Mass threshold")
pl.ylabel("Mass surface density")
pl.gca().set_xscale('log')
pl.gca().set_xlim(0.3, 40)
pl.savefig(paths.fpath("SimulatedCluster_SurfaceDensityRecoveryVsThreshold_Uniform.png"))
