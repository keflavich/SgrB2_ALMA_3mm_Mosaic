import numpy as np
from imf import imf
import paths
import copy
from astropy.table import Table,Column
from astropy import coordinates
from astropy import units as u
from astropy import constants
import regions
import latex_info

from constants import distance

import reproject
from astropy import wcs
from astropy.io import fits
# workaround for new regions API
arbitrary_wcs = wcs.WCS(fits.getheader(paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits')))

core_phot_tbl = Table.read(paths.tpath("continuum_photometry_withSIMBAD.ipac"),
                           format='ascii.ipac')
core_phot_tbl.add_column(Column(name='radius', data=np.zeros(len(core_phot_tbl))+np.nan))

# measured from core_flux_distributions
core_powerlaw_index = 1.94

kroupa = imf.Kroupa()

o_mmin = 8
hii_cutoff = 20

mmax = 150
over8fraction150 = (kroupa.m_integrate(o_mmin, mmax)[0] /
                    kroupa.m_integrate(kroupa.mmin, mmax)[0])
x = np.linspace(o_mmin,mmax,50000)
y = kroupa(x)
over8mean150 = (x*y).sum()/y.sum()
print("150mmax,  >8 fraction: {0}, mean = {1}".format(over8fraction150, over8mean150))


mmax = 200
over8fraction = (kroupa.m_integrate(o_mmin, mmax)[0] /
                 kroupa.m_integrate(kroupa.mmin, mmax)[0])
# from wikipedia, median of power-law is well-defined for more slopes than the
# mean, but we're interested in the mean so I just compute that numerically
# below
over8median = 2**(1/(1.94-1)) * o_mmin

# we want the mean, not the median, which is not analytic (or at least it's
# easiest to compute numerically)
x = np.linspace(o_mmin,mmax,50000)
y = kroupa(x)
over8mean = (x*y).sum()/y.sum()
print("200mmax,  >8 fraction: {0}, mean = {1}".format(over8fraction, over8mean))

nsources = len(core_phot_tbl)

print("Mass fraction M>8 = {0}".format(over8fraction))
print("Mean mass Mbar(M>8) = {0}".format(over8mean))
print("Mass of observed sources, assuming all are 8 msun = {0}".format(nsources*8))
print("Total Mass estimate if all sources are 8 msun = {0}".format(nsources*8/over8fraction))
print("Total Mass estimate if Mbar={1} = {0}".format(nsources*over8mean/over8fraction, over8mean))

# Comparison to Schmiedeke et al, 2016 tbl 2
clusters = regions.read_ds9(paths.rpath('schmiedeke_clusters.reg'))
clusters.append(regions.CircleSkyRegion(clusters[0].center,
                                        radius=1*u.deg,
                                        meta={'text':'Total'})
               )
#this can be used to determine n_cores and n_hii in a bigger region,
#but is incompatible with something else below...
# clusters.append(regions.CircleSkyRegion(clusters[0].center,
#                                         radius=35*u.arcsec,
#                                         meta={'text':'M Bigger'})
#                )
# Cluster M Bigger: N(cores)= 52 N(HII)= 49 counted mass=   2852.28 inferred mass=  11522.81 HII-only inferred mass:   16003.53 core-inferred mass=   7042.09


# add in DePree HII regions based on whether or not their names
# are already in the table, since we didn't count the larger HII regions
hii_regions = regions.read_ds9(paths.rpath('SgrB2_1.3cm_hiiRegions_masked_Done.reg'))
hii_regions = Table.read(paths.tpath("Schmiedeke2016_HIIregions_tableB1.txt"), format='ascii.fixed_width')
schmiedeke_dust_regions = Table.read(paths.tpath("Schmiedeke2016_dustsources_tableB2.txt"), format='ascii.fixed_width')
schmiedeke_summary_table = Table.read(paths.tpath("Schmiedeke2016_HIIregions_table2.txt"), format='ascii.fixed_width')

def obj_in_tbl(objname):
    for ii,name in enumerate(core_phot_tbl['name']):
        if str(objname).lower() in str(name).lower():
            return ii
    return False


for row in hii_regions:
    #if 'text' not in reg.meta:
    #    continue
    #nm = reg.meta['text'].strip("{}")
    nm = row['ID']
    if (hasattr(nm,'mask') and nm.mask):
        pass
        #print("Skipping row {0} because masked.".format(row))
    elif obj_in_tbl(nm) is False:
        if row['robs'] > 5:
            #print("Skipping row {0} because it's too large".format(row))
            continue
        coord = coordinates.SkyCoord(row['RA'], row['Dec'], frame='fk5', unit=(u.hour, u.deg))
        core_phot_tbl.add_row({'name': nm, 'SIMBAD_OTYPE':'HII',
                               'RA': coord.ra,
                               'Dec': coord.dec,
                               'color': 'green',
                               'Muno_xray_ID': '-',
                               'Caswell_Name': '-',
                               'radius': u.Quantity(row['robs'], unit=1000*u.au).to(u.pc),
                               #'Classification': 'HII',
                               #'RA': reg.center.ra[0],
                               #'Dec': reg.center.dec[0]
                              })
    else:
        match = obj_in_tbl(nm)
        core_phot_tbl[match]['radius'] = u.Quantity(row['robs'], unit=1000*u.au).to(u.pc).value
        #print('{0} found in table'.format(nm))
hii = core_phot_tbl['SIMBAD_OTYPE'] == 'HII'
core_coords = coordinates.SkyCoord(core_phot_tbl['RA'], core_phot_tbl['Dec'],
                                   frame='fk5')

# Mean of "cores"
x = np.linspace(o_mmin,hii_cutoff,50000)
y = kroupa(x)
over8lt20mean = (x*y).sum()/y.sum()
over8lt20fraction = (kroupa.m_integrate(o_mmin, hii_cutoff)[0] /
                     kroupa.m_integrate(kroupa.mmin, mmax)[0])
over8lt20representative = over8lt20mean/over8lt20fraction
print("Representative mass for 20>M>8: {0}".format(over8lt20representative))
# Mean of "HII regions"
x = np.linspace(hii_cutoff,mmax,50000)
y = kroupa(x)
over20mean = (x*y).sum()/y.sum()
over20fraction = (kroupa.m_integrate(hii_cutoff, mmax)[0] /
                  kroupa.m_integrate(kroupa.mmin, mmax)[0])
over20representative = over20mean/over20fraction
print("Representative mass for M>20: {0}".format(over20representative))

tbl = Table(names=['Name', '$N({\\rm cores})$', '$N({\\rm H\\textsc{ii}})$', '$M_{\\rm count}$',
                   '$M_{\\rm inferred}$', '$M_{\\rm inferred,max}$', '$M_{\\rm inferred, H\\textsc{ii}}$',
                   '$M_{\\rm inferred, cores}$', '$M_{\\rm count}^s$', '$M_{\\rm inf}^s$',
                   'SFR'],
            dtype=['S23', int, int, int, int, int, int, int, int, int, float])

for col in tbl.colnames:
    if 'M' in col:
        tbl[col].unit = u.Msun
tbl['SFR'].unit = u.Msun/u.yr
sgrb2_age_myr = 0.74
sgrb2_brick_age_myr = 0.43

total_picking_max = 0*u.M_sun

cluster_column = np.array(['--']*len(core_phot_tbl))

print("Mass fraction M>20 = {0}".format(over20fraction))
print("Mean mass M>20 = {0}".format(over20mean))
print("Mass fraction 8<M<20 = {0}".format(over8lt20fraction))
print("Mean mass 8<M<20 = {0}".format(over8lt20mean))
allclusters = np.zeros_like(hii, dtype='bool')
NandM = np.zeros_like(hii, dtype='bool')
for reg in clusters:
    mask = reg.contains(core_coords, arbitrary_wcs)
    nhii = (hii & mask).sum()
    ncores = ((~hii) & mask).sum()

    mass = ncores * over8lt20mean + nhii * over20mean
    # the fractions are fractions-of-total, so we're estimating the total mass
    # from each independent population assuming they're the same age
    hii_only_inferred_mass = nhii * over20mean / over20fraction
    core_inferred_mass = ncores * over8lt20mean / over8lt20fraction
    inferred_mass = (core_inferred_mass +
                     hii_only_inferred_mass) / 2.


    name = reg.meta.get('text', reg.meta.get('label')).strip("{}")
    print("Cluster {0:4s}: N(cores)={1:3d} N(HII)={2:3d} counted mass={3:10.2f}"
          " inferred mass={4:10.2f} HII-only inferred mass: {5:10.2f}"
          " core-inferred mass={6:10.2f}"
          .format(name, ncores, nhii, mass,
                  inferred_mass, hii_only_inferred_mass, core_inferred_mass))

    if name == 'Total':
        sst_mask = [-1]
    else:
        sst_mask = schmiedeke_summary_table['Name'] == 'Sgr B2({0})'.format(name)
        cluster_column[mask] = name

        allclusters += mask
        print("total stars in clusters: ",allclusters.sum())
        if name in ('M','N'):
            NandM += mask
            if name == "M":
                maskM = mask
            elif name == "N":
                maskN = mask

        rowmax = max([hii_only_inferred_mass, core_inferred_mass])*u.M_sun
        total_picking_max += rowmax

    row = [name,
           ncores,
           nhii,
           latex_info.round_to_n(mass,2)*u.M_sun,
           latex_info.round_to_n(inferred_mass, 2)*u.M_sun,
           latex_info.round_to_n(rowmax.value, 2)*u.M_sun,
           latex_info.round_to_n(hii_only_inferred_mass, 2)*u.M_sun,
           latex_info.round_to_n(core_inferred_mass, 2)*u.M_sun,
           schmiedeke_summary_table[sst_mask]['M∗ initial']*u.M_sun,
           schmiedeke_summary_table[sst_mask]['M∗ all']*1e3*u.M_sun,
           latex_info.round_to_n(inferred_mass/sgrb2_age_myr/1e6,2)*u.M_sun/u.yr,
          ]

    if name == 'Total':
        totalrow = row
        totalrow[5] = latex_info.round_to_n(total_picking_max.value, 2)*u.M_sun
    else:
        tbl.add_row(row)


# right now, the total_picking_max is only the clustered
clustered_picking_max = copy.copy(total_picking_max)
# but, we need to subtract off NE and S...
NEandS = allclusters & ~NandM


# now get the non-clustered total
mask = ~allclusters

nhii = (hii & mask).sum()
ncores = ((~hii) & mask).sum()

nhii_clustered = (hii & ~mask).sum()
ncores_clustered = ((~hii) & ~mask).sum()

nhii_NM = (hii & NandM).sum()
ncores_NM = ((~hii) & NandM).sum()

# not meaningful mass = ncores * over8lt20mean + nhii * over20mean
# the fractions are fractions-of-total, so we're estimating the total mass
# from each independent population assuming they're the same age
hii_only_inferred_mass = nhii * over20mean / over20fraction
core_inferred_mass = ncores * over8lt20mean / over8lt20fraction
inferred_mass = (core_inferred_mass +
                 hii_only_inferred_mass) / 2.
print("Clustered_picking_max: {0} = {1}".format(clustered_picking_max, total_picking_max))
total_picking_max += max([hii_only_inferred_mass, core_inferred_mass])*u.M_sun
print("Total picking max: {0}".format(clustered_picking_max))
totalrow[5] = latex_info.round_to_n(total_picking_max.value, 2)*u.M_sun

clustered_hii_only_inferred_mass = nhii_clustered * over20mean / over20fraction
clustered_core_inferred_mass = ncores_clustered * over8lt20mean / over8lt20fraction

NMclustered_hii_only_inferred_mass = nhii_NM * over20mean / over20fraction
NMclustered_core_inferred_mass = ncores_NM * over8lt20mean / over8lt20fraction
NMclustered_total = (hii & maskM).sum() * over20mean/over20fraction + ((~hii) & maskN).sum() * over8lt20mean/over8lt20fraction


tbl.add_row(['Unassociated',
             ncores,
             nhii,
             latex_info.round_to_n(mass,2)*u.M_sun,
             latex_info.round_to_n(inferred_mass, 2)*u.M_sun,
             latex_info.round_to_n(core_inferred_mass, 2)*u.M_sun,
             latex_info.round_to_n(hii_only_inferred_mass, 2)*u.M_sun,
             latex_info.round_to_n(core_inferred_mass, 2)*u.M_sun,
             -999,
             -999,
             latex_info.round_to_n(inferred_mass/sgrb2_age_myr/1e6,2)*u.M_sun/u.yr,
            ])
tbl.add_row(totalrow)
#tbl.add_row(['Total$_{max}$', -999, -999, -999,
#             latex_info.round_to_n(total_picking_max.value,2),
#             latex_info.round_to_n(total_picking_max.value,2), -999, -999,
#             -999, -999,
#             latex_info.round_to_n(total_picking_max.value / sgrb2_age_myr / 1e6,2)])


print("SFR({1} Myr) = {0}".format(totalrow[-1], sgrb2_age_myr))
sfrbrickage = (latex_info.round_to_n(totalrow[4].to(u.Msun).value/(sgrb2_brick_age_myr*1e6), 2)*(u.Msun/u.yr))
print("SFR({1} Myr) = {0}".format(sfrbrickage, sgrb2_brick_age_myr))
sfrdynage_totmax = latex_info.round_to_n(total_picking_max.value / sgrb2_age_myr / 1e6,2)
print("SFR({1} Myr) = {0}".format(sfrdynage_totmax, sgrb2_age_myr))
sfrbrickage_totmax = (latex_info.round_to_n(total_picking_max.to(u.Msun).value/(sgrb2_brick_age_myr*1e6), 2)*(u.Msun/u.yr))
print("SFR({1} Myr) = {0}".format(sfrbrickage_totmax, sgrb2_brick_age_myr))

with open(paths.texpath('sfr.tex'), 'w') as fh:
    fh.write("\\newcommand{{\\sfrdynage}}{{{0}\\xspace}}\n".format(totalrow[-1].to(u.Msun/u.yr).value))
    fh.write("\\newcommand{{\\sfrbrickage}}{{{0}\\xspace}}\n".format(sfrbrickage.to(u.Msun/u.yr).value))
    fh.write("\\newcommand{{\\sfrdynagemax}}{{{0}\\xspace}}\n".format(sfrdynage_totmax))
    fh.write("\\newcommand{{\\sfrbrickagemax}}{{{0}\\xspace}}\n".format(sfrbrickage_totmax.to(u.Msun/u.yr).value))


formats = {'SFR': lambda x: latex_info.strip_trailing_zeros(str(x)),
           '$M_{\\rm inf}^s$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm count}^s$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$N({\\rm cores})$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$N({\\rm H\\textsc{ii}})$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm count}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred, H\\textsc{ii}}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred, cores}$': lambda x: "{0}".format(x) if x != -999 else '-',
          }

latexdict = latex_info.latexdict.copy()
latexdict['header_start'] = '\label{tab:clustermassestimates}'
latexdict['caption'] = 'Cluster Masses'
latexdict['preamble'] = '\centering\n\\begin{minipage}{130mm}'
latexdict['tablefoot'] = ("\par\n"
                          "$M_{{\\rm count}}$ is the mass of directly counted protostars, "
                          "assuming each millimeter source is {0:0.1f} \msun, or "
                          "{1:0.1f} \msun "
                          "if it is also an \hii region.  "
                          "$M_{{\\rm inferred,cores}}$ and $M_{{\\rm inferred,\hii}}$ are the inferred "
                          "total stellar masses assuming the counted objects represent "
                          "fractions of the total mass {2:0.2f} (cores) and "
                          "{3:0.2f} (\hii regions).  $M_{{\\rm inferred}}$ is the average "
                          "of these two.  "
                          "$M_{{\\rm count}}^s$ and $M_{{\\rm inf}}^s$ are the counted and inferred "
                          "masses reported in \citet{{Schmiedeke2016a}}.  "
                          "The star formation rate is computed using $M_{{\\rm inferred}}$ and"
                          " an age $t=0.74$ Myr, "
                          "which is the time of the last pericenter passage in the "
                          "\citet{{Kruijssen2015a}} model."
                          "  The \emph{{Total}} column represents the total over the whole observed "
                          "region.  "
                          "The \emph{{Total}}$_{{\\rm max}}$ column takes the higher of $M_{{\\rm inferred,\hii}}$ "
                          "and $M_{{\\rm inferred,cores}}$ from each row and sums them.  "
                          "We have included \hii regions in the $N(\hii)$ counts  that are "
                          "\emph{{not}} included in our source table \\ref{{tab:photometry}} because "
                          "they are too diffuse, or because they are unresolved in "
                          "our data but were resolved in the \citet{{De-Pree2014a}} VLA data.  As a result, "
                          "the total source count is greater than the source count reported in Table \\ref{{tab:photometry}}. "
                          "Also, the unassociated \hii region count is incomplete; it is missing both diffuse \hii "
                          "regions and possibly unresolved hypercompact \hii regions, since there are no "
                          "VLA observations comparable to \citet{{De-Pree2014a}} in the unassociated regions."
                          .format(over8lt20mean, over20mean, over8lt20fraction,
                                  over20fraction)
                         )
tbl.write(paths.texpath('cluster_mass_estimates.tex'), format='ascii.latex',
          formats=formats,
          latexdict=latexdict, overwrite=True)



tbl.add_row(['Clustered with NE, S', ncores_clustered, nhii_clustered, -999,
             latex_info.round_to_n(clustered_picking_max.value,2),
             latex_info.round_to_n(clustered_picking_max.value,2),
             latex_info.round_to_n(clustered_hii_only_inferred_mass, 2),
             latex_info.round_to_n(clustered_core_inferred_mass, 2),
             -999, -999,
             latex_info.round_to_n(clustered_picking_max.value / sgrb2_age_myr / 1e6,2)])
tbl.add_row(['Clustered only M, N', ncores_NM, nhii_NM, -999,
             latex_info.round_to_n(NMclustered_total,2),
             latex_info.round_to_n(NMclustered_total,2),
             latex_info.round_to_n(NMclustered_hii_only_inferred_mass, 2),
             latex_info.round_to_n(NMclustered_core_inferred_mass, 2),
             -999, -999,
             latex_info.round_to_n(clustered_picking_max.value / sgrb2_age_myr / 1e6,2)])

latexdict['tablefoot'] = ("\\\\\n"
                          "Partial reproduction of Table 2 in \citet{{Ginsburg2018a}}. "
                          #"$M_{{count}}$ is the mass of directly counted protostars, "
                          #"assuming each millimeter source is {0:0.1f} \msun, or "
                          #"{1:0.1f} \msun "
                          #"if it is also an \hii region.  "
                          "$M_{{\\rm inferred,cores}}$ and $M_{{\\rm inferred,\hii}}$ are the inferred "
                          "total stellar masses assuming the counted objects represent "
                          "fractions of the total mass of {2:0.2f} (cores) and "
                          "{3:0.2f} (\hii regions).  "
                          "$M_{{\\rm inferred,max}}$ is the greater "
                          "of these two.  "
                          #"$M_{{inferred}}$ is the average "
                          #"of these two.  "
                          #"$M_{{count}}^s$ and $M_{{inf}}^s$ are the counted and inferred "
                          #"masses reported in \citet{{Schmiedeke2016a}}.  "
                          #"The star formation rate is computed using $M_{{inferred}}$ and"
                          #" an age $t=0.74$ Myr, "
                          #"which is the time of the last pericenter passage in the "
                          #"\citet{{Kruijssen2015a}} model."
                          "The \emph{{Total}} row represents the total over the whole observed "
                          "region.  "
                          "The two \emph{{Clustered}} rows show the total inferred mass of clusters "
                          "including all "
                          "four candidate clusters M, N, NE and S, then the mass of clusters "
                          "including only M and N."
                          "\n\end{{minipage}}"
                          #"The \emph{{Total}}$_{{max}}$ column takes the higher of $M_{{inferred,\hii}}$ "
                          #"and $M_{{inferred,cores}}$ from each row and sums them.  "
                          #"We have included \hii regions in the $N(\hii)$ counts  that are "
                          #"\emph{{not}} included in our source table \\ref{{tab:photometry}} because "
                          #"they are too diffuse, or because they are unresolved in "
                          #"our data but were resolved in the \citet{{De-Pree2014a}} VLA data.  As a result, "
                          #"the total source count is greater than the source count reported in Table \\ref{{tab:photometry}}. "
                          #"Also, the unassociated \hii region count is incomplete; it is missing both diffuse \hii "
                          #"regions and possibly unresolved hypercompact \hii regions, since there are no "
                          #"VLA observations comparable to \citet{{De-Pree2014a}} in the unassociated regions."
                          .format(over8lt20mean, over20mean, over8lt20fraction,
                                  over20fraction)
                         )

formats = {'SFR': lambda x: latex_info.strip_trailing_zeros(str(x)),
           '$M_{\\rm inf}^s$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm count}^s$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$N({\\rm cores})$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$N({\\rm H\\textsc{ii}})$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm count}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred,max}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred, H\\textsc{ii}}$': lambda x: "{0}".format(x) if x != -999 else '-',
           '$M_{\\rm inferred, cores}$': lambda x: "{0}".format(x) if x != -999 else '-',
          }
cols = ['Name', '$N({\\rm cores})$', '$N({\\rm H\\textsc{ii}})$', '$M_{\\rm inferred, cores}$', '$M_{\\rm inferred, H\\textsc{ii}}$', '$M_{\\rm inferred,max}$', ]

tbl[cols].write(paths.cfepath('cluster_mass_estimates_cfe.tex'),
                format='ascii.latex', formats=formats, latexdict=latexdict,
                overwrite=True)
tbl.write(paths.cfepath('cluster_mass_estimates_cfe.csv'), format='ascii.csv',
          formats=formats,
          overwrite=True)
tbl.write(paths.tpath('cluster_mass_estimates_cfe.csv'), format='ascii.csv',
          formats=formats,
          overwrite=True)



core_phot_tbl.add_column(Column(name='Cluster', data=cluster_column))

classification = Column(name='Classification',
                        data=["{0}{1}{2}{3} {4}"
                              .format(('S' if row['color'] == 'green' else 'W'),
                                      ("\_" if row['Muno_xray_ID'] == '-' else "X"),
                                      ("\_" if row['Caswell_Name'] == '-' else "M"),
                                      ("\_" if np.isnan(row['McGrath_V_H2O']) else "W"),
                                      # would be nice... ("\_" if row['HII_name'] == '-' else "H"),
                                      str(row['SIMBAD_OTYPE']))
                              for row in core_phot_tbl])
core_phot_tbl.add_column(classification)

core_phot_tbl.write(paths.tpath("continuum_photometry_withSIMBAD_andclusters.ipac"),
                    format='ascii.ipac', overwrite=True)

# cutoff analysis...

for cutoff in np.linspace(8,100,10):
    x = np.linspace(cutoff,mmax,50000)
    y = kroupa(x)
    over20mean = (x*y).sum()/y.sum()
    over20fraction = (kroupa.m_integrate(cutoff, mmax)[0] /
                      kroupa.m_integrate(kroupa.mmin, mmax)[0])
    print("Mass fraction M>{1} = {0}".format(over20fraction, cutoff))
    print("Mean mass M>{1} = {0}".format(over20mean, cutoff))
    print("Represented mass M>{1} = {0}".format(over20mean/over20fraction, cutoff))


print()
print(">x, <20")

for low_cutoff in [8,10,12,14,15]:
    x = np.linspace(low_cutoff,hii_cutoff,50000)
    y = kroupa(x)
    overXlt20mean = (x*y).sum()/y.sum()
    overXlt20fraction = (kroupa.m_integrate(low_cutoff, hii_cutoff)[0] /
                         kroupa.m_integrate(kroupa.mmin, mmax)[0])
    overXlt20representative = overXlt20mean/overXlt20fraction
    print("Mass fraction 20>M>{1} = {0}".format(overXlt20fraction, low_cutoff))
    print("Mean mass 20>M>{1} = {0}".format(overXlt20mean, low_cutoff))
    print("Represented mass 20>M>{1} = {0}".format(overXlt20mean/overXlt20fraction, low_cutoff))


# cloud mass estimate
datapath = '/Users/adam/work/sgrb2/alma/FITS/continuumdata'
colfilename = datapath+'/column_maps/scuba_col_herscheltem.fits'

colfile = fits.open(colfilename)[0]

sgrb2contfile = fits.open(paths.Fpath('merge/continuum/SgrB2_selfcal_full_TCTE7m_selfcal5_ampphase_taylorterms_multiscale_deeper_mask2.5mJy.image.tt0.pbcor.fits'))
mywcs = wcs.WCS(colfile.header)
pix_area = wcs.utils.proj_plane_pixel_area(mywcs)*u.deg**2

observed_region,_ = reproject.reproject_interp(sgrb2contfile, colfile.header)

colmap = (u.Quantity(colfile.data, u.cm**-2)-5e22*u.cm**-2)

total_mass = (colmap[np.isfinite(observed_region)] * (pix_area * distance**2).to(u.cm**2, u.dimensionless_angles()) * 2.8*u.Da).to(1e6*u.M_sun).sum()
print("Total cloud mass: {0}".format(total_mass))
total_mass_fgsub = ((colmap[np.isfinite(observed_region)]-5e22*u.cm**-2) * (pix_area * distance**2).to(u.cm**2, u.dimensionless_angles()) * 2.8*u.Da).to(1e6*u.M_sun).sum()
print("Total cloud mass, 5e22 foreground subtracted: {0}".format(total_mass_fgsub))

sorted_colmap = np.sort(colmap[np.isfinite(observed_region)])
sorted_massmap = (sorted_colmap * (pix_area * distance**2).to(u.cm**2, u.dimensionless_angles()) * 2.8*u.Da).to(u.M_sun)
massmap_cumsum = np.cumsum(sorted_massmap*(sorted_massmap>0))
sorted_massmap_fgsub = ((sorted_colmap-5e22*u.cm**-2) * (pix_area * distance**2).to(u.cm**2, u.dimensionless_angles()) * 2.8*u.Da).to(u.M_sun)
massmap_cumsum_fgsub = np.cumsum(sorted_massmap_fgsub*(sorted_massmap_fgsub>0))

import pylab as pl
pl.figure(1).clf()
ax = pl.gca()
ax.loglog(sorted_colmap, massmap_cumsum[-1]-massmap_cumsum, label='No Foreground Subtraction')
ax.loglog(sorted_colmap-5e22*u.cm**-2, massmap_cumsum_fgsub[-1]-massmap_cumsum_fgsub, label='$N=5\\times10^{22}$ cm$^{-2}$ Foreground Subtraction')
pl.axis([7e21,1e25,5e4,2e6])
pl.grid()
ax.set_yticks([1e5, 5e5, 1e6])
ax.set_yticklabels(['$10^5$',r'$5\times10^5$','$10^6$'])
pl.legend(loc='best')
pl.xlabel("Column Density [N($H_2$) cm$^{-2}$]")
pl.ylabel("Cumulative Mass at $N>N_{x}$ [M$_{\odot}$]")
pl.savefig(paths.fpath('cumulative_mass_plot.pdf'), bbox_inches='tight')

reff = ((pix_area * np.arange(sorted_colmap.size-1, -1, -1) / np.pi * distance**2)**(1/2.)).to(u.pc, u.dimensionless_angles())
neff = (massmap_cumsum / (4/3*np.pi*reff**3) / (2.8*u.Da)).to(u.cm**-3)
rhoeff = (massmap_cumsum / (4/3*np.pi*reff**3)).to(u.g*u.cm**-3)
tff = (3 * np.pi / (32 * constants.G * rhoeff))**0.5

tbl2 = Table.read(paths.tpath("continuum_photometry_plusbackground.ipac"), format='ascii.ipac',)

# this is to try to make an argument about whether M_*(N>Nx) / tff(N>Nx) is getting higher or lower...
pl.figure(2).clf()
ax = pl.gca()
ax.loglog(sorted_colmap, tff.to(u.yr))
ax.hlines([0.74e6, 0.43e6,], 1e22, 1e25)
ax.vlines(tbl2['ScubaHTemColumn'][(np.array([0.9,0.5])*len(tbl2)).astype('int')], 1e3, 1e6,)


"""
Result as of 3/24/2017:
Cluster M   : N(cores)= 10 N(HII)= 37 counted mass=   1803.50 inferred mass=   6719.27 HII-only inferred mass:   12084.30 core-inferred mass=   1354.25
Cluster N   : N(cores)= 10 N(HII)=  4 counted mass=    301.72 inferred mass=   1330.33 HII-only inferred mass:    1306.41 core-inferred mass=   1354.25
Cluster NE  : N(cores)=  4 N(HII)=  0 counted mass=     47.87 inferred mass=    270.85 HII-only inferred mass:       0.00 core-inferred mass=    541.70
Cluster S   : N(cores)=  3 N(HII)=  1 counted mass=     81.41 inferred mass=    366.44 HII-only inferred mass:     326.60 core-inferred mass=    406.27

Updated 4/24/2017:
Cluster M   : N(cores)= 10 N(HII)= 55 counted mass=   2622.65 inferred mass=   9658.70 HII-only inferred mass:   17963.15 core-inferred mass=   1354.25
Cluster N   : N(cores)= 10 N(HII)=  7 counted mass=    438.24 inferred mass=   1820.23 HII-only inferred mass:    2286.22 core-inferred mass=   1354.25
Cluster NE  : N(cores)=  4 N(HII)=  0 counted mass=     47.87 inferred mass=    270.85 HII-only inferred mass:       0.00 core-inferred mass=    541.70
Cluster S   : N(cores)=  3 N(HII)=  1 counted mass=     81.41 inferred mass=    366.44 HII-only inferred mass:     326.60 core-inferred mass=    406.27


still missing 2 in NE, 5, in M, 1 in S, and have one too many in N
"""
