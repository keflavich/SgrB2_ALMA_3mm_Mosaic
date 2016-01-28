import numpy as np
from casa import table as tb

def goodenough_field_solutions(tablename, minsnr=5, maxphasenoise=np.pi/4.,
                               pols=[0]):
    """
    After an initial self-calibration run, determine which fields have good
    enough solutions.  This only inspects the *phase* component of the
    solutions.

    Parameters
    ----------
    tablename : str
        The name of the calibration table (e.g., phase.cal)
    minsnr : float
        The minimum *average* signal to noise ratio for a given field
    maxphasenoise : float
        The maximum average phase noise permissible for a given field in
        radians
    pols : list
        The list of polarizations to include in the heuristics

    Returns
    -------
    An array of field IDs.  This will need to be converted to a list of strings
    for use in CASA tasks
    """
    tb.open(tablename)
    solns = tb.getcol('CPARAM')
    fields = tb.getcol('FIELD_ID')
    snr = tb.getcol('SNR')

    all_angles = []
    all_snrs = []

    ufields = np.unique(fields)

    for field in ufields:
        sel = fields==field
        angles = np.angle(solns[:,:,sel])
        all_angles.append(angles)
        all_snrs.append(snr[:,:,sel])

    all_angles = np.array(all_angles)
    all_snrs = np.array(all_snrs)
    # not sure what 2nd column of these is...

    # first is mean across pols, then mean across all data points
    good_enough = ((all_snrs[:,pols,0,:].mean(axis=1).mean(axis=1) > minsnr) &
                   (all_angles[:,pols,0,:].std(axis=2).mean(axis=1) < maxphasenoise))

    return ufields[good_enough]

