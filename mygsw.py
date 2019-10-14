import os
import numpy as np
from gsw._utilities import Bunch
import gsw


"""
These functions are all from the deprecated all python version
of gibbs sea water. They are included here for the sole purpose of being 
able to recreate the geostrophic streamfunction along isopycnals function

"""

def interp_S_T(S, T, z, znew, P=None):
    """
    Linear interpolation of ndarrays *S* and *T* from *z* to *znew*.
    Optionally interpolate a third ndarray, *P*.
    *z* must be strictly increasing or strictly decreasing.  It must
    be a 1-D array, and its length must match the last dimension
    of *S* and *T*.
    *znew* may be a scalar or a sequence.
    It is assumed, but not checked, that *S*, *T*, and *z* are
    all plain ndarrays, not masked arrays or other sequences.
    Out-of-range values of *znew*, and *nan* in *S* and *T*,
    yield corresponding *nan* in the output.
    The basic algorithm is from scipy.interpolate.
    """

    isscalar = False
    if not np.iterable(znew):
        isscalar = True
        znew = [znew]
    znew = np.asarray(znew)
    inverted = False
    if z[1] - z[0] < 0:
        inverted = True
        z = z[::-1]
        S = S[..., ::-1]
        T = T[..., ::-1]
        if P is not None:
            P = P[..., ::-1]
    if (np.diff(z) <= 0).any():
        raise ValueError("z must be strictly increasing or decreasing")
    hi = np.searchsorted(z, znew)
    hi = hi.clip(1, len(z) - 1).astype(int)
    lo = hi - 1
    z_lo = z[lo]
    z_hi = z[hi]
    S_lo = S[lo]
    S_hi = S[hi]
    T_lo = T[lo]
    T_hi = T[hi]
    zratio = (znew - z_lo) / (z_hi - z_lo)
    Si = S_lo + (S_hi - S_lo) * zratio
    Ti = T_lo + (T_hi - T_lo) * zratio
    if P is not None:
        Pi = P[lo] + (P[hi] - P[lo]) * zratio
    if inverted:
        Si = Si[..., ::-1]
        Ti = Ti[..., ::-1]
        if P is not None:
            Pi = Pi[..., ::-1]
    outside = (znew < z.min()) | (znew > z.max())
    if np.any(outside):
        Si[..., outside] = np.nan
        Ti[..., outside] = np.nan
        if P is not None:
            Pi[..., outside] = np.nan
    if isscalar:
        Si = Si[0]
        Ti = Ti[0]
        if P is not None:
            Pi = Pi[0]
    if P is None:
        return Si, Ti
    return Si, Ti, Pi

class Cache_npz(object):
    def __init__(self):
	    self._cache = dict()
	    self._default_path = os.path.join(os.path.dirname(__file__), 'data')

    def __call__(self, fname, datadir=None):
	    if datadir is None:
		    datadir = self._default_path
	    fpath = os.path.join(datadir, fname)
	    try:
		    return self._cache[fpath]
	    except KeyError:
		    pass
	    d = np.load(fpath)
	    ret = Bunch(d)
	    self._cache[fpath] = ret
	    return ret

_npz_cache = Cache_npz()

def enthalpy_SSO_0(p):
    """
    Enthalpy at SSO and CT(T=0) (75-term equation)
    This function calculates enthalpy at the Standard Ocean Salinty, SSO,
    and at a Conservative Temperature of zero degrees C, as a function of
    pressure, p, in dbar, using a streamlined version of the 75-term
    computationally-efficient expression for specific volume, that is, a
    streamlined version of the code "enthalpy(SA,CT,p)".
    Parameters
    ----------
    p : array_like
        pressure [dbar]
    Returns
    -------
    enthalpy_SSO_0 : array_like
                     enthalpy at (SSO, CT=0, p)
                     [J kg :sup:`-1`]
    """

    h006 = -2.1078768810e-9
    h007 =  2.8019291329e-10

    z = p * 1e-4
    db2Pascal = 1e4

    dynamic_enthalpy_SSO_0_p = z * (9.726613854843870e-4 + z * (
        -2.252956605630465e-5 + z * (2.376909655387404e-6 + z * (
            -1.664294869986011e-7 + z * (-5.988108894465758e-9 + z * (
                h006 + h007 * z))))))

    return dynamic_enthalpy_SSO_0_p * db2Pascal * 1e4

def read_data(fname, datadir=None):
    """
    Read variables from a numpy '.npz' file into a minimal class providing
    attribute access.  A cache is used to avoid re-reading the same file.
    """
    return _npz_cache(fname, datadir=datadir)


def interp_ref_cast(spycnl, A="gn"):
    """
    Translation of:
    [SA_iref_cast, CT_iref_cast, p_iref_cast] = interp_ref_cast(spycnl, A)
    interp_ref_cast            linear interpolation of the reference cast
    ==========================================================================
    This function interpolates the reference cast with respect to the
    interpolating variable "spycnl".  This reference cast is at the location
    188E,4N from the reference data set which underlies the Jackett &
    McDougall (1997) Neutral Density computer code.  This function finds the
    values of SA, CT and p on this reference cast which correspond to the
    value of isopycnal which is passed to this function from the function
    "geo_strf_isopycnal_CT".  The isopycnal could be either gamma_n or
    sigma_2. If A is set to any of the following 's2','S2','sigma2','sigma_2'
    the interpolation will take place in sigma 2 space, any other input
    will result in the programme working in gamma_n space.
    VERSION NUMBER: 3.0 (14th April, 2011)
    REFERENCE:
    Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable
    for the world<92>s oceans. Journal of Physical Oceanography, 27, 237-263.
    FIXME? Do we need argument checking here to handle masked arrays,
    etc.?  I suspect not, since I don't think this is intended to be
    user-callable, but is instead used internally by user-callable
    functions.
    Note: The v3.03 matlab code is incorrectly using approximate numbers
    for the gamma_n case, even when the sigma_2 case is in effect.
    That bug is fixed here.
    """

    if A.lower() in ["s2", "sigma2", "sigma_2"]:
	    A = "s2"
    gsw_data = read_data("gsw_data_v3_0.npz")
    SA_ref = gsw_data.SA_ref_cast
    CT_ref = gsw_data.CT_ref_cast
    p_ref = gsw_data.p_ref_cast
    if A == "s2":
	    zvar_ref = gsw_data.sigma_2_ref_cast
    else:
	    zvar_ref = gsw_data.gamma_n_ref_cast
    zvar_new = spycnl
    Si, Ci, Pi = interp_S_T(SA_ref, CT_ref, zvar_ref, zvar_new, P=p_ref)
    shallower = spycnl <= zvar_ref[0]
    deeper = spycnl >= zvar_ref[-1]
    if shallower.any():
	    Si[shallower] = SA_ref[0]
	    Ci[shallower] = CT_ref[0]
	    Pi[shallower] = p_ref[0]
    if deeper.any():
	    Si[deeper] = SA_ref[-1]
	    Ci[deeper] = CT_ref[-1]
	    Pi[deeper] = p_ref[-1]
    return Si, Ci, Pi

def geo_strf_isopycnal(SA,CT,p,p_ref,Neutral_Density,p_Neutral_Density):
    p = np.abs(np.asarray(p))

    dyn_height = gsw.geo_strf_dyn_height(SA,CT,p,p_ref)
    
    p_Neutral_Density, idxs = np.unique(p_Neutral_Density,return_index=True)
    Neutral_Density = Neutral_Density[idxs]

    filt = np.where(np.isin(p,p_Neutral_Density))

    nsdyn_height = dyn_height[filt]
    nstemps = CT[filt]
    nssals = SA[filt]

    db2Pa = 1e4
    sa_iref_cast,ct_iref_cast,p_iref_cast = interp_ref_cast(Neutral_Density,"s2")
    cp0 = 3991.86795711963
    #print("##################")
    #print("py iref_cast: ",p_iref_cast,ct_iref_cast,sa_iref_cast)
    #print("py nssal and nstemps: ",nssals,nstemps,ns)
    part1 = 0.5 *db2Pa*(p_Neutral_Density-p_iref_cast)*(gsw.specvol(nssals,nstemps,p_Neutral_Density)-gsw.specvol(sa_iref_cast,ct_iref_cast,p_Neutral_Density))
    part2 = 0
    part3 = (-0.225e-15)*(db2Pa*db2Pa)*(nstemps-ct_iref_cast)*(p_Neutral_Density-p_iref_cast)*(p_Neutral_Density-p_iref_cast)
    part4 = nsdyn_height - enthalpy_SSO_0(p_Neutral_Density)
    part5 = gsw.enthalpy(sa_iref_cast,ct_iref_cast,p_Neutral_Density) -cp0*ct_iref_cast
    #print("specvol delta", (gsw.specvol(nssals,nstemps,ns)-gsw.specvol(sa_iref_cast,ct_iref_cast,ns)))
    #print("py part1: ",part1+part2)
    #print("py part2: ",part3)
    #print("dyn height", nsdyn_height)
    #print("SS0 enthalpy",  mygsw.enthalpy_SSO_0(ns))
    #print("py part3: ",part4+part5)
    #print("result:", part1+part2+part3+part4+part5)
    #print("##################")
    return part1 + part2 + part3 + part4 + part5
