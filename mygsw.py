import os
import numpy as np

"""""""""""""
These functions are all from the deprecated all python version
of gibbs sea water. They are included here for the sole purpose of being 
able to recreate the geostrophic streamfunction along isopycnals function

"""""""""""""

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

class Bunch(dict):
    """
    A dictionary that also provides access via attributes.

    Additional methods update_values and update_None provide
    control over whether new keys are added to the dictionary
    when updating, and whether an attempt to add a new key is
    ignored or raises a KeyError.

    The Bunch also prints differently than a normal
    dictionary, using str() instead of repr() for its
    keys and values, and in key-sorted order.  The printing
    format can be customized by subclassing with a different
    str_ftm class attribute.  Do not assign directly to this
    class attribute, because that would substitute an instance
    attribute which would then become part of the Bunch, and
    would be reported as such by the keys() method.

    To output a string representation with
    a particular format, without subclassing, use the
    formatted() method.
    """

    str_fmt = "{0!s:<{klen}} : {1!s:>{vlen}}\n"

    def __init__(self, *args, **kwargs):
        """
        *args* can be dictionaries, bunches, or sequences of
        key,value tuples.  *kwargs* can be used to initialize
        or add key, value pairs.
        """
        dict.__init__(self)
        self.__dict__ = self
        for arg in args:
            self.update(arg)
        self.update(kwargs)

    def __str__(self):
        return self.formatted()

    def formatted(self, fmt=None, types=False):
        """
        Return a string with keys and/or values or types.

        *fmt* is a format string as used in the str.format() method.

        The str.format() method is called with key, value as positional
        arguments, and klen, vlen as kwargs.  The latter are the maxima
        of the string lengths for the keys and values, respectively,
        up to respective maxima of 20 and 40.
        """
        if fmt is None:
            fmt = self.str_fmt

        items = list(self.items())
        items.sort()

        klens = []
        vlens = []
        for i, (k, v) in enumerate(items):
            lenk = len(str(k))
            if types:
                v = type(v).__name__
            lenv = len(str(v))
            items[i] = (k, v)
            klens.append(lenk)
            vlens.append(lenv)

        klen = min(20, max(klens))
        vlen = min(40, max(vlens))
        slist = [fmt.format(k, v, klen=klen, vlen=vlen) for k, v in items]
        return ''.join(slist)

    def from_pyfile(self, filename):
        """
        Read in variables from a python code file.
        """
        # We can't simply exec the code directly, because in
        # Python 3 the scoping for list comprehensions would
        # lead to a NameError.  Wrapping the code in a function
        # fixes this.
        d = dict()
        lines = ["def _temp_func():\n"]
        with open(filename) as f:
            lines.extend(["    " + line for line in f])
        lines.extend(["\n    return(locals())\n",
                      "_temp_out = _temp_func()\n",
                      "del(_temp_func)\n"])
        codetext = "".join(lines)
        code = compile(codetext, filename, 'exec')
        exec(code, globals(), d)
        self.update(d["_temp_out"])
        return self

    def update_values(self, *args, **kw):
        """
        arguments are dictionary-like; if present, they act as
        additional sources of kwargs, with the actual kwargs
        taking precedence.

        One reserved optional kwarg is "strict".  If present and
        True, then any attempt to update with keys that are not
        already in the Bunch instance will raise a KeyError.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
            newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in newkw.items() if k in self])
        self.update(dsub)

    def update_None(self, *args, **kw):
        """
        Similar to update_values, except that an existing value
        will be updated only if it is None.
        """
        strict = kw.pop("strict", False)
        newkw = dict()
        for d in args:
            newkw.update(d)
        newkw.update(kw)
        self._check_strict(strict, newkw)
        dsub = dict([(k, v) for (k, v) in newkw.items()
                                if k in self and self[k] is None])
        self.update(dsub)

    def _check_strict(self, strict, kw):
        if strict:
            bad = set(kw.keys()) - set(self.keys())
            if bad:
                bk = list(bad)
                bk.sort()
                ek = list(self.keys())
                ek.sort()
                raise KeyError(
                    "Update keys %s don't match existing keys %s" % (bk, ek))



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
