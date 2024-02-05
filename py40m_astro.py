""" Astronomy-related routines.

TODO:
- Update docs
- Update style
- Fix pyslalib-dependent code
"""

from warnings import warn

import ephem
import numpy as np
from scipy.integrate import quad

import py40m_calendar as pcal

#from pyslalib import slalib

OVRO_LONGITUDE='-118:16:53.83'
OVRO_LATITUDE='37:13:55.7'
OVRO_ELEVATION=1236 # meters

VLA_LONGITUDE='-107:37:03.819'
VLA_LATITUDE='34:04:43.497'
VLA_ELEVATION=2124 # meters


# Observers
# OVRO site coordinates from http://www.astro.caltech.edu/~tjp/ovroman/introduction.html
def observer_ovro40m(*args,**kwargs):
    ovro40m=ephem.Observer(*args,**kwargs)
    ovro40m.long=ephem.degrees(OVRO_LONGITUDE)
    ovro40m.lat=ephem.degrees(OVRO_LATITUDE)
    ovro40m.elevation=OVRO_ELEVATION
    return ovro40m

def observer_vla(*args,**kwargs):
    vla=ephem.Observer(*args,**kwargs)
    vla.long=ephem.degrees(VLA_LONGITUDE)
    vla.lat=ephem.degrees(VLA_LATITUDE)
    vla.elevation=VLA_ELEVATION
    return vla


def ovro40m_azza(body, time, pressure=None, temp=None):
    """Compute the AZ and ZA of a specified ephem.body.

    Returns data in degrees.

    Body can be string name of a planet or 'sun' or 'moon' or
    an ephem.body object.

    Computes apparent coordinates, including a correction for
    atmospheric refraction.  Pass in pressure (mbar) and temp (deg C) if
    you want to use other than the default (1010mbar, 25 deg C).  If
    pressure is 0, no refraction correction will be used.

    This set of routines seems to agree to ~ 1 arcsec with NASA/JPL's
    HORIZONS portal, at least for solar position.

    """
    (az,el) = ovro40m_azel(body,time,pressure,temp)
    return (az*180.0/np.pi, 90.0-el*180.0/np.pi)


# def ovro40m_azel2(ra, dec, utc, pressure=None, temp=None):
#     mjd=pcal.get_mjd(utc)
#     dtt = slalib.sla_dtt(mjd) # offset to convert UTC to terrestrial (dynamical) time
#     # get apparent ra/dec
#     ra,da = slalib.sla_map(rm=ra, dm=dec,
#                            pr=0, pd=0, # no proper motion
#                            px=0, # no parallax
#                            rv=0, # no radial velocity
#                            eq=2000.0,
#                            date=mjd+dtt)
#     # now get observed
#     ao,zo,ho,do,ro = slalib.sla_aop(rap=ra, dap=da,
#                                     date=mjd,
#                                     dut=0, # ignoring for now
#                                     elongm=float(ephem.degrees(OVRO_LONGITUDE)),
#                                     phim=float(ephem.degrees(OVRO_LATITUDE)),
#                                     hm=OVRO_ELEVATION,
#                                     xp=0, yp=0, # ignoring polar motion
#                                     tdk=293.15, # assume this
#                                     pmb=865.0, # assume this
#                                     rh=0.2, # assume this
#                                     wl=2.0e4, # 2cm in um
#                                     tlr=0)
#     return ao, np.pi/2-zo


def ovro40m_azel(body, time, pressure=None, temp=None):
    """Compute the az and elevation of a specified ephem.body.

    Returns data as ephem.angle objects.

    Body can be string name of a planet or 'sun' or 'moon' or
    an ephem.body object.

    Computes apparent coordinates, including a correction for
    atmospheric refraction.  Pass in pressure (mbar) and temp (deg C) if
    you want to use other than the default (1010mbar, 25 deg C).  If
    pressure is 0, no refraction correction will be used.

    This set of routines seems to agree to ~ 1 arcsec with NASA/JPL's
    HORIZONS portal, at least for solar position.

    The time input should be an ephem.date, or a datetime or a string
    that can be converted to one by py40m.calendar.get_datetime().

    """
    if isinstance(body, str):
        body=body.lower()
        if body=='mercury':
            body=ephem.Mercury()
        elif body=='venus':
            body=ephem.Venus()
        elif body=='earth':
            body=ephem.Earth()
        elif body=='mars':
            body=ephem.Mars()
        elif body=='jupiter':
            body=ephem.Jupiter()
        elif body=='saturn':
            body=ephem.Saturn()
        elif body=='uranus':
            body=ephem.Uranus()
        elif body=='neptune':
            body=ephem.Neptune()
        elif body=='pluto':
            body=ephem.Pluto()
        elif body=='sun':
            body=ephem.Sun()
        elif body=='moon':
            body=ephem.Moon()
        else:
            raise Exception('Unknown body name')
    else:
        try:
            ra=body[0]
            dec=body[1]
        except TypeError:
            pass
        else:
            body=body_from_radec(ra,dec)


    observer=observer_ovro40m()

    # convert time into an ephem.date
    observer.date=get_ephem_date(time)

    if pressure is not None:
        observer.pressure=pressure
    if temp is not None:
        observer.temp=temp
    body.compute(observer)
    return (body.az, body.alt)


def body_from_radec(ra, dec, epoch='2000'):
    """Create and compute a pyephem body from ra and dec, in degrees.

    RA and dec should be ephem.Angle objects.  Epoch may be specified if
    desired, otherwise defaults to '2000' (should be a string).

    """
    b = ephem.FixedBody()
    b._ra=ra
    b._dec=dec
    b._epoch=epoch
    b.compute(ephem.date('2011/04/20 00:00:00'), epoch=epoch)
    # Note that we use a completely arbitrary date for the compute.
    # This is to (hopefully) guarantee repeatability.  It seems that
    # internals within ephem.body vary with time, even for
    # fully-specified coordinates.  This particularly affects sky
    # separation computations of identical angles, which sometimes
    # seem to give infinite separations.

    return b


def radec_from_azel(az, el, lst, lat=ephem.degrees(OVRO_LATITUDE)):
    """Compute the RA and dec corresponding to the azimuth and elevation of an observer.

    Defaults to OVRO latitude.  Note: no refraction correction is applied.

    """
    dec_rad = np.arcsin(np.sin(lat)*np.sin(el) + np.cos(lat)*np.cos(el)*np.cos(az))
    ha_rad = np.arctan2( -np.sin(az)*np.cos(el),
                          (np.sin(el)-np.sin(lat)*np.sin(dec_rad))/np.cos(lat) )

    dec = ephem.degrees(dec_rad)
    ra_rad = lst-ha_rad
    if ra_rad < 0:
        ra_rad += 2*np.pi
    ra = ephem.hours(ra_rad)

    return (ra, dec)


# def radec_from_azel2(az, el, utc, long=float(ephem.degrees(OVRO_LONGITUDE)),
#                      lat=float(ephem.degrees(OVRO_LATITUDE)), elev=OVRO_ELEVATION):
#     mjd=pcal.get_mjd(utc)
#     # get apparent coordinates
#     rap,dap = slalib.sla_oap(type_bn='a',
#                              ob1=az, ob2=np.pi/2-el,
#                              date=mjd,
#                              dut=0, # ignore for now, should use MJD-based LUT
#                              elongm=long,
#                              phim=lat,
#                              hm=elev,
#                              xp=0, yp=0, # ignore polar motion, could also get from IERS table
#                              tdk=293.15,
#                              pmb=865, # assumed, typical value in mbar
#                              rh=0.2, # assumed, typical value in fraction
#                              wl=2.0e4, # 2cm in um
#                              tlr=0) # use example 0.0065 K/m instead?
#     # now get mean J2000 coordinates
#     dtt = slalib.sla_dtt(mjd) # offset to convert UTC to terrestrial (dynamical) time
#     rm,dm = slalib.sla_amp(ra=rap, da=dap,
#                            date=mjd+dtt,
#                            eq=2000)
#     if rm >= 2*np.pi:
#         rm -= 2*np.pi
#     return ephem.hours(rm), ephem.degrees(dm)


def sky_sep(ra1, dec1, ra2, dec2):
    """Compute the sky separation between two positions.

    ra and dec should be ephem.Angle objects.  Both should be referred
    to the same epoch.  Computations assume J2000, but this should not
    be relevant to sky computations.

    Return value is an ephem.degrees angle.

    """
    #thresh=0.05/3600 # threshold below which coordinates are identical, in degrees
    ## ephem.separation buggy for some identical coordinates
    #if (abs(ra1-ra2)*180/np.pi<thresh) and (abs(dec1-dec2)*180/np.pi<thresh):
    #    return ephem.degrees(0)

    #b1=body_from_radec(ra1, dec1)
    #b2=body_from_radec(ra2, dec2)
    return ephem.separation((ra1,dec1), (ra2,dec2))


def parallactic_angle(ra, dec, observer=None, date=None):
    """Compute the parallactic angle of a source.

    Pass an ephem.observer, or it will default to OVRO.  If date is
    not set, it will default to now() for the default observer or
    whatever is already set if one is passed in.

    """
    if observer is None:
        observer=observer_ovro40m()
    if date is not None:
        observer.date=get_ephem_date(date)

    lat = observer.lat
    t_sid = observer.sidereal_time()

    b = body_from_radec(ra, dec)

    return ephem.degrees(
        parallactic_angle_lst(b.ra, b.dec,
                              t_sid, lat))



# def parallactic_angle_slalib(ra, dec, lst, lat):
#     pa = slalib.sla_pa(lst-ra, dec, lat)
#     return pa


def parallactic_angle_lst(ra, dec, lst, lat):
    """Compute the parallactic angle given the LST and latitude.

    Sign convention is that parallactic angle > 0 when an object is
    west of the meridian.

    Everything in radians, including LST.  Returns radians.

    """
    ha=lst-ra
    return np.arctan2(np.cos(lat)*np.sin(ha),
                   np.cos(dec)*np.sin(lat)-np.sin(dec)*np.cos(lat)*np.cos(ha))


def equatorial_to_galactic(ra, dec, epoch="2000"):
    """Convert equatorial coordinates to galactic coordinates.

    RA and dec should be ephem.Angle objects (e.g., use ephem.hours()
    and ephem.degrees).  Epoch defaults to J2000 (specify only the
    year).

    Coordinates can also be specified as strings in any format
    understood by the ephem.Angle creation routines.

    Return value is long, lat, each an ephem.Angle object that prints
    as degrees.

    """
    equ=ephem.Equatorial(ra, dec, epoch=epoch)
    gal=ephem.Galactic(equ)
    return (gal.long, gal.lat)

def galactic_to_equatorial(long, lat, epoch="2000"):
    gal=ephem.Galactic(long, lat)
    equ=ephem.Equatorial(gal, epoch=epoch)
    return (equ.ra, equ.dec)


# helper function for get_epoch_name
def _get_angle_to_precision(x, precision, include_decimal=True):
    """Return a string corresponding to x at given precision.

    Precision is the number of decimal places after the point. If
    precision=0, no decimal place is included.

    Always includes 2 zero-padded digits before the decimal.

    Discards sign.

    """
    factor = 10**precision

    # The next lines seem to better preserve precision by scaling
    # before subtracting.  I suspect this does not fully guarantee
    # there will not be representational rounding errors, but it
    # provides some protection.
    int_part = int(x)
    frac_part = int(x*factor - int_part*factor)

    int_part = abs(int_part)
    frac_part = abs(frac_part)

    if int_part > 99:
        raise ValueError('Cannot represent %f with two integer digits.'%(x))

    int_str = "%02d"%(int_part)
    if precision > 0:
        if include_decimal:
            dec_str = '.'
        else:
            dec_str = ''
        frac_str = dec_str+"%0*d"%(precision,frac_part)
    else:
        frac_str=""
    return int_str+frac_str


def get_epoch_name(ra, dec,
                   epoch='2000', prefix=None,
                   ra_precision=None,
                   dec_precision=None,
                   input_epoch='2000'):
    """Compute the name for a source given its coordinates.

    Defaults to a J2000 name of the form 'JHHMM+DDMM'

    Epoch of the input ra and dec is specified by input_epoch.  The
    coordinates must be either ephem.angle objects or float radians.

    Specify the epoch as a year. Default is 2000.

    If a prefix is given, this will be appended before the coordinate
    portion of the name.  If nothing is specified, 'J' will be used
    for an epoch of 2000, no prefix for 1950.  To get no prefix, pass
    the empty string.

    The ra_precision and dec_precision fields give the precisions to
    use for the RA and Dec output portions. The defaults depend on the
    epoch. If precision is 0, then the format will be HHMM or DDMM.
    Precision > 0 will add decimal places, so if ra_precision=1 then
    the RA will be output as HHMM.m.  If it's < 0, then digits will be
    removed, so dec_precision=-1 will cause Dec to be output as DD.d
    (note the switch from DDMM to decimal degrees).  Precision < -1 is
    not implemented.  In the case of precision=-1, no decimal point is used.

    The default for J2000 is ra_precision=dec_precision=0. For B1950,
    the default is ra_precision=0, dec_precision=-1.

    """
    if ra < 0:
        raise ValueError('Negative RA not allowed (got %f).'%(ra))

    if ra_precision is not None and ra_precision < -1:
        raise ValueError('ra_precision < -1 not supported.')

    if dec_precision is not None and dec_precision < -1:
        raise ValueError('dec_precision < -1 not supported.')

    if epoch == '2000':
        if prefix is None:
            prefix='J'
        if ra_precision is None:
            ra_precision=0
        if dec_precision is None:
            dec_precision=0

    elif epoch == '1950':
        if prefix is None:
            prefix=''
        if ra_precision is None:
            ra_precision=0
        if dec_precision is None:
            dec_precision=-1

    coords = ephem.Equatorial(ra, dec, epoch=input_epoch)
    coords = ephem.Equatorial(coords, epoch=epoch)
    ra,dec = coords.get()

    ra_hms = (1.0/15.0) * (180.0/np.pi) * ra
    dec_dms = (180.0/np.pi) * dec

    if ra_precision == -1:
        ra_h_str = _get_angle_to_precision(ra_hms, 1, include_decimal=False)
        ra_ms_str=""
    else:
        ra_h = int(ra_hms)
        ra_ms = 60*ra_hms - 60*ra_h

        ra_h_str = "%02d"%(ra_h)
        ra_ms_str = _get_angle_to_precision(ra_ms, ra_precision)

    if dec_precision == -1:
        dec_d_str = _get_angle_to_precision(dec_dms, 1, include_decimal=False)
        dec_ms_str=""
    else:
        dec_d = int(dec_dms)
        dec_ms = 60*dec_dms - 60*dec_d
        dec_d = abs(dec_d)
        dec_ms = abs(dec_ms)

        dec_d_str = "%02d"%(dec_d)
        dec_ms_str = _get_angle_to_precision(dec_ms, dec_precision)

    if dec >= 0:
        s = ra_h_str+ra_ms_str + "+" + dec_d_str+dec_ms_str
    else:
        s = ra_h_str+ra_ms_str + "-" + dec_d_str+dec_ms_str

    return prefix+s


def log_nuFnu(nu, Fnu, err_Fnu):
    """Compute log nuFnu and its error.

    Computes nu*F_nu and its error, appropriate for plotting on log-log SED plots.

    Arguments
    nu - observation frequency, in Hz
    Fnu - observed flux in Jy
    err_Fnu - flux uncertainty in Jy


    Returns a tuple, R.  R[0]=log nuFnu in erg/s/cm^2.
    R[1]=the logarithmic error, in erg/s/cm^2.

    """
    log_nuFnu_Jy = np.log10(nu)+np.log10(Fnu)
    log_nuFnu_cgs = log_nuFnu_Jy - 23.0 # 1 Jy = 1e-23 erg/s/cm^2
    err_log_nuFnu_cgs = (err_Fnu/Fnu)/np.log(10.0)

    return (log_nuFnu_cgs, err_log_nuFnu_cgs)


#===========================================================
# Cosmology routines

# constants from WMAP 5yr
#     (Komatsu et al, 2009, ApJ, 180, 330)
c_m_per_s=299792458.0
c_km_per_s=c_m_per_s/1000.0
c_cm_per_s=c_m_per_s*100.0

WMAP5_LCDM_PARAMS={
    "H0":70.5, # km/s/Mpc
    "sigma_H0":1.3, # km/s/Mpc
    "Omega_b":0.0456,
    "sigma_Omega_b":0.0015,
    "Omega_c":0.228,
    "sigma_Omega_c":0.013,
    "Omega_M":0.0456+0.228, # Omega_b+Omega_c
    "sigma_Omega_M":np.sqrt(0.0015**2+0.013**2), # sigma(Omega_b+Omega_c)
    "Omega_Lambda":0.726,
    "sigma_Omega_Lambda":0.015
    }

LCDM_PARAMS=WMAP5_LCDM_PARAMS

def _integrate(f, x0, x1, max_err_ratio=1e-6):
    """Integrates using scipy.integrate.quad with warning on large error."""
    (result, err) = quad(f,x0,x1)
    if err/result >= max_err_ratio:
        warn("Warning: large integration error possible (ratio={0})".format(err/result))
    return result

#-----------------------------------------------------------
# comoving coord distance and helpers
#
def _comoving_coord_distance_value(z,Omega_M,Omega_L,H0):
    """Compute the CCD itself.

    Helper function.  Do not use outside the module.

    """
    return (1/H0) * _integrate(
        lambda zz: 1.0/(np.sqrt(Omega_L+Omega_M*((1+zz)**3))),
        0, z)

def _comoving_coord_distance_partials(z,Omega_M,Omega_L,H0,ccd):
    """Compute the partial derivs of CCD wrt its inputs.

    Helper function.  Do not use outside the module.

    Return value is (dCCD/dz, dCCD/dOmega_M, dCCD/dOmega_Lambda, dCCD/dH0)

    """
    dXdz=(1/H0) * 1.0/np.sqrt(Omega_L+Omega_M*(1+z)**3.0)
    dXdH0= -ccd / H0
    # here, we're not assuming Omega_L=1-Omega_M, account for this on your own.
    dXdOmega_M = -(1/(2*H0)) * _integrate(
        lambda zz: ((1+zz)**3.0) / ((1.0 + Omega_M*((1+zz)**3.0 - 1.0))**1.5),
        0, z)
    dXdOmega_L = -(1/(2*H0)) * _integrate(
        lambda zz: 1.0 / ((1.0 + Omega_M*((1+zz)**3.0 - 1.0))**1.5),
        0, z)
    return(dXdz,dXdH0,dXdOmega_M,dXdOmega_L)


def comoving_coord_distance_flatMD(z, H0=LCDM_PARAMS["H0"]/c_km_per_s):
    """Compute the comoving coordinate distance in Mpc in a flat MD universe."""
    return (1-1.0/np.sqrt(1+z))*2/H0

def comoving_coord_distance_LCDM(z,
                                 dz=0,
                                 Omega_M=LCDM_PARAMS["Omega_M"],
                                 dOmega_M=LCDM_PARAMS["sigma_Omega_M"],
                                 H0=LCDM_PARAMS["H0"]/c_km_per_s,
                                 dH0=LCDM_PARAMS["sigma_H0"]/c_km_per_s,
                                 q0=None):
    """Compute the comoving coordinate distance in Mpc for a flat, Lambda-CDM universe.

    Arguments:

    z - redshift for which to calculate comoving coordinate distance
    dz - uncertainty in z

    Omega_M - matter contribution to the universe as a fraction of
    critical density.  Omega_Lambda will be taken as 1-Omega_M to
    ensure flatness.
    dOmega_M - uncertainty in Omega_M

    H0 - Hubble rate today, in 1/[length unit].  This assumes c=1.
    Pass in H0/c if your H0 is in velocity/[length unit] where c is
    in your preferred units.
    dH0 - uncertainty in H0

    q0 - deceleration parameter. If this is provided, overrides Omega_M

    Return value: (dist, dist_err)

    dist is the comoving coordinate distance in [length unit].  This
    is defined as

               c   /\z                dz'
    Chi(z) = ----  |  ---------------------------------------
              H0  \/0  sqrt(Omega_Lambda + Omega_M*(1+z')^3)


    dist_err is the 1-sigma uncertainty in the comoving coordinate
    distance.  This is computed assuming the universe is flat, i.e.,
    Omega_M+Omega_Lambda=1 exactly, i.e., dOmega_Lambda=-dOmega_M.

    For integration, this routine uses scipy.integrate.quad.
    See the scipy documentation for the method.

    Parameters default to the WMAP 5-year values:
    (Komatsu et al, 2009, ApJ, 180, 330)

      H0 = 70.5 +/- 1.3 km / s / Mpc
         => H0/c = (2.352 +/- 0.043) x 10^-4 / Mpc
      Omega_M = Omega_CDM+Omega_Baryon = 0.274 +/- 0.013
      Omega_L = 1-Omega_M

    """
    if q0 is not None:
        Omega_M=(q0+1)*2./3.

    Omega_L = 1.0-Omega_M
    ccd=_comoving_coord_distance_value(z,Omega_M,Omega_L,H0)

    # partial derivs of comoving-coord-dist X
    (dXdz,
     dXdH0,
     dXdOmega_M,
     dXdOmega_L)=_comoving_coord_distance_partials(z,Omega_M,Omega_L,H0,ccd)

    dccd = np.sqrt(
        (dXdH0 * dH0) **2.0
        + (dXdz * dz) **2.0
        + ((dXdOmega_M-dXdOmega_L) * dOmega_M) ** 2.0)
    # Mix dXdOmega_M and dXdOmega_L on last line to constrain to flat universe.
    # The alternative would be to allow the universe to curve and to include
    # the curvature in the calculations...
    return (ccd, dccd)


def angular_diameter_distance_LCDM(z,
                                   dz=None,
                                   Omega_M=LCDM_PARAMS["Omega_M"],
                                   dOmega_M=LCDM_PARAMS["sigma_Omega_M"],
                                   H0=LCDM_PARAMS["H0"]/c_km_per_s,
                                   dH0=LCDM_PARAMS["sigma_H0"]/c_km_per_s):
    """Compute the angular diameter distance for a given z in flat LCDM universe.

    See comoving_coord_distance_LCDM() help for details.  Cosmological parameters
    default to WMAP 5yr results.

    If dz is None (default), errors will not be calculated and a float will be
    returned.  If it is not None, a tuple (value, error) will be returned.

    """
    Omega_L = 1.0-Omega_M
    ccd=_comoving_coord_distance_value(z,Omega_M,Omega_L,H0)

    # partial derivs of comoving-coord-dist X
    (dXdz,
     dXdH0,
     dXdOmega_M,
     dXdOmega_L)=_comoving_coord_distance_partials(z,Omega_M,Omega_L,H0,ccd)

    # compute angular diam dist
    D = ccd/(1+z)

    # and its uncertainty if requested
    if dz is not None:
        # now compute angular diameter distance partials
        dDdz = dXdz/(1+z) - ccd/((1+z)**2)
        dDdH0 = dXdH0/(1+z)
        dDdOmega_M = dXdOmega_M/(1+z)
        dDdOmega_L = dXdOmega_L/(1+z)

        sigma_D=np.sqrt(
            (dDdH0 * dH0) **2.0
            + (dDdz * dz) **2.0
            + ((dDdOmega_M-dDdOmega_L) * dOmega_M) ** 2.0)

    if dz is None:
        return D
    else:
        return(D, sigma_D)


def luminosity_distance_LCDM(z,
                             dz=None,
                             Omega_M=LCDM_PARAMS["Omega_M"],
                             dOmega_M=LCDM_PARAMS["sigma_Omega_M"],
                             H0=LCDM_PARAMS["H0"]/c_km_per_s,
                             dH0=LCDM_PARAMS["sigma_H0"]/c_km_per_s):
    """Compute the luminosity distance for a given z in flat LCDM universe.

    See comoving_coord_distance_LCDM() help for details.  Cosmological
    parameters default to WMAP 5yr results.

    If dz is None (default), errors will not be calculated and a float
    will be returned.  If it is not None, a tuple (value, error) will
    be returned.

    """
    Omega_L = 1.0-Omega_M
    ccd=_comoving_coord_distance_value(z,Omega_M,Omega_L,H0)

    # partial derivs of comoving-coord-dist X
    (dXdz,
     dXdH0,
     dXdOmega_M,
     dXdOmega_L)=_comoving_coord_distance_partials(z,Omega_M,Omega_L,H0,ccd)

    # compute luminosity dist
    D = ccd*(1+z)

    # and its uncertainty if requested
    if dz is not None:
        # now compute luminosity distance partials
        dDdz = dXdz*(1+z) + ccd
        dDdH0 = dXdH0*(1+z)
        dDdOmega_M = dXdOmega_M*(1+z)
        dDdOmega_L = dXdOmega_L*(1+z)

        sigma_D=np.sqrt(
            (dDdH0 * dH0) **2.0
            + (dDdz * dz) **2.0
            + ((dDdOmega_M-dDdOmega_L) * dOmega_M) ** 2.0)

    if dz is None:
        return D
    else:
        return(D, sigma_D)

def luminosity_distance_pen_approx(z,
                                   dz=None,
                                   Omega_M=LCDM_PARAMS["Omega_M"],
                                   dOmega_M=LCDM_PARAMS["sigma_Omega_M"],
                                   H0=LCDM_PARAMS["H0"]/c_km_per_s,
                                   dH0=LCDM_PARAMS["sigma_H0"]/c_km_per_s):
    """Compute luminosity distance using Pen (1999) approximation.

    Pen (1999) produced an analytical approximation to the luminosity
    distance for a flat pressureless matter cosmology with a cosmological
    constant.

    Errors not currently propagated. Returns (D, 0) for compatibility
    with the numerically integrated version.

    Claims:
     1) exact for Omega_M -> 1- or ->0+
     2) relative error -> 0 for z->inf
     3) For 0.2<Omega_M<1, rel error < 0.4%
     4) For any parameters, relative error < 4%

    Note: Pen calls Omega_M "Omega_0"

    """
    def _eta(a, Omega_0):
        """Conformal time approximation."""
        s=((1-Omega_0)/Omega_0)**(1./3)
        return 2*np.sqrt(s**3+1)*( a**-4 - 0.1540*(s/a**3) + 0.4304*(s/a)**2
                                   + 0.19097*(s**3/a) + 0.069941*s**4 )**-0.125
    a=1.0/(1+z)
    dL=(1.0/H0)*(1+z)*( _eta(1,Omega_M)-_eta(a, Omega_M) )
    return dL


#-----------------------------------------------------------
#
def luminosity_from_flux(S_Jy, z, alpha=0.0, useLCDM=True, **kwargs):
    """Compute specific luminosity given a flux density.

    Returns the K-corrected luminosity corresponding to the observed
    frequency in the emitter frame, assuming a power law spectrum with
    spectral index alpha.

    Flux density should be in Jy. Returns W/Hz.

    Extra arguments (**kwargs) passed along to
    luminosity_distance_LCDM. If you change H0, ensure it is in units
    of km/s/Mpc.

    Does *not* scale by 4pi (i.e., this is per steradian, not
    isotropic).

    If useLCDM is True (default), uses the LCDM cosmology. Otherwise
    uses a flat matter-dominated cosmology.

    """
    S_mks=S_Jy*1.0e-26
    if useLCDM:
        (r, r_err)=comoving_coord_distance_LCDM(z, **kwargs) # Mpc
    else:
        r=comoving_coord_distance_flatMD(z, **kwargs)
        r_err=0

    r *= 1e6 * 3.0857e16 # m
    r_err *= 1e6 * 3.0857e16 # m
    return r*r * (1+z)**(1-alpha) * S_mks

def isotropic_luminosity_from_flux(*args, **kwargs):
    """Compute isotropic specific luminosity given a flux density."""
    return 4*np.pi*luminosity_from_flux(*args, **kwargs)

#-----------------------------------------------------------
# Sidereal time conversions
def lst_to_datetime(lst, date, long=OVRO_LONGITUDE):
    """Convert local sidereal time to a datetime for a specific date.

    lst is an ephem.Angle object representing the LST to check (this can
    also be a string "hh:mm:ss.ssss")

    date is an datetime.date object specifying the date

    long is the longitude of interest.  Defaults to OVRO 40m's longitude
    of -118:16:53.83.

    """

    obs=ephem.Observer()
    obs.long=long
    obs.lat="00:00:00" # arbitrary
    obs.elevation=0

    obs.date=date.strftime("%Y/%m/%d") # sets to midnight on spec. date

    # The observer is now set to midnight on the date of interest.
    # Since LST is defined as when a particular RA transits the meridian,
    # we can find the UTC by finding when a body at arbitrary declination
    # transits our observer.
    return obs.next_transit(body_from_radec(lst,"00:00:00")).datetime()

def datetime_to_lst(datetime, long=OVRO_LONGITUDE):
    """Convert a datetime to LST."""
    obs=ephem.Observer()
    obs.long=ephem.degrees(long)
    obs.lat=ephem.degrees("45:00:00") # arbitrary
    obs.elevation=0
    obs.date=datetime.strftime("%Y/%m/%d %H:%M:%S")
    return obs.sidereal_time()


#-----------------------------------------------------------
# Convert a string or datetime to a #$@%@#$ ephem.date
def get_ephem_date(d):
    # convert time into an ephem.date
    if isinstance(d, ephem.date):
        return d
    else:
        dt = pcal.get_datetime(d)
        return ephem.date(dt.strftime('%Y/%m/%d %H:%M:%S'))

#-----------------------------------------------------------
# Compute airmass
def airmass(el):
    """Compute airmass at apparent elevation el (degrees).

    Uses the interpolation formula from Kasten and Young (1989,
    Applied Optics, 28, 22, 4375) which is very slightly different
    from just using sec(za).  Almost surely an unnecessary correction,
    but this should be more accurate for el<20 degrees.

    """
    a=0.50572
    b=6.07995 # degrees
    c=1.6364

    return 1.0/(np.sin(el*np.pi/180)+a*(el+b)**(-c))


#-----------------------------------------------------------
def beam_radec(target_radec, obs_time, ant=True, beam_sep=12.95, cos_el_thresh=0.1):
    """Compute the RA and DEC of the OVRO 40m antenna and reference beams.

    Computes the RA and DEC of the antenna and reference beams when
    observing a position (target_ra, target_dec) at specified
    obs_time, using the small-angle approximation for beam offsets.
    If ant is True, the target is placed in the antenna beam; if
    False, it is placed in the reference beam.  Beams are positioned
    beam_sep arcmin apart in azimuth with the reference beam in the
    negative cross-elevation direction from the antenna beam.

    Positions should be specified as ephem.angles, obs_time as a
    datetime.datetime object.

    Returns ((ant_ra, ant_dec), (ref_ra, ref_dec)) where all are
    ephem.angle objects (hours for ra, degrees for dec).

    Raises ValueError if the position is too close to the NCP for
    reliable results. If you really want to push this, you can adjust
    cos_el_thresh closer to 0 to allow higher declinations, but that's
    probably a bad idea.

    """
    (target_ra,target_dec)=target_radec # unpack

    if abs(np.cos(target_dec)) < cos_el_thresh:
        raise ValueError("Requested beam locations too near the pole. Adjust cos_el_thresh if you know what you're doing and really want to do this, but beam locations are probably inaccurate.")
    # convert from arcmin to an ephem.angle
    beam_sep = ephem.degrees(beam_sep*np.pi/180.0/60.0)

    par_ang = parallactic_angle(target_ra, target_dec, date=obs_time)
    if ant:
        ant_ra = target_ra
        ant_dec = target_dec
        ref_ra = target_ra + beam_sep*np.cos(par_ang)/np.cos(target_dec)
        ref_dec = target_dec - beam_sep*np.sin(par_ang)
    else:
        ant_ra = target_ra - beam_sep*np.cos(par_ang)/np.cos(target_dec)
        ant_dec = target_dec + beam_sep*np.sin(par_ang)
        ref_ra = target_ra
        ref_dec = target_dec

    return ((ephem.hours(ant_ra),ephem.degrees(ant_dec)),
            (ephem.hours(ref_ra),ephem.degrees(ref_dec)))


def beam_radec2(target_radec, obs_time, ant=True, beam_sep=12.95):
    (ra,dec)=target_radec
    # convert from arcmin to an ephem.angle
    beam_sep = ephem.degrees(beam_sep*np.pi/180.0/60.0)

    az,el = ovro40m_azel2(ra,dec, obs_time)
    az_offset = beam_sep / np.cos(el)
    if ant:
        ant_az=az
        ref_az=az-az_offset
    else:
        ant_az=az+az_offset
        ref_az=az

    (ant,ref)= (radec_from_azel2(ant_az, el, obs_time),
                radec_from_azel2(ref_az, el, obs_time))

    if airmass(el) > 3:
        print (el)

    if ref[1] < 0.24:
        print ('fnop')

    return ant,ref



def polar_to_cartesian_astro(radius, theta):
    """Convert polar coordinates to cartesian x/y coordinates.

    Uses the astronomy convention that angles are measured north through east,
    with 0 degrees due north. Positive x direction is east, positive y is north.

    Theta is measured in degrees.

    """
    # angle is north through east
    x=radius*np.sin(theta*np.pi/180.)
    y=radius*np.cos(theta*np.pi/180.)

    return x,y


def brightness_temp_vlbi(Sc, th_maj, th_min, freq, z):
    """Compute brightness temperature of gaussian core.

    Sc - core flux density (Jy)
    th_maj, th_min - major and minor axes of gaussian (mas)
    freq - observation frequency (GHz)
    z - redshift
    """