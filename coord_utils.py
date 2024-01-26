"""
Basic functions for coordinate/time transformations and source properties
    - Observer is OVRO 40m telescope
    - Assumes all the coordinates are J2000.

TODO:
    - change obs_range to tuple
    - change how check_observability checks range
"""
import telescope_data as tel_data
import numpy as np
import ephem

rad_to_deg = 180 / np.pi
deg_to_rad = np.pi / 180

def radec_zaaz(ra, dec, lst):
    """Zenith/Azimuth for a source given RA/DEC and observer position.

    Parameters
    ----------
    ra : float
        RA decimal degrees
    dec : float
        DEC decimal degrees
    lst : float
        Decimal hours

    Returns
    -------
    za, az : float
        ZA/AZ coordinates decimal degrees

    Notes
    -----
    This simple recipe doesn't consider observer elevation nor refraction
    effects, etc. It should be off by some arcminutes, this is probably not a
    problem for the optimizer/scheduler. Of course it is not good enough to
    point a telescope!
    """
    # Wrap around to [0, 24) range
    lst = lst % 24
    # Observer coordinates
    long = tel_data.lon
    lat = tel_data.lat
    # Local Hour Angle.
    # 15 on ra is to convert to decimal hours
    # 15 on total converts from hours to degrees
    lha = (lst - ra / 15) * 15
    # convert everything to radians
    long = long * np.pi / 180
    lat = lat * np.pi / 180
    lha = lha * np.pi / 180
    dec = dec * np.pi / 180
    # za in radians
    za = np.pi / 2 - \
        np.asin(np.cos(lha) * np.cos(dec) * np.cos(lat) + \
                      np.sin(dec) * np.sin(lat))
    # az in radians. -180 < atan2 < 180. But 0 < az < 360
    az = np.arctan2(-np.sin(lha) * np.cos(dec), \
                          -np.sin(lat) * np.cos(lha) * np.cos(dec) + \
                          np.sin(dec) * np.cos(lat))
    # to degrees
    za = za * 180 / np.pi
    az = az * 180 / np.pi
    # 0 < az < 360
    if az < 0:
        az = az + 360
    return za, az

def dis_zaaz(za1, az1, za2, az2):
    """Angular distance in ZAAZ coordinates.

    Parameters
    ----------
    za1, az1 : float
        Source coordinates in decimal degrees
    za2, az2 : float
        Target coordinates in decimal degrees

    Returns
    -------
    dis : float
        Angular distance in degrees
    """
    # use dot product of two unitary vectors
    # convert to spherical coordinates and radians
    theta1 = za1 * np.pi / 180
    phi1 = (360 - az1) * np.pi / 180
    theta2 = za2 * np.pi / 180
    phi2 = (360 - az2) * np.pi / 180
    # calculate vectors
    r1 = np.array([np.sin(theta1) * np.cos(phi1),\
                    np.sin(theta1) * np.sin(phi1),\
                    np.cos(theta1)])
    r2 = np.array([np.sin(theta2) * np.cos(phi2),\
                    np.sin(theta2) * np.sin(phi2),\
                    np.cos(theta2)])
    # get angle between vectors in degrees
    dis = np.acos(min(np.sum(r1 * r2), 1) ) * 180 / np.pi
    return dis

def dis_radec(ra1, dec1, ra2, dec2):
    """Angular distance in RADEC coordinates.
    ra, dec and final distance are in degrees

    Parameters
    ----------
    ra1, dec1 : float
        Source coordinates in decimal degrees
    ra2, dec2 : float
        Target coordinates in decimal degrees

    Returns
    -------
    dis : float
        Angular distance in degrees
    """
    # convert to radians
    ra1 = ra1 * np.pi / 180.0
    dec1 = dec1 * np.pi / 180.0
    ra2 = ra2 * np.pi / 180.0
    dec2 = dec2 * np.pi / 180.0
    # calculate the vectors
    r1 = np.array([np.cos(dec1) * np.cos(ra1),\
                    np.cos(dec1) * np.sin(ra1),\
                    np.sin(dec1)])

    r2 = np.array([np.cos(dec2) * np.cos(ra2),\
                    np.cos(dec2) * np.sin(ra2),\
                    np.sin(dec2)])
    # get angle between vectors in degrees
    dis = np.acos( min(np.sum(r1 * r2), 1) ) * 180. / np.pi
    return dis

def jd_to_djd(jd):
    """Convert from JD to DJD

    Parameters
    ----------
    jd : float
        Julian Date
    
    Returns
    -------
    djd : float
        Dublin Julian Date
    """
    djd = jd - 2415020
    return djd

def sun_pos(jd):
    """Sun position in ZAAZ coordinates given JD.

    Parameters
    ----------
    jd : float
        Julian Date
    
    Returns
    -------
    ra, dec : float
        RA/DEC coordinates decimal degrees

    Notes
    -----
    Uses PyEphem to get position.
    """
    # convert to Dublin Julian Day as used by PyEphem
    djd = jd_to_djd(jd)
    # get postion of Sun
    Sun = ephem.Sun(djd)
    # Angles in decimal degrees for ra and dec
    ra = Sun.ra * rad_to_deg
    dec = Sun.dec * rad_to_deg
    return ra, dec

def moon_pos(jd):
    """Moon position in ZAAZ coordinates given JD.

    Parameters
    ----------
    jd : float
        Julian Date
    
    Returns
    -------
    ra, dec : float
        RA/DEC coordinates decimal degrees

    Notes
    -----
    Uses PyEphem to get position.
    """
    # convert to Dublin Julian Day as used by PyEphem
    djd = jd_to_djd(jd)
    # get postion of Moon
    Moon = ephem.Moon(djd)
    # Angles in decimal degrees for ra and dec
    ra = Moon.ra * rad_to_deg
    dec = Moon.dec * rad_to_deg
    return ra, dec

def obsrange_from_observable(observable, lsts):
    """ Take the observable and lsts vectors and returns the obs_range

    Parameters
    ----------
    observable : array_like
        list with 0/1 that represents the observability for a given sky
        position (ra/dec) at the corresponding lst in vector lsts
    lsts : array_like
        Decimal hours

    Returns
    -------
    obs_range : list
        list with pair of values with the extremes of the intervals.
    """
    # Generate observability ranges from the observable vector 
    obs_range = []
    i = 0
    N = len(lsts)
    while i < N:
        # Find first lst of observable range
        if observable[i] == 1:
            # lower limit for range
            lst_l = lsts[i]
            # default upper limit, in case it is last element
            lst_u = lsts[i]
            i += 1
            while i < N:       
                # Find first zero element
                if observable[i] == 0:
                    # upper limit is previous element
                    lst_u = lsts[i - 1]
                    # save observability range
                    obs_range.append([lst_l, lst_u])
                    break
                # If last element on array close obs_range
                elif i == N - 1:
                    # upper limit
                    lst_u = lsts[i]
                    # save observability range
                    obs_range.append([lst_l, lst_u])
                # try next element
                i += 1
        # try next element
        i += 1
    return obs_range

def source_obsrange(ra, dec, zmin, zmax, azmin=0, azmax=360, t_resol=5):
    """ Calculate the range of lst for which a given sky position is in the 
    range zmin < za < zmax.

    Parameters
    ----------
    ra, dec : float
        RA/DEC coordinates decimal degrees
    zmin, zmax : float
        Minimum/maximum zenith angle
    azmin, azmax : float, optional
        Minimum/maximum azimuth angle
    t_resol : float, optional
        Time resolution in minutes

    Returns
    -------
    obs_range : list
        list with pair of values with the extremes of the intervals.
    """  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Create and array of sample lsts to calculate position
    # resolution of 5 minutes
    n_lst = np.ceil(24 * 60 / t_resol)
    lsts = np.linspace(0, 24, n_lst)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Calculate the position of the source on that sample lsts
    ZA = [radec_zaaz(ra, dec, lst)[0] for lst in lsts]
    AZ = [radec_zaaz(ra, dec, lst)[1] for lst in lsts]
    #~~~~~~~~~~~~~~~~~~~~~~~
    # Generate a vector with:
    # 0: Source is not observable for that lst
    # 1: Source is observable
    observable = []
    for i in range(len(ZA)):
        za = ZA[i]
        az = AZ[i]

        if (zmin <= za and za <= zmax) and (azmin <= az and az <= azmax):
            observable.append(1)
        else:
            observable.append(0)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Generate observability ranges from the observable vector  
    obs_range = obsrange_from_observable(observable, lsts)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # plot the results for inspection
    # print obs_range
    # pylab.plot(lsts, za)
    # pylab.plot(lsts, numpy.array(observable) * 60.)
    # pylab.show()    
    return obs_range

def check_observability(obs_range, lst_obs):
    """ Check if a given lst is contained on obs_range.

    Parameters
    ----------
    obs_range : list
        list with pair of values with the extremes of the intervals.
    lst_obs : float
        Decimal hours

    Returns
    -------
    observable : bool
        True if all objects are observable
    """
    # wrap around to [0, 24) range
    lst_obs = lst_obs % 24
    # default observability value
    observable = False
    # Loop through the observing ranges and check is lst_obs is in at least one
    for range in obs_range:
        if range[0] <= lst_obs and lst_obs <= range[1]:
            observable = True
            break
    return observable