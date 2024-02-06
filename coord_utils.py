"""
Basic functions for coordinate/time transformations and source properties
    - Observer is OVRO 40m telescope
    - Assumes all the coordinates are J2000.

TODO:
    - change obs_range to tuple
    - change how check_observability checks range
    - change closest_range to tuple
    - make lst_range_to_obs_range return tuples
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
        np.arcsin(np.cos(lha) * np.cos(dec) * np.cos(lat) + \
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
    # Create and array of sample lsts to calculate position
    # resolution of 5 minutes
    n_lst = np.ceil(24 * 60 / t_resol)
    lsts = np.linspace(0, 24, n_lst)
    # Calculate the position of the source on that sample lsts
    ZA = [radec_zaaz(ra, dec, lst)[0] for lst in lsts]
    AZ = [radec_zaaz(ra, dec, lst)[1] for lst in lsts]
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
    # Generate observability ranges from the observable vector
    obs_range = obsrange_from_observable(observable, lsts)
    # Plot the results for inspection
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
        List with pair of values with the extremes of the intervals.
    lst_obs : float
        Decimal hours

    Returns
    -------
    observable : bool
        True if all objects are observable
    """
    # wrap around to [0, 24) range
    lst_obs = lst_obs % 24
    # Loop through the observing ranges and check is lst_obs is in at least one
    for range in obs_range:
        if range[0] <= lst_obs <= range[1]:
            return True
    return False

def closest_obsrange(obs_range, lst_obs):
    """Determine the closest observing range

    This is used to determine the observing range that goes into schedules
    start and stop ranges.

    Given an observing range tries:

    1) See if lst_obs is inside one of the ranges in obs_range and return that

    2) If 1) doesn't work, finds the range that is inmediately after
    returns it

    Parameters
    ----------
    obs_range : list
        List with pair of values with the extremes of the intervals.
    lst_obs : float
        Decimal hours

    Returns
    -------
    closest_range : list
        List with a pair of values containing the closest observing range

    Notes
    -----
    Be careful about ranges that finish in 24.00, those might continue on
    0.0 and went for some extra time

    A simmilar situation for the case when the range starts at 0.0.
    """
    # wrap around to [0, 24) range
    lst_obs = lst_obs % 24
    # ----------------------------------------
    def extend_range(obs_range, range):
        """Internal function

        This function extends a range if there is a contiguous one
        """
        # lower and upper range limits
        range_l = range[0]
        range_u = range[1]
        # check if window continues at 00 but doesn't start at 00
        if range_u == 24 and range_l != 0:
            # loop through ranges
            for alt_range in obs_range:
                # if range starting at 00, join with previous range
                if alt_range[0] == 0:
                    range_u = alt_range[1]
        # check if window continues before 00 but doesn't end at 24.0
        elif range_l == 0 and range_u != 24:
            # loop through ranges
            for alt_range in obs_range:
                # if range finishing at 24.0, join with previous range
                if alt_range[1] == 24:
                    range_l = alt_range[0]
        return range_l, range_u
    # ----------------------------------------
    # try to find range that contains lst_obs
    for range in obs_range:
        # lst in range
        if range[0] <= lst_obs and lst_obs <= range[1]:
            # extend the range if necessary
            range_l, range_u = extend_range(obs_range, range)
            # save range
            closest_range = [range_l, range_u]
            return closest_range
    # if lst_obs not inside a given range find next range
    for range in obs_range:
        # this is the next range
        if lst_obs <= range[0]:
            # extend the range if necessary
            range_l, range_u = extend_range(obs_range, range)
            # save range
            closest_range = [range_l, range_u]
            return closest_range
    # if lst_obs doesn't have next range it could have previous
    for range in obs_range:
        # this is the previous range
        if range[1] <= lst_obs:
            # extend the range if necessary
            range_l, range_u = extend_range(obs_range, range)
            # save range
            closest_range = [range_l, range_u]
            return closest_range
    # default return value
    return []

def wait_time_for_observability(obs_range, lst_obs):
    """ Return the time to wait before source is on observability range.

    Parameters
    ----------
    obs_range : list
        List with pair of values with the extremes of the intervals.
    lst_obs : float
        Decimal hours

    Returns
    -------
    d_lst_wait : float
        Wait time returned is in decimal hours
    """
    # wrap around for [0, 24) range
    lst_obs = lst_obs % 24
    # check if source is observable, if so return 0
    if check_observability(obs_range, lst_obs):
        d_lst_wait = 0.0
    # find the wait time
    else:
        # start time is first element of closest obsrange
        lst_start = closest_obsrange(obs_range, lst_obs)[0]
        # two cases for wait time, be careful
        if lst_obs <= lst_start:
            d_lst_wait = lst_start - lst_obs
        else:
            d_lst_wait = (24.0 - lst_obs) + lst_start
    return d_lst_wait

def intersect_obs_ranges(obs_ranges, t_resol=0.01):
    """ Given a list of observing ranges, calculate their intersection.

    Parameters
    ----------
    obs_range : list
        List with pair of values with the extremes of the intervals.
    t_resol : float, optional
        Time resolution in minutes. Its default value is 0.01 minutes.

    Returns
    -------
    obs_range : list
        list with a pair of values with the extremes of the intervals.
    """
    # make lsts vector with given time resolution
    n_lst = np.ceil(24 * 60 / t_resol)
    lsts = np.linspace(0, 24, n_lst)
    # calculate observable vectors for intersection
    observable = []
    for lst in lsts:
        # initialize observability flag, observable unlees opposite proved
        observable_at_lst = True
        # goes through the list of observing ranges
        for obs_range in obs_ranges:
            # check observability for particular range
            if not check_observability(obs_range, lst):
                observable_at_lst = False
        # Test abservability for all regions
        if observable_at_lst:
            observable.append(1)
        else:
            observable.append(0)
    # get the combined obs_range
    obs_range = obsrange_from_observable(observable, lsts)
    return obs_range

def lst_range_to_obs_range(lst_start, lst_stop):
    """ Given an lst_start an lst_stop returns the equivalent obs_range

    Parameters
    ----------
    lst_start, lst_stop : float
        LST range in hours.

    Returns
    -------
    obs_range : list
        list with pair of values with the extremes of the intervals.
    """
    # convert to [0, 24) range
    lst_start = lst_start % 24.0
    lst_stop = lst_stop % 24.0
    # two cases for obs_range
    if lst_start <= lst_stop:
        obs_range = [[lst_start, lst_stop]]
    else:
        obs_range = [[0.0, lst_stop], [lst_start, 24.0]]
    return obs_range

def dec_to_sex(dec, out_type='list', dec_plac=0, floor=True, time=True):
    """ Convert a decimal number to sexagesimal format

    Parameters
    ----------
    dec : float
        Decimal number
    out_type : string, default='list'
        Type of output required:
            out_type = 'list',    [dd, mm, ss]
            out_type = 'string'   'dd:mm:ss'
    dec_plac : int, default=0
        Number of decimal places for seconds
    floor : bool, default=True
        Take the floor for seconds if true
    time : bool, default=True
        Set 0-23 range for hh

    Returns
    -------
    sex : list or string
        Number in the specified sexagesimal format
    """
    # determine sign of the number
    if dec < 0:
        sign = -1
        dec = -dec
    else:
        sign = 1
    # calculate components
    dd = np.floor(dec)
    mm = np.floor( (dec - dd) * 60 )
    ss = (dec - dd - mm / 60.0) * 3600
    # Take the floor for seconds
    if floor:
        ss = np.floor(ss)
    # Use 0-23 range for time
    if time:
        dd = dd % 24
    # Check the output type
    if out_type == 'string':
        # format depend on number of decimal places
        if dec_plac == 0:
            format = '%02d:%02d:%02.0f'
        else:
            format =\
                '%02d:%02d:%0' + str(dec_plac + 3) + '.' + str(dec_plac) + 'f'
        sex = format % (dd, mm, ss)
        # add sign if necessary
        if sign == -1:
            sex = '-' + sex
    elif out_type == 'list':
        sex = [dd, mm, ss]
        # add sign if necessary
        if sign == -1:
            if dd != 0:
                sex[0] = -sex[0]
            elif mm != 0:
                sex[1] = -sex[1]
            else:
                sex[2] = -sex[2]
    return sex

def sex_to_dec(sex):
    """ Convert a sexagesimal number to decimal

    Parameters
    ----------
    sex : list
        List containing the sexagesimal components of the number

    Returns
    -------
    dec : float
        Number converted to decimal.

    Notes
    -----
    To make the conversion consider that the zero goes on the first nonzero
    component of the sexigesimal representation of the number.
    """

    # get the negative sign if any
    if sex[0] < 0 or sex[1] < 0 or sex[2] < 0:
        sign = -1
    else:
        sign = 1
    dec = sign * (abs(sex[0]) + abs(sex[1]) / 60 + abs(sex[2]) / 3600)
    return dec

def parallactic_angle_lst(ra, dec, lst, lat=tel_data.lat):
    """Compute the parallactic angle given the LST and latitude.

    Sign convention is that parallactic angle > 0 when an object is
    west of the meridian.

    Parameters
    ----------
    ra, dec : float
        RA and Dec in degrees.
    lst : float
        LST in decimal hours.
    lat: float, optional
        Telescope latitude in decimal degrees. Defaults to the value set in the
        telescope_data module.

    Returns
    -------
    pa_deg : float
        Angle in degrees.
    """
    #first need to convert everything to radians
    ra_rad = ra * deg_to_rad
    dec_rad = dec * deg_to_rad
    #get lst to right range
    while(lst >= 24):
        lst -= 24
    lst_rad = lst * 15 * deg_to_rad
    lat_rad = lat * deg_to_rad
    ha = lst_rad - ra_rad
    pa = np.arctan2(np.cos(lat_rad) * np.sin(ha),
                   np.cos(dec_rad) * np.sin(lat_rad)
                   - np.sin(dec_rad) * np.cos(lat_rad) * np.cos(ha))
    #then convert back to degrees for returning
    pa_deg = pa * rad_to_deg
    return pa_deg