""" Telescope constants used in the scheduler.

Check the telescope_data.ipynb notebook for more info.

TODO:
- Replace with objects instead of constants (requires code refactor)
- Check move_in_azimuth
- Change za to el (?)
"""

from numpy import sign

# Telescope coordinates
lon = -(118 + 16/60 + 53.83/3600)   # degrees
lat = 37 + 13/60 + 55.7/3600        # degrees
elevation = 1236                    # meters

# Slew parameters
# Old control system
slew_speed = 4.58 / 3600    # hours/degree
t_delay = 11.22 / 3600      # hours
# Based on data from Jan 2015 until July 2015
# v_X: speed in deg/sec
# t_X: delay in sec
v_el, t_el = 0.235347359169, 13.0936769942      # Elevation
v_az_p, t_az_p = 0.220658677426, 13.1007239104  # Azimuth, positive
v_az_n, t_az_n = 0.216049390308, 13.4735934952  # Azimuth, negative

# Measurement time
# July 14, 2015
time_doflux = 65.2 / 3600   # hours
time_dopoint = 575 / 3600   # hours
time_cal = 575 / 3600       # hours

# Source and region observability limits
zmin_s, zmax_s = 20, 60
zmin_r, zmax_r = 20, 60
# For sources below low declination limit
zmax_low_dec_s = 65
zmax_low_dec_r = 65
# Minimum source declination
dec_min = -20.0

# Sun and Moon avoidance
sun_avoid_angle = 10    # degree
moon_avoid_angle = 10   # degree

# Telescope slew model and observing time
def obs_time(type):
    """ Returns the time on source as function of type.

    Parameters
    ----------
    type : int
        Type of observation as defined on point key of sources list

    Returns
    -------
    obs_time : int
        Time for that observation in hours

    Notes
    -----
    It is not quite complete since I still have to add the capability for 
    handling low flux pointing calibrators that have a longer integrattion time
    than normal ones

    Add capability to handle cases not in the current list
    """
    if type == 0:   # normal source
        return time_doflux
    elif type == 1: # First pointing calibrator
        return time_dopoint
    elif type == 2: # Pointing calibrator
        return time_dopoint
        
    return

def move_in_azimuth(az_i, az_f): 
    """ Slew time with azimuth wrap.
    
    Starts with az_i in range [-89, 335] and az_f in [0, 360] and chose
    shortest posible path between them considering the azimuth wrap. It returns
    az_f in the range [-89, 335]

    Parameters
    ----------
    az_i : float
        Initial azimuth
    az_f : float
        Final azimuth

    Returns
    -------
    az_f : float
        Final azimuth with wrap

    Notes
    -----
    The telescope azimuth goes from -89 to 335 degrees.

    If we are using az in the [0, 360] range we have an ambiguity when the az
    is in the range [271, 335]

    In summary:
    First turn [0, 271] and [335, 360]
    Second or first turn [271, 335]
    
    wmax, March 22, 2012
    """
    # Check that coordinates make sense
    assert (-89 <= az_i <= 335),\
           f'za_i = {az_i}, not in az wrap coordinates [-89, 335]'
    assert (0 <= az_f <= 360),\
           f'za_f = {az_f}, not in az geometric coordinates [0, 360]'
    # Puts az_f in its possible values depending on azimuth wrap
    # First turn with no ambiguity
    if (0 <= az_f <= 271):
        return az_f
        # az_turn = 1
    # First turn negative azimuth
    elif (335 <= az_f <= 360):
        return az_f - 360
        # az_turn = 1
    # First of second turn depending on which is closer
    else:
        # First turn
        az_f_first = az_f - 360
        d_az_first = abs(az_f_first - az_i)
        # Second turm
        az_f_second = az_f
        d_az_second = abs(az_f_second - az_i)
        # Find shortest
        if d_az_first <= d_az_second:       
            return az_f_first
            # az_turn = 1
        else:
            return az_f_second
            # az_turn = 2

def slew_time(za_i, az_i, za_f, az_f):
    """ Given initial and final telescope position returns the slew time in 
    hours.

    The positions already have the information on the azimuth wrap, i.e. they
    are in the range [-89, 335] which the telescope really has.

    The telescope model is defined in usints of deg/sec and sec but the final 
    slewing time should be in hours.

    Parameters
    ----------
    za_i : float
        Zenith angle of first source in decimal degrees
    az_i : float
        Azimuth of first source in decimal degrees
    za_f : float
        Zenith angle of second source in decimal degrees
    az_f : float
        Azimuth of second source in decimal degrees

    Returns
    -------
    slew_time : float
        Slewing time between sources in hours

    Notes
    -----
    wmax, March 22, 2012
    """
    # Zenith distance
    d_za = abs(za_f - za_i)
    # Azimuth distance
    d_az = az_f - az_i
    # Get slew times for zenith. Same speed in both senses.
    za_slew_time_secs = d_za / v_el + t_el
    # Get slew time for azimuth. Sense dependent model.
    if d_az >= 0:
        # Azimuth slew in positive direction
        az_slew_time_secs = d_az / v_az_p + t_az_p
    else:
        # Azimuth slew in negative direction
        az_slew_time_secs = abs(d_az) / v_az_n + t_az_n
    # The longest is returned
    slew_time_secs = max(az_slew_time_secs, za_slew_time_secs)   
    slew_time_hours = slew_time_secs / 3600
    return slew_time_hours

def slew_time_nowrap(za_i, az_i, za_f, az_f):
    """ Given initial and final telescope position returns the slew time in 
    hours

    The telescope model is defined in usints of deg/sec and sec but the final 
    slewing time should be in hours.

    This model has different speeds for different axis and directions in the
    case of the azimuth. It's still using sky coordinates instead of telescope
    coordinates so it's still inaccurate in that sense.

    Parameters
    ----------
    za_i : float
        Zenith angle of first source in decimal degrees
    az_i : float
        Azimuth of first source in decimal degrees
    za_f : float
        Zenith angle of second source in decimal degrees
    az_f : float
        Azimuth of second source in decimal degrees

    Returns
    -------
    slew_time : float
        Slewing time between sources in hours

    Notes
    -----
    August 16, 2010
    This new model has different speeds for different axis and directions in the 
    case of the azimuth.
    It still using sky coordinates instead of telescope coordinates so it still 
    inacurate in that sense.
    """
    # zenith distance
    d_za = abs(za_f - za_i)
    # azimuth distance
    dd_az = az_f - az_i
    # the distance has to be less than 180 deg
    if abs(dd_az) > 180:
        d_az = sign(-dd_az)*(360 - abs(dd_az))    
    else:
        d_az = dd_az
    # Check which axis dominates slew
    if d_za > abs(d_az):
        # Zenith angle slew
        slew_time_secs = d_za / v_el + t_el
    elif d_az >= 0:
        # Azimuth slew in positive direction
        slew_time_secs = d_az / v_az_p + t_az_p
    elif d_az < 0:
        # Azimuth slew in negative direction
        slew_time_secs = abs(d_az) / v_az_n + t_az_n
    slew_time_hours = slew_time_secs / 3600
    return slew_time_hours

def slew_time_old(za_i, az_i, za_f, az_f):
    """ Given initial and final telescope position returns the slew time in 
    hours

    Parameters
    ----------
    za_i : float
        Zenith angle of first source in decimal degrees
    az_i : float
        Azimuth of first source in decimal degrees
    za_f : float
        Zenith angle of second source in decimal degrees
    az_f : float
        Azimuth of second source in decimal degrees

    Returns
    -------
    slew_time : float
        Slewing time between sources in hours

    Notes
    -----
    This model was used before August 16, 2010
    """
    # zenith distance
    d_za = abs(za_i - za_f)
    # azimuth distance
    dd_az = abs(az_i - az_f)
    # the distance has to be less than 180 deg
    if dd_az > 180:
        d_az = 360 - dd_az 
    else:
        d_az = dd_az
    return max(d_za, d_az) * slew_speed + t_delay