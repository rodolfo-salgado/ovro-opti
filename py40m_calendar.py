""" Calendar and clock related routines.

TODO:
- Update docs
- Update style
"""

import datetime
from time import mktime
import dateutil.parser as parser
import numpy

#-----------------------------------------------------------
def julian_to_date(j):
    """Convert Julian day number to a date.

    Follows section 4 of Hatcher 1984 QJRAS 25, 53--55.
    Does *not* work properly with fractional day numbers.
    """
    g=int(int((j-4479.5)/36524.25) * 0.75 + 0.5) - 37
    N=j+g
    A=int(N/365.25)-4712
    dp=int(numpy.mod(N-59.25, 365.25)) # d'
    M=numpy.mod(int((dp+0.5)/30.6)+2, 12)+1 # 1=Jan
    D=int(numpy.mod(dp+0.5, 30.6))+1
    return datetime.date(A, M, D)

#-----------------------------------------------------------
def mjd_to_date(mjd):
    """Convert modified Julian day to a date.

    The algorithm in julian_to_date does not handle fractional days
    properly, so we drop the fractional part explicitly.  Ignore the
    0.5 that seems to be missing in the conversion to JD, it has no
    bearing on the conversion of the day part of the date.
    """
    return julian_to_date(int(mjd+2400001))

#-----------------------------------------------------------
def mjd_to_datetime(mjd):
    """Convert MJD to a date/time using fractional MJD.

    Does not handle leap seconds.
    """
    # the numpy.floor() should be unnecessary on next line, included
    # it for clarity
    whole_days=int(numpy.floor(mjd))
    frac_days=mjd-int(mjd)

    h = frac_days*24.0
    frac_h = h-int(h)
    h = int(h)
    m = frac_h*60.0
    frac_m = m-int(m)
    m = int(m)
    s=int(round(frac_m*60.0))

    # Ok, we might have rounded to 60 seconds.  Correct this and
    # propagate up as necessary.  This does not handle leap seconds.
    while s >= 60:
        s-=60
        m+=1
        while m >= 60:
            m-=60
            h+=1
            while h >= 24:
                h-=24
                whole_days+=1

    date=mjd_to_date(whole_days)
    time=datetime.time(h, m, s)

    return datetime.datetime.combine(date, time)


#-----------------------------------------------------------
def date_to_julian(dt):
    """Convert a datetime.date object to a Julian day.

    Follows section 3 of Hatcher 1984 QJRAS 25, 53--55.
    """
    A=dt.year
    M=dt.month # 1-12
    D=dt.day

    Aprime = A-int((12.0-M)/10.0)
    Mprime = (M-3)%12

    y=int(365.25*(Aprime+4712))
    d=int(30.6*Mprime + 0.5)
    N=y+d+D+59

    g=int( 0.75 * int((Aprime/100.0)+49) ) - 38
    return N-g

#-----------------------------------------------------------
def date_to_mjd(dt):
    """Convert a datetime.date object to an MJD day."""
    return date_to_julian(dt)-2400001

#-----------------------------------------------------------
def mjd_to_jd(mjd):
    return mjd+2400000.5

#-----------------------------------------------------------
def jd_to_mjd(jd):
    return jd-2400000.5

#-----------------------------------------------------------
def jd_to_datetime(jd):
    return mjd_to_datetime(jd_to_mjd(jd))

#-----------------------------------------------------------
def datetime_to_mjd(dt):
    """Convert a datetime.datetime object to a fractional MJD date."""
    d = date_to_mjd(dt)
    if isinstance(dt, datetime.datetime):
        delta_t = dt-dt.replace(hour=0, minute=0, second=0, microsecond=0)
        if delta_t.days != 0:
            raise Exception('More than a day since midnight!')
        d += (delta_t.seconds + delta_t.microseconds/1.0e6)/86400.0
    return d

#-----------------------------------------------------------
def timestr_to_datetime(s):
    """Try to convert a string to a datetime."""
    dt=parser.parse(s)
    return(dt)

#-----------------------------------------------------------
def timestr_to_mjd(s):
    """Try to convert a string to MJD.  You probably want to use get_mjd() instead."""
    dt=parser.parse(s)
    return(datetime_to_mjd(dt))

#-----------------------------------------------------------
def get_mjd(d, minimum_mjd=10000):
    """Figure out the MJD for date d.  Always return a float.

    d can be an integer or float MJD, a string MJD, a parseable date
    string, or a datetime.date or datetime.datetime object.

    If input is a numerical value and < minimum_mjd, it will be
    interpreted as a decimal year.

    """
    check_min=False

    if isinstance(d,float) or isinstance(d,int):
        check_min=True
        mjd=float(d)

    elif isinstance(d,str):
        try:
            mjd=float(d)
        except ValueError:
            mjd=datetime_to_mjd(parser.parse(d))
        else:
            check_min=True

    elif isinstance(d, datetime.date):  # datetime.datetime is also a datetime.date so this handles both
        return float(datetime_to_mjd(d))

    if check_min and mjd<minimum_mjd: # interpet it as a decimal year if so
        year=int(mjd) # get the year
        this_year=datetime.datetime(year,1,1)
        year_length=datetime.datetime(year+1,1,1)-this_year # find number of days (to account for leap yrs)
        mjd=datetime_to_mjd(this_year+datetime.timedelta(days=((mjd-year)*year_length.days)))

    return mjd


#-----------------------------------------------------------
def get_datetime(d):
    """Figure out the datetime for a date d.

    d can be a date, datetime, or string that can be parsed by timestr_to_datetime().

    """
    if isinstance(d,datetime.datetime):
        return d
    elif isinstance(d,datetime.date):
        return datetime.datetime.combine(d,datetime.time(0))
    elif isinstance(d, str):
        return timestr_to_datetime(d)

#-----------------------------------------------------------
def parse_hms_time(s):
    """Parse a time string in either hh:mm:ss or hh:mm format.

    Raises a ValueError if it fails.
    """
    try:
        time = datetime.datetime.strptime(s, '%H:%M:%S').time()
    except ValueError:
        time = datetime.datetime.strptime(s, '%H:%M').time()

    return time

#-----------------------------------------------------------
def parse_dump_timestamp(s):
    """Parse a date/time string from a cmbprog dump file.

    Format is yyyy/ddd/hhmmss where ddd is the day of the year.
    """
    (year,doy,tm) = (s.strip()).split('/')
    year=int(year)
    doy=int(doy)

    # convert DoY to real date
    date = datetime.date(year, 1, 1) + datetime.timedelta(days=(doy-1))

    hh = int(tm[0:2])
    mm = int(tm[2:4])
    ss = int(tm[4:6])

    time = datetime.time(hh, mm, ss)

    dt = datetime.datetime.combine(date, time)
    return dt


#-----------------------------------------------------------
def parse_40m_timestamp(s, year):
    """Parse a date/time string from a 40m telescope log file.

    Format is ddd/hhmmss where ddd is the day of the year.
    """
    (doy, tm) = (s.strip()).split('/')
    doy=int(doy)

    # convert DoY to real date
    date = datetime.date(year, 1, 1) + datetime.timedelta(days=(doy-1))

    hh = int(tm[0:2])
    mm = int(tm[2:4])
    ss = int(tm[4:6])

    time = datetime.time(hh, mm, ss)

    dt = datetime.datetime.combine(date, time)
    return dt


#-----------------------------------------------------------
def parse_csvfile_timestamp(s):
    """Parse date/time string from our standard CSV file format.

    Format is YYYY-MM-DD HH:MM:SS.
    """
    return datetime.datetime.strptime(s, '%Y-%m-%d %H:%M:%S')

#-----------------------------------------------------------
def parse_sql_timestamp(s):
    """Parse the SQL format time string into a datetime.

    Currently this is the same as our csvfile standard.

    """
    return parse_csvfile_timestamp(s)

#-----------------------------------------------------------
def date_to_40m_log(d):
    """Convert datetime into name for a 40m log file for that day.

    Assumes you've got the right day (remember it's in UTC!)
    Accepts datetimes just as well.

    """
    return d.strftime('%d%b%y').lower()

def datetime_to_40m_log(d):
    """Passes through to date_to_40m_log.

    Don't use this, date_to_40m_log works just fine for datetimes.
    This is just here for backward compatibility.

    """
    return date_to_40m_log(d)

#-----------------------------------------------------------
def parse_40m_log_name(s):
    """Convert 40m log format date (eg, 13jan09) into datetime.date object."""
    return datetime.datetime.strptime(s, '%d%b%y').date()


#-----------------------------------------------------------
def daterange(begin, end, delta = datetime.timedelta(1)):
    """Form a range of dates and iterate over them.

    Arguments:
    begin -- a date (or datetime) object; the beginning of the range.
    end   -- a date (or datetime) object; the end of the range.
    delta -- (optional) a timedelta object; how much to step each iteration.
             Default step is 1 day.

    Code origin: http://code.activestate.com/recipes/574441/
    (Python Recipe #574441, author: Michael Cornelius)
    """
    if not isinstance(delta, datetime.timedelta):
        delta = datetime.timedelta(delta)

    ZERO = datetime.timedelta(0)

    if begin < end:
        if delta <= ZERO:
            raise StopIteration
        test = end.__gt__
    else:
        if delta >= ZERO:
            raise StopIteration
        test = end.__lt__

    while test(begin):
        yield begin
        begin += delta

#-----------------------------------------------------------
def date_list_to_intervals(datelist):
    """Convert a list of dates to an equivalent set of intervals.

    Returns tuples (start, stop) where the intervals should be
    start <= date < stop.  Datelist should be datetime.date or
    datetime.datetime objects.

    """
    if len(datelist)==0:
        return []

    datelist.sort()
    intervals=[]
    cur_start=datelist[0]
    for d,nd in zip(datelist[:-1],datelist[1:]):
        if isinstance(d, datetime.datetime):
            d=d.date()
        if isinstance(nd, datetime.datetime):
            nd=nd.date()

        if (nd-d) > datetime.timedelta(1):
            intervals.append((cur_start, d+datetime.timedelta(1)))
            cur_start=nd
    # above loop will never include the last date to close
    # an interval.  Add it.
    intervals.append((cur_start, datelist[-1]+datetime.timedelta(1)))

    return intervals


#-----------------------------------------------------------
def _divide_date_interval(start, end, max_dt):
    """Helper function for divide_date_intervals.

    Divides a single interval into multiple, no longer than
    max_dt (datetime.timedelta).

    """
    intervals=[]
    while (end-start) > max_dt:
        intervals.append( (start, start+max_dt) )
        start=start+max_dt
    intervals.append((start, end))
    return intervals


#-----------------------------------------------------------
def divide_date_intervals(datelist, max_dt):
    """Divide a list of date intervals into chunks no longer than max_dt.

    Divides a list of intervals into shorter intervals.

    Arguments:

    datelist - list of (datetime,datetime) tuples representing start
    and end dates for each interval.  The convention is datelist[0] <=
    date < datelist[1].  This list will not be modified.

    max_dt - datetime.timedelta representing the maximum time change to
    permit in a single interval.

    Return:

    list of tuples in same format as datelist, representing the
    finer-grained intervals.

    """

    newlist=[]

    for interval in datelist:
        newlist.extend(_divide_date_interval(interval[0],interval[1],max_dt))
    return newlist


#-----------------------------------------------------------
def cmbprog_time_to_datetime(ref_datetime, time):
    """Convert a CMBPROG time member to a datetime.

    CMBPROG stores its 'time' members as seconds since a reference
    time.  The reference time is the time of the first procedure
    loaded.  To convert this to an absolute time, you need to find the
    date/time of that procedure and pass it in as ref_datetime (should
    be a datetime structure).  This routine goes through the necessary
    simple but fiddly arithmetic to convert the time offset into a
    real datetime.
    """
    t_0 = mktime(ref_datetime.timetuple()) # get seconds since 1970-01-01
    return datetime.datetime.fromtimestamp(t_0+time)


#-----------------------------------------------------------
def fermi_MET_to_datetime(met):
    """Convert MET to a datetime.

    Fermi uses "Mission-Elapsed-Time" (MET) as a measure of time
    in its data.  MET is defined as seconds since 2001-01-01 00:00:00 UTC.

    """
    return datetime.datetime(2001,1,1,0,0,0)+datetime.timedelta(seconds=met)

#-----------------------------------------------------------
def mjd_to_lst_day(mjd):
    # equation from http://www.vla.nrao.edu/astro/guides/dopset/
    return 6572.1572917 + 1.002737909350759 * mjd

def lst_day_to_mjd(lstd):
    return (lstd-6572.1572917)/1.002737909350759

#-----------------------------------------------------------
def get_decimal_year(d):
    def _decyear(dd):
        dt=mjd_to_datetime(get_mjd(dd))
        this_year=datetime.datetime(year=dt.year,month=1,day=1)
        next_year=datetime.datetime(year=dt.year+1,month=1,day=1)
        year_length=(next_year-this_year).total_seconds()
        return dt.year + (dt-this_year).total_seconds()/year_length

    try:
        return numpy.array([_decyear(dd) for dd in d])
    except TypeError:
        return _decyear(d)
