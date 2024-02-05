"""
Sky optimizer
Given a list of regions with coordinates, delta_t_obs, and observability
ranges, returns a sorted list in such a way that all the regions are
observed taking a minimum time.
to do
-----
    - Add 3C286 for daily observations
        - Extend this to any other sources that need faster cadence
    - Improve optimization by changing the order in which observable sources
    are observed to reduce slewing time
    - Fill gaps in the observations automatically
Walter Max-Moerbeck, March 4, 2009.
"""
import pickle
import numpy
import random
# import pylab
# import math
import time
# import copy
import telescope_data as tel_data
import coord_utils as cu
# import ipynb.fs.full.region_optimizer as ro
import ipynb.fs.full.py40m_astro as pastro

def load_variables(savefile):
    """ Loads the sources and regions variables previously generated by
    region_optmizer.py
    """
    with open(savefile) as file:
        sources = pickle.load(file)
        regions = pickle.load(file)
    return sources, regions

def save_results(regions_sorted, regions_sorted_lst, savefile):
    """ Save the variables for use by the write schedule module
    regions_sorted:    Sorted list with the region numbers
    regions_sorted_lst: Sorted list with lsts of observation
    """
    with open(savefile, 'w') as file:
        pickle.dump(regions_sorted, file)
        pickle.dump(regions_sorted_lst, file)
    return

def get_region_by_number(regions, number):
    """ Return the region with the given Healpix number

    This number will not change as new sources are added to the schedule which
    makes it convenient to index them.

    Notes
    -----
    Since the region splitting, number is actually a string, because
    splitted regions get a string ID, e.g. '105a', '105b'.
    """
    # default value
    region = None
    # loop through the regions until finding it
    for reg in regions:
        if str(reg['number']) == str(number):
            region = reg
    return region

def sort_regions_by_lst(regions_order, regions_order_lst, lst_start):
    """ Sort regions by LST. Pop out the region that's there both at the beginning
    the end of the list (new sorting routine which allows any start LST time)
    Returns a region_order and lst vector.
    regions_order: original order of the regions
    regions_order_lst: original lst order of the regions

    NOTES:
    - Assumes that the wanted LST range is between 0 and 72h
    - The final LST list will not be optimized yet because it includes
    2 of the same regions and the other one is just popped out.
    Therefore need to simulate the order to get the final LST list.

    2012-04-20 / thovatta

    """
#    print "regions before sorting:", regions_order
    subtracted_lst=[]
    for lst in regions_order_lst:
        #account for schedules where initial lst start > wanted lst start
        if regions_order_lst[0] > (lst_start+0.5):
            lst -= 24
#        print 'lst before sorting', lst
        if lst < lst_start:
            lst += int(regions_order_lst[-1]/24)*24+24
#        print 'lst after sorting', lst
        subtracted_lst.append(lst)
    # Sort the array
    regions_sorted = list(numpy.array(regions_order)[numpy.array(subtracted_lst).argsort()])
    # sorted_lst = list(numpy.array(subtracted_lst)[numpy.array(subtracted_lst).argsort()])
    # print "regions after sorting:", regions_sorted, "at lst", sorted_lst
    return regions_sorted

def position_last_source_on_region(region, sources, lst):
    """ Get the ZA/AZ for the last source on a region at a given lst

    Returns za, az in degrees
    """
    # get name for last source in region
    last_source = region['sources'][region['order'][-1]]
    # get ra/dec for the source
    ra, dec = sources[last_source]['ra'], sources[last_source]['dec']
    # calculate za/az coordinates for lst
    za, az = cu.radec_zaaz(ra, dec, lst)
    return za, az

def calculate_region_obstime(region, sources, lst_start, az_t):
    """ Calculates the time taken to observe a region by taking the az wrap into
    account. Modification of the path_obstime routine in region_optimizer.py
    inputs:
    region = current region to be observed
    sources = list of sources
    lst_start = start lst of the region
    az_t = current az position of the telescope
    2012-04-23 / thovatta
    """
    # current lst
    lst = lst_start
    # add slew and observing time for different sources
    for i in range(len(region['order']) - 1):
        # get name of the sources
        a = region['sources'][region['order'][i]]
        b = region['sources'][region['order'][i+1]]
#        print 'a =', a
#        print 'region = ', region
#        print 'sources[a]=', sources[a]
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # time to observe first source and move on next one
        # get az/za for the two sources
        # Add observing time from source type, pointing or not
        lst += tel_data.obs_time(sources[a]['pointing'])
        # First source at end of the observation
        za_a, az_a = cu.radec_zaaz(sources[a]['ra'], sources[a]['dec'], lst)
        az_a = tel_data.move_in_azimuth(az_t, az_a)
        # Estimate the slew time using the position of second source at
        # end of observation first source taking the az wrap into account
        za_b, az_b = cu.radec_zaaz(sources[b]['ra'], sources[b]['dec'], lst)
        az_b = tel_data.move_in_azimuth(az_a, az_b)
        # add slew time
        lst += tel_data.slew_time(za_a, az_a, za_b, az_b)
        # store the current az of the telescope to az b
        az_t = az_b
    # Add observation time for last source
    lst += tel_data.obs_time(sources[region['sources'][region['order'][-1]]]['pointing'])
    # return the path length
    return lst - lst_start

def report_total_time(report):
    """ Take and observation simulation report and return total obseving time
    in hours
    wmax, March 23, 2012
    """
    # Times
    t_obs = [row[3] for row in report]
    t_slew = [row[4] for row in report]
    t_wait = [row[5] for row in report]
    # Total
    t_total = sum(t_obs) + sum(t_wait) + sum(t_slew)
    return t_total

def report_obs_time(report):
    """ Take an observation simulation report and return total obs time in hours
    thovatta / 20120419
    """
    # Times
    t_obs = [row[3] for row in report]
    return sum(t_obs)

def report_slew_time(report):
    """ Take an observation simulation report and return total slew time in hours
    thovatta / 20120419
    """
    # Times
    t_slew = [row[4] for row in report]
    return sum(t_slew)

def report_wait_time(report):
    """ Take an observation simulation report and return total wait time in hours
    thovatta / 20120419
    """
    # Times
    t_wait = [row[5] for row in report]
    return sum(t_wait)

def report_stats(report, printing=True):
    """ Generate a series of statistics using the report from
    simulate_regions_observation.
    """
    # coordinates
    za = [row[1] for row in report]
    az = [row[2] for row in report]
    # times
    t_obs = [row[3] for row in report]
    t_slew = [row[4] for row in report]
    t_wait = [row[5] for row in report]
    lst = [row[6] for row in report]
    if printing:
        # print to standard output a report with the observations
        print ('')
        print ('--------------------------------------------------')
        print ('Schedule summary')
        print ('number lst za   az   t_obs  t_slew  t_wait')
        print ('        h  deg  deg  h      h       h')
        for region_line in report:
            print ('%s\t%.5f\t%.2f\t%.2f\t%.5f\t%.5f\t%.5f' % (region_line[0],
                                                              region_line[6],
                                                              region_line[1],
                                                              region_line[2],
                                                              region_line[3],
                                                              region_line[4],
                                                              region_line[5]))
        # print basic stats
        print ('')
        print ('--------------------------------------------------')
        print ('total time for cycle = ',\
            sum(t_obs) + sum(t_wait) + sum(t_slew), ' hours')
        print ('t_obs = %.1f, t_slew = %.1f, t_wait = %.1f' % (sum(t_obs),
                                                          sum(t_slew),
                                                          sum(t_wait)))
    return za, az, t_obs, t_slew, t_wait, lst

def simulate_regions_observation(regions,
                                 regions_order,
                                 lst_start,
                                 sources,
                                 wait=False,
                                 za_t=0.0,
                                 az_t=180.0):
    """ Simulate a given region order

    wait=True    Wait for regions to be observable before moving to that
                 position
    wait=False   Just move to that position, even if not observable

    NOTE:
    Unlike the case for single region simulation. In this case the
    observation times can be quite long (~1 hour), so the positions before and
    after observation can be different.

    Between observations it is assumed that telescope is parked at the last
    ZAAZ

    Added that outputs also the lst time of a region into the report
    2012-04-21 / thovatta
    """
    # initialize report
    report = []
    # initialize time
    lst = lst_start % 24.0
    # loop through all the regions
    for i in range(len(regions_order)):
        # get current region
        curr_region = get_region_by_number(regions, regions_order[i])
        # Check that region is observable or wait==False
        # if not observable add t_wait
        if cu.check_observability(curr_region['obs_range'], lst) or wait is False:
            t_wait = 0.0
        elif wait:
            t_wait =\
                cu.wait_time_for_observability(curr_region['obs_range'], lst)
        # add wait time
        lst += t_wait
        # Get position at begining observation
        za_c, az_c = cu.radec_zaaz(curr_region['ra'], curr_region['dec'], lst)
        #----------------------------------------
        # New code for azimuth wrap
        # Considering current position of telescope convert geometric to
        # telescope coordinates with azimuth wrap incorporated
        az_c = tel_data.move_in_azimuth(az_t, az_c)
        #----------------------------------------
        # Get slew time from previous telescope position and add it to lst
        t_slew = tel_data.slew_time(za_t, az_t, za_c, az_c)
        lst += t_slew
        # copy the lst time as the observing lst
        obs_lst = lst
        # Add observation time
        t_obs = curr_region['obstime']
        lst += t_obs
        # update telescope position at end of observation
        # this is position for last source at the end
        za_ls, az_ls = position_last_source_on_region(curr_region, sources, lst)
        #----------------------------------------
        # New code for azimuth wrap
        # Considering current position of telescope convert geometric to
        # telescope coordinates with azimuth wrap incorporated
        az_ls = tel_data.move_in_azimuth(az_t, az_ls)
        #----------------------------------------
        za_t, az_t = za_ls, az_ls
        # add report line
        report.append([curr_region['number'],
                       za_c, az_c, t_obs, t_slew, t_wait, obs_lst])
    return report

def simulate_regions_final(regions,
                           regions_order,
                           lst_start,
                           sources,
                           wait=False,
                           za_t=0.0,
                           az_t=180.0):
    """ Simulate a given region order
    wait=True    Wait for regions to be observable before moving to that
    wait=False   Just move to that position, even if not observable
    NOTES:
    Unlike the case for single region simulation. In this case the
    observation times can be quite long (~1 hour), so the positions before and
    after observation can be different.
    Between observations it is assumed that telescope is parked at the last
    ZAAZ
    Added that outputs also the lst time of a region into the report

    2012-04-21 / thovatta

    This routine also calculates the real observing time of a region because it depends on the
    hour angle when the region is observed (some regions have az wrap in them)
    2012-04-22 / thovatta
    """
    # initialize report
    report = []
    # initialize time
    lst = lst_start #% 24.0
    # loop through all the regions
    for i in range(len(regions_order)):
        # get current region
        curr_region = get_region_by_number(regions, regions_order[i])
        # Check that region is observable or wait==False
        # if not observable add t_wait
        if cu.check_observability(curr_region['obs_range'], lst) or wait is False:
            t_wait = 0.0
        elif wait:
            t_wait =\
                cu.wait_time_for_observability(curr_region['obs_range'], lst)
        # add wait time
        lst += t_wait
        # Change compared to other simulation function: loop through the sources
        # in the region to get the correct slew and observing times for the region
        # First get the slew time to the region
        za_c, az_c = cu.radec_zaaz(curr_region['ra'], curr_region['dec'], lst)
        #----------------------------------------
        # New code for azimuth wrap
        # Considering current position of telescope convert geometric to
        # telescope coordinates with azimuth wrap incorporated
        az_c = tel_data.move_in_azimuth(az_t, az_c)
        #----------------------------------------
        # Get slew time from previous telescope position and add it to lst
        t_slew = tel_data.slew_time(za_t, az_t, za_c, az_c)
        lst += t_slew
        # copy the lst time as the observing lst
        obs_lst = lst
        # Get the observing time of that region
        t_obs = calculate_region_obstime(curr_region, sources, lst, az_c)
        #t_obs original from regions to compare
        print(f'Region: {curr_region['number']}', \
            'Obs time: {curr_region[obs_time]} (saved)',
            '{t_obs} (real)', sep='\t')
        lst += t_obs
        # update telescope position at end of observation
        # this is position for last source at the end
        za_ls, az_ls = position_last_source_on_region(curr_region, sources, lst)
        # telescope coordinates with azimuth wrap incorporated
        az_ls = tel_data.move_in_azimuth(az_t, az_ls)
        za_t, az_t = za_ls, az_ls
        # add report line
        report.append([curr_region['number'],
                       za_c, az_c, t_obs, t_slew, t_wait, obs_lst])
    return report

def order_regions_slew_time(regions,
                                sources,
                                lst_start=0,
                                lst_obs_win=12./60.,
                                za_t=0.0,
                                az_t=180.0,
                                delta_lst=5./60.,
                                daily_regions_list=\
                              [{'number': '136', #3C286 custom observing times
                                'obs_range': [[10.0, 11.7],[15.1,17.2]]},
                               {'number': '137', #DR21
                                'obs_range': [[16.25,18.8],[22.5,24.0]]},
                               {'number': '138', #3C48
                                'obs_range': [[5.3, 6.5],[20.5,22.7]]},
                               {'number': '142', #3C147
                                'obs_range': [[2.0, 4.0],[7.25,9.25]]},
                               {'number': '139', #NGC7027
                                'obs_range': [[16.75,19.3],[22.8,24.0]]},
                               {'number': '140', #3C295
                                'obs_range': [[10.5, 12.75],[15.5,18.75]]},
                               {'number': '141', #3C161
                                'obs_range': [[7.5, 9.25]]}],
                                sun_jd=0):
    """ Order regions minimizing the slew time only and starting at each of the 135 regions
        messy code which orders the sources based on smalles slewing time and returns the mimimun found.
        Daily calibration regions are accounted for in a not-so-clever way.
        This routine also accounts for the sun_angle by eliminating regions
        within 10 degrees from the sun on a give MJD. Returns the mimimun order, minimum time and minimum
        order_lst and the lst_start of the schedules
        2012-04-20 / thovatta
        #Modified to have all polarization calibrators to be observed once per day.
    """
    # First deal with daily regions and add observation counter to daily regions
    for reg in daily_regions_list:
        reg['n_obs'] = 0
        reg['lst'] = -30
    # daily_regions are eliminated from list, they are treated specially
    daily_regions_numbers = [reg['number'] for reg in daily_regions_list]
    regions_sun = [reg for reg in regions if\
                        reg['number'] not in daily_regions_numbers]
    #Check if the region is within 10 degrees from the sun in the middle of the schedule range
    #if a sun_jd is given
    regions_zero_lst = []
    if sun_jd != 0:
        sun_ra, sun_dec = cu.sun_pos(sun_jd)
        print ('sun ra dec=', sun_ra, sun_dec)
        for reg in regions_sun:
            sep=pastro.sky_sep(pastro.ephem.degrees(reg['ra']*numpy.pi/180.0),pastro.ephem.degrees(reg['dec']*numpy.pi/180.0),pastro.ephem.degrees(sun_ra*numpy.pi/180.0), pastro.ephem.degrees(sun_dec*numpy.pi/180.0)) * 180.0 / numpy.pi
            print ('region', reg['number'], 'sun_angle', sep, reg['ra'], reg['dec'])
            if sep > 10:
                regions_zero_lst.append(reg['number'])
    else:
        regions_zero_lst.append(reg['number'] for reg in regions_sun)
    print ('regions after sun angle', regions_zero_lst)
    print ('lst_start =', lst_start)
    initial_start_lst = lst_start #needed for checking if a or b schedules are being generated
    #check which regions are observable at LST 0 in order to select the one which has the best total time from those
    min_total_time = 10000
    min_order = []
    min_order_lst = []
    min_lst_start = 0
    lst_start=0
    while lst_start < 24:
        regions_order = []
        for region_number in regions_zero_lst:
        # get the region corresponding to region_number
            region = get_region_by_number(regions, region_number)
            if region is None:
                continue
            # find if region is observable at this lst and add it and the slew time to lists
            if cu.check_observability(region['obs_range'], lst_start) and\
                        cu.check_observability(region['obs_range'],
                                               lst_start + lst_obs_win + region['obstime']):
                regions_order.append(region['number'])
        #do the same for daily regions and see if they are observable at the start_lst
        for region_number in daily_regions_numbers:
        # get the region corresponding to region_number
            region = get_region_by_number(regions, region_number)
            if region is None:
                continue
            # find if region is observable at this lst and add it and the slew time to lists
            if cu.check_observability(region['obs_range'], lst_start) and\
                        cu.check_observability(region['obs_range'],
                                               lst_start + lst_obs_win + region['obstime']):
                regions_order.append(region['number'])
        print('regions observable at lst =', lst_start, 'are:', regions_order)
        # by starting at each region, go through other regions selecting always the shortest slew time
        # loop through the regions
        for region_number in regions_order:
            regions_added = []
            regions_added_lst = []
            non_grouped_regions = regions_zero_lst[:]
            # added = 0 #checking if time for 3C286 full track has been added
            # get the region corresponding to the region_number and define the lst to be the lst when this region becomes observable
            start_region = get_region_by_number(regions, region_number)
            #lst = start_region['obs_range'][0][0]
            lst = lst_start
            # simulate observation to check feasibility
            # region position at begining of observation
            za_r, az_r =\
                  cu.radec_zaaz(start_region['ra'], start_region['dec'], lst)
            # Consider azimuth wrap
            az_r = tel_data.move_in_azimuth(az_t, az_r)
            # slew time from previous telescope position
            t_slew_region = tel_data.slew_time(za_t, az_t, za_r, az_r)
            # observation time
            #t_obs_region = start_region['obstime']
            t_obs_region = calculate_region_obstime(start_region, sources, lst+t_slew_region, az_r)
            #update lst
            lst = lst + t_slew_region
            #lst_start = lst - t_slew_region
            #if region_number == 45:
            # print 'start lst =', lst_start, 'lst=', lst, 'slew = ', t_slew_region, 'range =', start_region['obs_range'], 't_obs = ', t_obs_region
            #the first region must be observable so add into regions_added and take away from regions left
            #check if this is cal region
            if region_number in daily_regions_numbers:
                regions_added.append(start_region['number'])
                regions_added_lst.append(lst)
                for reg in daily_regions_list:
                    if reg['number'] == region_number:
                        reg['n_obs'] += 1
                        reg['lst'] = lst
                lst += t_obs_region
            else:
                #region_pop_index =\
                #    numpy.where(numpy.array(non_grouped_regions) == region_number)[0] # old code before region splitting
                region_pop_index = non_grouped_regions.index(region_number)
                regions_added.append(start_region['number'])
                regions_added_lst.append(lst)
                non_grouped_regions.pop(region_pop_index)
                lst += t_obs_region
            # update telescope position at end of observation
             # this is position for last source at the end
            za_ls, az_ls =\
                    position_last_source_on_region(start_region,
                                                   sources,
                                                   lst)
            # Considering current position of telescope convert geometric to
            # telescope coordinates with azimuth wrap incorporated
            az_ls = tel_data.move_in_azimuth(az_t, az_ls)
            za_t, az_t = za_ls, az_ls
            # go through the remaining regions by adding the smallest slew time one as the next region
            while non_grouped_regions:
                # #Code that was used to see if 3C286 long track can somehow be automatically added. Could not get to work.
                # #check if LST is close to 8:30 LST when 3C286 should be started, only if start_lst != 0
                # print("initial start lst=", initial_start_lst, "added=", added, "lst=",lst)
                # if initial_start_lst != 0 and added == 0 and ((lst > 8 and lst < 9) or (lst > 32 and lst < 33)):
                #     print("adding time for 3C286 at lst", lst)
                #     lst += 10 #add the amount of time required to observe 3C286
                #     added = 1
                #     print("new lst=", lst)
                observable_regions = []
                observable_regions_slew = []
                # loop through the regions
                for region_index in non_grouped_regions:
                        # get the region corresponding to region_number
                    region = get_region_by_number(regions, region_index)
                    t_obs_region = calculate_region_obstime(region, sources, lst, az_t)
                        # find if region is observable at this lst and add it and the slew time to lists
                    if cu.check_observability(region['obs_range'], lst) and\
                                cu.check_observability(region['obs_range'],
                                                       lst + lst_obs_win + t_obs_region):
                            # region position at beginning observation
                        za_r, az_r =\
                                    cu.radec_zaaz(region['ra'], region['dec'], lst)
                            # consider Az wrap
                        az_r = tel_data.move_in_azimuth(az_t, az_r)
                            # slew time from previous telescope position
                        t_slew_region = tel_data.slew_time(za_t, az_t, za_r, az_r)
                        observable_regions_slew.append(t_slew_region)
                        observable_regions.append(region_index)
                print("observable regions at lst", lst, "are", observable_regions)
                # if no region is observable, advance time by a bit and goes to next cycle
                if len(observable_regions) == 0:
                    lst += delta_lst
                    continue
                # sort regions by slew time and get the region with shortest slew time
                index = numpy.array(observable_regions_slew).argsort()
                observable_regions = list(numpy.array(observable_regions)[index])
                region_index = observable_regions[0]
                region = get_region_by_number(regions, region_index)
                #region_pop_index =\
                #          numpy.where(numpy.array(non_grouped_regions) == region_index)[0] # old code before region splitting
                region_pop_index = non_grouped_regions.index(region_index)
                 # check that region is observable at the end of observation
                za_r, az_r =\
                          cu.radec_zaaz(region['ra'], region['dec'], lst)
                 # Consider azimuth wrap
                az_r = tel_data.move_in_azimuth(az_t, az_r)
                # slew time from previous telescope position
                t_slew_region = tel_data.slew_time(za_t, az_t, za_r, az_r)
                # observation time
                # t_obs_region = region['obstime']
                t_obs_region = calculate_region_obstime(region, sources, lst+t_slew_region, az_r)
                if cu.check_observability(region['obs_range'],
                                      lst + t_slew_region + t_obs_region):
                    # observe the source
                    regions_added.append(region['number'])
                    regions_added_lst.append(lst)
                    lst += t_obs_region + t_slew_region
                    non_grouped_regions.pop(region_pop_index)
                    # print("observing region ", region['number'])
                       # update telescope position at end of observation
                       # this is position for last source at the end
                    za_ls, az_ls =\
                         position_last_source_on_region(region,
                                                        sources,
                                                        lst)
                        # Considering current position of telescope convert geometric to
                        # telescope coordinates with azimuth wrap incorporated
                    az_ls = tel_data.move_in_azimuth(az_t, az_ls)
                    za_t, az_t = za_ls, az_ls
                        #check if a daily cal observation could be inserted here. Check that no cal region has been observed within the last
                        # 12 hours and that there are less than 4 added calibrators
                    for reg in daily_regions_list:
                        # get region information
                        reg_day = get_region_by_number(regions, reg['number'])
                        # FIX start:
                        if reg_day is None:
                            continue
                        # FIX end
                        if cu.check_observability(reg['obs_range'], lst) and\
                                     cu.check_observability(reg['obs_range'],lst + lst_obs_win + reg_day['obstime']) and\
                                     lst - reg['lst'] > 20 and\
                                     reg['n_obs'] <= lst / 24:
                            region = reg_day
                                  # check that region is observable at the end of observation
                            za_r, az_r =\
                                      cu.radec_zaaz(region['ra'], region['dec'], lst)
                                  # Consider azimuth wrap
                            az_r = tel_data.move_in_azimuth(az_t, az_r)
                                  # slew time from previous telescope position
                            t_slew_region = tel_data.slew_time(za_t, az_t, za_r, az_r)
                                  # observation time
                            t_obs_region = region['obstime']
                                  # check if the region is observable
                            if cu.check_observability(region['obs_range'],
                                                       lst + t_slew_region + t_obs_region):
                                regions_added.append(region['number'])
                                regions_added_lst.append(lst)
                                reg['n_obs'] += 1
                                reg['lst'] = lst
                                lst += t_obs_region + t_slew_region
                                       # update telescope position at end of observation                                                                                                                    # this is position for last source at the end
                                za_ls, az_ls =\
                                     position_last_source_on_region(region,
                                                                    sources,
                                                                    lst)
                                       # Considering current position of telescope convert geometric to                                                                                                     # telescope coordinates with azimuth wrap incorporated
                                az_ls = tel_data.move_in_azimuth(az_t, az_ls)
                                za_t, az_t = za_ls, az_ls
                                # print("observed calibrator region ", region['number'])
                                # break #break out of the loop not to include two regions in a row
                else:
                    lst += delta_lst
            # region_number=1 #test purposes!
            # # Add the first region into the list so that the schedule makes a full circle and lst start can be shifted
            # regions_added.append(region_number)
            # regions_added_lst.append(lst)
            sorted_regions = sort_regions_by_lst(regions_added, regions_added_lst, initial_start_lst)
            report =simulate_regions_final(regions,
                                            sorted_regions,
                                            initial_start_lst,
                                            sources,
                                            wait=True,
                                            za_t=za_t,
                                            az_t=az_t)
            t_total_i = report_total_time(report)
            # obs = report_obs_time(report)
            slew = report_slew_time(report)
            # wait = report_wait_time(report)
            za, az, t_obs, t_slew, t_wait, sorted_regions_lst = report_stats(report,printing=False)
            total_time = lst - lst_start
            print ('starting with region ', region_number, 'at lst', lst_start, 'total time = ',
                   total_time, 'total simulated time =', t_total_i,  'slew time', slew)
              #re-initialize daily regions
            for reg in daily_regions_list:
                  #print 'daily observations region', reg['number'], 'nobs and lst =', reg['n_obs'], reg['lst']
                reg['n_obs'] = 0
                reg['lst'] = -30
             #check if the total time is less than the mimimum time and if so, copy these into the min time lists
            if (t_total_i < min_total_time):
                min_total_time = t_total_i
                min_order = sorted_regions[:]
                min_order_lst = sorted_regions_lst[:]
                min_lst_start = initial_start_lst
        lst_start += 1
    return min_order, min_total_time, min_order_lst, min_lst_start

def modify_path_keep_cal(regions,
                regions_order,
                lst_start,
                sources,
                wait=False,
                za_t=0.0,
                az_t=180.0,
                cal_list = ['136', '137', '138', '139', '140', '141', '142']):
    """ Similar to modify_path but keeps the calibration regions intact. Needs a list of regions not to swap
        as an input.
        thovatta / 20120420
    """
    # Make a copy of regions order to modify without touching original
    regions_order_opt = regions_order[:]
    # Simulate orginal path and get total lenght
    report = simulate_regions_observation(regions,
                                          regions_order,
                                          lst_start,
                                          sources,
                                          wait=True,
                                          za_t=za_t,
                                          az_t=az_t)
    t_total_i = report_total_time(report)
    # Make a copy of current best solution
    regions_order_mod = regions_order_opt[:]
    cal_index = []
    #go through the regions and cal_list to find the corresponding indices
    for i in range(len(regions_order)):
        for cal in cal_list:
            if regions_order[i] == cal:
                cal_index.append(i)
                #print 'cal region', cal, 'at index', i
    # Swap the order of two regions leaving first region fixed and make sure that the swap numbers are not any
    # of the cal sources or the first and last source in the list which are set by the LST start time
    i_swap = [0,0]
    cal_check = 1
    cal_check_source = 0
    while cal_check:
        i_swap = random.sample(range(1, len(regions_order_mod)), 2)
        for i in cal_index:
            if i_swap[0] == i or i_swap[1] == i:
                cal_check_source = 1
        if i_swap[0] == 0 or i_swap[1] == 0 or\
                    i_swap[0] == len(regions_order_mod) or i_swap[1] == len(regions_order_mod) or\
                    i_swap[0] == i_swap[1] or\
                    cal_check_source == 1:
            cal_check = 1
            cal_check_source = 0
        else:
            cal_check = 0
    #print 'iswap = ', i_swap
    regions_1 = regions_order_mod[i_swap[0]]
    regions_order_mod[i_swap[0]] = regions_order_mod[i_swap[1]]
    regions_order_mod[i_swap[1]] = regions_1
    # Simulate and evaluate time
    report  = simulate_regions_observation(regions,
                                           regions_order_mod,
                                           lst_start,
                                           sources,
                                           wait=True,
                                           za_t=za_t,
                                           az_t=az_t)
    t_total = report_total_time(report)
    #print 'iswap vals = ', i_swap
    delta = t_total - t_total_i
    return regions_order_mod, t_total, delta

def genetic_algorithm_sky(regions, sources, order_opt, lst_start, tam_poblacion, prob_mutacion, num_generaciones):
    """Esta función está implementada para aplicar algoritmo genético,
    buscando minimizar t_wait, t_slew y t_obs.

    Antonia Bravo Rojo, Dic 20, 2023.
    """
    t0 = time.time()
    # report = simulate_regions_final(regions, order_opt, lst_start, sources, wait=True)   #entrega [region_number,za_c, az_c, t_obs, t_slew, t_wait, obs_lst]
    # minimum_obstime = report_obs_time(report)
    # Creating a random starting population
    unique_order_opt = list(dict.fromkeys(order_opt))
    # N = len(order_opt)
    M = len(unique_order_opt) #136 porque no considera las regiones 110 y 125.
    poblacion = []
    for _ in range(tam_poblacion):
        updated_regions, updated_time, delta =\
                                        modify_path_keep_cal(regions,
                                                       unique_order_opt,
                                                       lst_start,
                                                       sources,
                                                       wait=True,
                                                       za_t=0.0,
                                                       az_t=180.0)
        poblacion.append(updated_regions)
    for generacion in range(num_generaciones):
        for i in range(len(poblacion)):
            new_order = poblacion[i]
            new_report = simulate_regions_final(regions, new_order, lst_start, sources, wait=True)
            new_t_obs = report_obs_time(new_report)
            new_t_wait = report_wait_time(new_report)
            new_t_slew = report_slew_time(new_report)
            #if obstime is reduced we have new optimum and others new conditions
            #if new_t_obs < minimum_obstime:
            #if new_t_obs < minimum_obstime and (new_t_slew + new_t_wait)<15:
            if (new_t_obs + new_t_slew + new_t_wait)<88:
                # minimum_obstime = new_t_obs
                order_opt = new_order
        nueva_poblacion = order_opt[:]
        while len(nueva_poblacion) < tam_poblacion:
            padre1, padre2 = random.sample(poblacion[:int(0.5*tam_poblacion)], 2)
            punto_cruce = random.randint(0, M)
            hijo = padre1[:punto_cruce] + [x for x in padre2 if x not in padre1[:punto_cruce]]
            if random.random() < prob_mutacion:
                idx1, idx2 = random.sample(range(0,M), 2)
                hijo[idx1], hijo[idx2] = hijo[idx2], hijo[idx1]
            nueva_poblacion.append(hijo)
            nuevo_orden = nueva_poblacion[-1]
            nuevo_report = simulate_regions_final(regions, nuevo_orden, lst_start, sources, wait=True)
            nuevo_t_obs = report_wait_time(nuevo_report)
            nuevo_t_wait = report_wait_time(nuevo_report)
            nuevo_t_slew = report_slew_time(nuevo_report)
            #if nuevo_t_obs < minimum_obstime:
            #if nuevo_t_obs < minimum_obstime and (nuevo_t_slew + nuevo_t_wait)<15:
            if (nuevo_t_obs + nuevo_t_slew + nuevo_t_wait)<88:
                # minimum_obstime = nuevo_t_obs
                order_opt = nuevo_orden
    #print( 'Final length = ', minimum_obstime)
    print( 'optimum travel order =  ', order_opt)
    print( 'estamos trabajando con n° de fuentes igual a ', len(order_opt))
    print( '  \n')
    final_report = simulate_regions_final(regions, order_opt, lst_start, sources, wait=True)
    t_obs = report_obs_time(final_report)
    t_wait = report_wait_time(final_report)
    t_slew = report_slew_time(final_report)
    print( 't_obs = ', t_obs)
    print( 't_slew = ', t_slew)
    print( 't_wait = ', t_wait)
    print( '  \n')
    tf = time.time()
    print('time taken:', (tf-t0)/60, 'min')
    return order_opt