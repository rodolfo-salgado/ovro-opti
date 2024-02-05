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
import pylab
import math
import time
import copy
import telescope_data as tel_data
import coord_utils as cu
import ipynb.fs.full.region_optimizer as ro
import ipynb.fs.full.py40m_astro as pastro