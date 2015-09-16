#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Abhay Devasthale and Martin Raspaud

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmach.it>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Loader for ascat, netcdf format.
   The driver works for netcdf format of ASCAT soil moisture swath data downloaded from 
   here: http://navigator.eumetsat.int/discovery/Start/DirectSearch/DetailResult.do?f%28r0%29=EO:EUM:DAT:METOP:SOMO12
   rename the CONFIG file mpop/mpop/etc/metop.ascat.cfg.template to metop.cfg to read the ASCAT data
"""
import numpy as np
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import os

from netCDF4 import Dataset


def load(satscene):
    """Load ascat data.
    """

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "number": satscene.number,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname,
              "time_slot": satscene.time_slot,
              "time": satscene.time_slot.strftime('%Y%m%d%H%M%S')
              }
    filename = os.path.join(
        conf.get("ascat-level2", "dir"),
        satscene.time_slot.strftime(conf.get("ascat-level2",
                                             "filename",
                                             raw=True)) % values)
    # Load data from netCDF file
    ds = Dataset(filename, 'r')
    for chn_name in satscene.channels_to_load:
        # Read variable corresponding to channel name
        data = np.ma.masked_array(
            ds.variables[chn_name][:], np.isnan(ds.variables[chn_name][:]))
        satscene[chn_name] = data
    lons = ds.variables['longitude'][:]
    lats = ds.variables['latitude'][:]

    # Set scene area as pyresample geometry object
    try:
        from pyresample import geometry
        satscene.area = geometry.SwathDefinition(lons=lons, lats=lats)
    except ImportError:
        # pyresample not available. Set lon and lats directly
        satscene.area = None
        satscene.lat = lats
        satscene.lon = lons
