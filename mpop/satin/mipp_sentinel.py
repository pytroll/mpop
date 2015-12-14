#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Helge Pfeiffer <rhp@dmi.dk>
#   Lars Orum Rasmussen <ras@dmi.dk>

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

"""Sentinel-1 reader
"""

import ConfigParser
import os

from mipp.read_geotiff import read_geotiff
from mipp.xsar import S1A

from datetime import datetime
import numpy as np

from glob import glob

from pyresample import geometry
import mpop
from mpop import CONFIG_PATH
import logging

LOG = logging.getLogger(__name__)

from mpop.plugin_base import Reader

class SentinelGRDChannel(mpop.channel.GenericChannel):

    def __init__(self, name='unknown', resolution='unknown'):
        mpop.channel.GenericChannel.__init__(self)
        self._is_loaded = False
        self.name = name
        self.resolution = resolution
        self.data = None
        self.shape = None
        self._projectables = []
        self._projectables.append(name)

    def is_loaded(self):
        return self._is_loaded

    def set_loaded(self):
        self._is_loaded = not self._is_loaded

    def project(self, coverage):
        """Project what can be projected in the product.
        """
        import copy
        res = copy.copy(self)
        res.data = coverage.project_array(self.data)
        return res

class GeoTiffReader(Reader):

    pformat = "mipp_sentinel"

    def load(self, satscene, **kwargs):

        LOG.debug('channels to load: ' + str(satscene.channels_to_load))
        conf = ConfigParser.ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
        options = {}

        for option, value in conf.items(satscene.instrument_name + "-level2",
                                        raw=True):
            options[option] = value

        options["resolution"] = kwargs.get("resolution", 'unknown')
        values = {"orbit": satscene.orbit}

        path_template = datetime.strftime(satscene.time_slot, options["dir"]) % values

        dirlist = glob(path_template)
        
        if len(dirlist) != 1:
            raise IOError("Couldn't identify a unique measurments directory!, " +
                          "from path_template '%s'" % path_template)

        dirname = dirlist[0]
        if not os.path.exists(dirname):
            raise IOError('Directory ' + str(dirname) + ' does not exist')

        # Read meta data
        mda = S1A.read_metadata(dirname)
        channels_available = set(mda.channels.keys())

        ##filelist = glob(os.path.join(dirname, options["filename"]))
        ##if len(filelist) == 0:
        ##    LOG.warning('No files found!')

        channels_to_load = satscene.channels_to_load.intersection(channels_available)
        
        # Loading of channels
        LOG.debug('available channels to load: ' + str(channels_to_load))
        for channel_name in channels_to_load:
            channel_file = mda.channels[channel_name]
            LOG.debug("Load channel: '%s' %s" % (channel_name, str(channel_file)))
            lons, lats, data = self.load_channel(channel_file)
    
            chn = SentinelGRDChannel(channel_name, mda.pixel_spacing[0])
            chn.area = geometry.SwathDefinition(lons=lons, lats=lats)
            chn.data = np.ma.masked_array(data)
            chn.shape = data.shape
            chn.set_loaded()
            satscene[channel_name] = chn

        satscene.info['manifest'] = mda


    def load_channel(self, filename):
        """Load one sentinel channel file"""
        from geotiepoints.basic_interpolator import BasicSatelliteInterpolator

        params, data = read_geotiff(filename)
        tie_lons = params['tiepoints']['lons']
        tie_lats = params['tiepoints']['lats']
        tie_cols = params['tiepoints']['cols']
        tie_rows = params['tiepoints']['rows']

        interpolator = BasicSatelliteInterpolator(tie_cols, 
                                                  tie_rows, 
                                                  tie_lats, 
                                                  tie_lons)
        lats, lons = interpolator.interpolate()
        
        return lons, lats, data
