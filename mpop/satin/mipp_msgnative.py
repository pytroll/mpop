#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Northaholic <northaholic@icloud.com>

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

"""Interface to Eumetcast level 1.5 Nativeformat. Uses the MIPP reader.
"""

import ConfigParser
import os

from datetime import datetime
from mpop.plugin_base import Reader
from mpop import CONFIG_PATH

import logging
LOG = logging.getLogger(__name__)

from mipp.native.MSG import NativeImage as NatImg

try:
    # Work around for on demand import of pyresample. pyresample might depend
    # on scipy.spatial which memory leaks on multiple imports
    is_pyresample_loaded = False
    from pyresample import geometry
    from mpop.projector import get_area_def
    is_pyresample_loaded = True
except ImportError:
    LOG.warning("pyresample missing. Can only work in satellite projection")


class NativeReader(Reader):

    pformat = "mipp_native"

    def load(self, satscene, calibrate=1, filename=None, **kwargs):
        """Load the data"""

        conf = ConfigParser.ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
        options = {}
        for option, value in conf.items(satscene.instrument_name + "-level2"):
            options[option] = value

        for section in conf.sections():
            if(section.startswith(satscene.instrument_name) and
               not (section == "satellite") and
               not section[:-1].endswith("-level") and
               not section.endswith("-granules")):
                options[section] = conf.items(section)

        for option, value in conf.items("satellite"):
            options[option] = value

        # Set the time from the filename:
        # FIXME!
        satscene.time_slot = datetime(2013, 11, 9, 12, 0)

        # Build an area on the fly from the mipp metadata
        proj_params = options["proj4_params"].split(" ")
        proj_dict = {}
        for param in proj_params:
            key, val = param.split("=")
            proj_dict[key] = val

        # Instanciate from mipp:
        image = NatImg(str(satscene.satname) +
                       str(satscene.number), filename=filename)

        for chn in satscene.channels_to_load:
            chobj = getattr(image, chn.lower())
            satscene[chn] = chobj.data

            shape = chobj.data.shape

            # Get the unit using the _get_unitnames function!
            satscene[chn].info['units'] = 'K'  # FIXME!
            satscene[chn].info['satname'] = satscene.satname
            satscene[chn].info['satnumber'] = satscene.number
            satscene[chn].info['instrument_name'] = satscene.instrument_name
            satscene[chn].info['time'] = satscene.time_slot

            if is_pyresample_loaded:
                # Build area_def on-the-fly
                satscene[chn].area = geometry.AreaDefinition(
                    satscene.satname + satscene.instrument_name +
                    str(image.area_extent) +
                    str(shape),
                    "On-the-fly area",
                    proj_dict["proj"],
                    proj_dict,
                    shape[1],
                    shape[0],
                    image.area_extent)
            else:
                LOG.info("Could not build area, pyresample missing...")
