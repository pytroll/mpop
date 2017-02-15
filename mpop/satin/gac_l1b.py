#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Abhay Devasthale and Martin Raspaud

# Author(s):

#   Abhay Devasthale <abhay.devasthale@smhi.se>
#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>
#   Sajid Pareeth <sajid.pareeth@fmch.it>

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
"""Read a gac file.
Reads L1b GAC data from KLM series of satellites (NOAA-15 and later) and does most of the computations.
Format specification can be found here:
http://www.ncdc.noaa.gov/oa/pod-guide/ncdc/docs/klm/html/c8/sec83142-1.htm

"""

import glob
import logging
import os
from ConfigParser import ConfigParser

import numpy as np

from mpop import CONFIG_PATH
from pygac.gac_klm import KLMReader
from pygac.gac_pod import PODReader

LOGGER = logging.getLogger(__name__)


def load(satscene, *args, **kwargs):
    """Read data from file and load it into *satscene*.
    A possible *calibrate* keyword argument is passed to the AAPP reader.
    Should be 0 for off (counts), 1 for default (brightness temperatures and
    reflectances), and 2 for radiances only.
    """
    del args

    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    options = {}
    for option, value in conf.items(satscene.instrument_name + "-level2",
                                    raw=True):
        options[option] = value

    if kwargs.get("filename") is not None:
        options["filename"] = kwargs["filename"]
        options["dir"] = None

    options["calibrate"] = kwargs.get("calibrate", True)

    LOGGER.info("Loading instrument '%s'" % satscene.instrument_name)

    try:
        CASES[satscene.instrument_name](satscene, options)
    except KeyError:
        raise KeyError("Unknown instrument '%s'" % satscene.instrument_name)


def load_avhrr(satscene, options):
    """Read avhrr data from file and load it into *satscene*.
    """
    if "filename" not in options:
        raise IOError("No filename given, cannot load.")

    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "number": satscene.number,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname}

    if options["dir"] is None:
        filename = options["filename"]
    else:
        filename = os.path.join(
            satscene.time_slot.strftime(options["dir"]) % values,
            satscene.time_slot.strftime(options["filename"]) % values)

        file_list = glob.glob(filename)

        if len(file_list) > 1:
            raise IOError("More than one l1b file matching!")
        elif len(file_list) == 0:
            raise IOError("No l1b file matching!: " + filename)

        filename = file_list[0]

    LOGGER.debug("Loading from " + filename)

    with open(filename) as fdes:
        data = fdes.read(3)
    if data in ["CMS", "NSS", "UKM", "DSS"]:
        reader = KLMReader
        chn_dict = AVHRR3_CHANNEL_NAMES
    else:
        reader = PODReader
        chn_dict = AVHRR_CHANNEL_NAMES

    chns = satscene.channels_to_load & set(chn_dict.keys())
    LOGGER.info("Loading channels " + str(sorted(list(chns))))

    if len(chns) == 0:
        return

    scene = reader()
    scene.read(filename)
    scene.get_lonlat()
    scene.adjust_clock_drift()
    channels = scene.get_calibrated_channels()

    # scene.navigate()
    try:
        from pyresample import geometry
    except ImportError as ex_:

        LOGGER.debug("Could not load pyresample: " + str(ex_))

        satscene.lat = scene.lats
        satscene.lon = scene.lons
    else:
        satscene.area = geometry.SwathDefinition(lons=scene.lons,
                                                 lats=scene.lats)
        area_name = ("swath_" + satscene.fullname + "_" +
                     str(satscene.time_slot) + "_" + str(scene.lats.shape))
        satscene.area.area_id = area_name
        satscene.area.name = "Satellite projection"
        satscene.area_id = area_name

    for chn in chns:
        data = channels[:, :, chn_dict[chn]]
        if np.ma.count(data) > 0:
            satscene[chn].data = np.ma.masked_invalid(data, copy=False)
            satscene[chn].area = satscene.area


AVHRR3_CHANNEL_NAMES = {"1": 0, "2": 1, "3A": 2, "3B": 3, "4": 4, "5": 5}
AVHRR_CHANNEL_NAMES = {"1": 0, "2": 1, "3": 2, "4": 3, "5": 4}

CASES = {"avhrr/1": load_avhrr, "avhrr/2": load_avhrr, "avhrr/3": load_avhrr, }
