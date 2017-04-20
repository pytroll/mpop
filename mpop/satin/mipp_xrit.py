#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010, 2011, 2013, 2014, 2015.

# SMHI,
# Folkborgsvägen 1,
# Norrköping,
# Sweden

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Esben S. Nielsen <esn@dmi.dk>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This file is part of mpop.

# mpop is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# mpop is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# mpop.  If not, see <http://www.gnu.org/licenses/>.

"""Interface to Eumetcast level 1.5 HRIT/LRIT format. Uses the MIPP reader.
"""
import ConfigParser
import fnmatch
import logging
import os

from pyproj import Proj

from mipp import CalibrationError, ReaderError, xrit
from mpop import CONFIG_PATH
from mpop.plugin_base import Reader
from mpop.satin.helper_functions import area_def_names_to_extent

LOGGER = logging.getLogger(__name__)

try:
    # Work around for on demand import of pyresample. pyresample depends
    # on scipy.spatial which memory leaks on multiple imports
    IS_PYRESAMPLE_LOADED = False
    from pyresample import geometry
    from mpop.projector import get_area_def
    IS_PYRESAMPLE_LOADED = True
except ImportError:
    LOGGER.warning("pyresample missing. Can only work in satellite projection")


class XritReader(Reader):

    '''Class for reading XRIT data.
    '''
    pformat = "mipp_xrit"

    def load(self, *args, **kwargs):
        load(*args, **kwargs)


def load(satscene, calibrate=True, area_extent=None, area_def_names=None,
         **kwargs):
    """Read data from file and load it into *satscene*. The *calibrate*
    argument is passed to mipp (should be 0 for off, 1 for default, and 2 for
    radiances only).
    """

    conf = ConfigParser.RawConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    options = {}
    for option, value in conf.items(satscene.instrument_name + "-level2"):
        options[option] = value

    for section in conf.sections():
        if(section.startswith(satscene.instrument_name) and
           not (section == "satellite") and
           # not section[:-1].endswith("-level") and
           not section.endswith("-granules")):
            options[section] = dict(conf.items(section))

    filenames = kwargs.get('filename')

    CASES.get(satscene.instrument_name, load_generic)(satscene,
                                                      options,
                                                      calibrate,
                                                      area_extent,
                                                      area_def_names,
                                                      filenames)


def load_generic(satscene, options, calibrate=True, area_extent=None,
                 area_def_names=None, filenames=None):
    """Read imager data from file and load it into *satscene*.
    """

    os.environ["PPP_CONFIG_DIR"] = CONFIG_PATH

    LOGGER.debug("Channels to load from %s: %s" % (satscene.instrument_name,
                                                   satscene.channels_to_load))

    # Compulsory global attributes
    satscene.info["title"] = (satscene.satname.capitalize() + satscene.number +
                              " satellite, " +
                              satscene.instrument_name.capitalize() +
                              " instrument.")
    satscene.info["institution"] = "Original data disseminated by EumetCast."
    satscene.add_to_history("HRIT/LRIT data read by mipp/mpop.")
    satscene.info["references"] = "No reference."
    satscene.info["comments"] = "No comment."

    from_area = False

    if satscene.end_time is not None:
        time_slot = satscene.time_slot, satscene.end_time
    else:
        time_slot = satscene.time_slot

    if area_extent is None and satscene.area is not None:
        if not satscene.area_def:
            satscene.area = get_area_def(satscene.area_id)
        area_extent = satscene.area.area_extent
        from_area = True

    area_converted_to_extent = False

    for chn in satscene.channels_to_load:
        use_filenames = False
        # Sort out filenames
        if filenames is not None:
            for section in options.keys():
                if section.endswith('-level1'):
                    break
            try:
                pattern_pro = eval(options[section].get('filename_pro'))
            except TypeError:
                pattern_pro = None
            try:
                pattern_epi = eval(options[section].get('filename_epi'))
            except TypeError:
                pattern_epi = None
            pattern = eval(options[section].get('filename'))

            epilogue = None
            prologue = None
            image_files = []

            if pattern_epi is not None:
                glob_epi = satscene.time_slot.strftime(
                    pattern_epi) % ({'segment': "EPI".ljust(9, '_'),
                                     'channel': chn + '*'})
            else:
                glob_epi = 'eggs_and_spam'

            if pattern_pro is not None:
                glob_pro = satscene.time_slot.strftime(
                    pattern_pro) % ({'segment': "PRO".ljust(9, '_'),
                                     'channel': chn + '*'})
            else:
                glob_pro = 'eggs_and_spam'

            glob_img = satscene.time_slot.strftime(
                pattern) % ({'segment': "*", 'channel': chn + '*'})

            for filename in filenames:
                if fnmatch.fnmatch(os.path.basename(filename), glob_img):
                    image_files.append(filename)
                elif pattern_pro is not None and fnmatch.fnmatch(
                    os.path.basename(filename),
                        glob_pro):
                    prologue = filename
                elif pattern_epi is not None and fnmatch.fnmatch(
                    os.path.basename(filename),
                        glob_epi):
                    epilogue = filename
            if len(image_files) == 0 and prologue is None and epilogue is None:
                use_filenames = False
            else:
                use_filenames = True

        if from_area:
            try:
                if use_filenames:
                    metadata = xrit.sat.load_files(prologue,
                                                   image_files,
                                                   epilogue,
                                                   platform_name=satscene.fullname,
                                                   only_metadata=True)
                else:
                    metadata = xrit.sat.load(satscene.fullname,
                                             time_slot,
                                             chn,
                                             only_metadata=True)
                if(satscene.area_def.proj_dict["proj"] != "geos" or
                   float(satscene.area_def.proj_dict["lon_0"]) !=
                   metadata.sublon):
                    raise ValueError("Slicing area must be in "
                                     "geos projection, and lon_0 should match "
                                     "the satellite's position.")
            except ReaderError, err:
                # if channel can't be found, go on with next channel
                LOGGER.error(str(err))
                continue

        # Convert area definitions to maximal area_extent
        if not area_converted_to_extent and area_def_names is not None:
            try:
                if use_filenames:
                    metadata = xrit.sat.load_files(prologue,
                                                   image_files,
                                                   epilogue,
                                                   platform_name=satscene.fullname,
                                                   only_metadata=True)
                else:
                    metadata = xrit.sat.load(satscene.fullname,
                                             time_slot,
                                             chn,
                                             only_metadata=True)
            except ReaderError as err:
                LOGGER.warning(str(err))
                continue
            # if area_extent is given, assume it gives the maximum
            # extent of the satellite view
            if area_extent is not None:
                area_extent = area_def_names_to_extent(area_def_names,
                                                       metadata.proj4_params,
                                                       area_extent)
            # otherwise use the default value (MSG3 extent at
            # lon0=0.0), that is, do not pass default_extent=area_extent
            else:
                area_extent = area_def_names_to_extent(area_def_names,
                                                       metadata.proj4_params,
                                                       default_extent=None)

            if area_extent is None:
                LOGGER.info('Could not derive area_extent from area_def_names')

            area_converted_to_extent = True

        try:
            if use_filenames:
                image = xrit.sat.load_files(prologue,
                                            image_files,
                                            epilogue,
                                            platform_name=satscene.fullname,
                                            mask=True,
                                            calibrate=calibrate)
            else:
                image = xrit.sat.load(satscene.fullname,
                                      time_slot,
                                      chn,
                                      mask=True,
                                      calibrate=calibrate)
            if area_extent:
                metadata, data = image(area_extent)
            else:
                metadata, data = image()
        except CalibrationError:
            LOGGER.warning(
                "Loading non calibrated data since calibration failed.")
            if use_filenames:
                image = xrit.sat.load_files(prologue,
                                            image_files,
                                            epilogue,
                                            platform_name=satscene.fullname,
                                            mask=True,
                                            calibrate=False)
            else:
                image = xrit.sat.load(satscene.fullname,
                                      time_slot,
                                      chn,
                                      mask=True,
                                      calibrate=False)
            if area_extent:
                metadata, data = image(area_extent)
            else:
                metadata, data = image()

        except ReaderError as err:
            # if channel can't be found, go on with next channel
            LOGGER.warning(str(err))
            continue

        satscene[chn] = data

        satscene[chn].info['units'] = metadata.calibration_unit
        satscene[chn].info['sublon'] = metadata.sublon
        satscene[chn].info['satname'] = satscene.satname
        satscene[chn].info['satnumber'] = satscene.number
        satscene[chn].info['instrument_name'] = satscene.instrument_name
        satscene[chn].info['time'] = satscene.time_slot

        # Build an area on the fly from the mipp metadata
        proj_params = getattr(metadata, "proj4_params").split(" ")
        proj_dict = {}
        for param in proj_params:
            key, val = param.split("=")
            proj_dict[key] = val

        if IS_PYRESAMPLE_LOADED:
            # Build area_def on-the-fly
            satscene[chn].area = geometry.AreaDefinition(
                satscene.satname + satscene.instrument_name +
                str(metadata.area_extent) +
                str(data.shape),
                "On-the-fly area",
                proj_dict["proj"],
                proj_dict,
                data.shape[1],
                data.shape[0],
                metadata.area_extent)
        else:
            LOGGER.info("Could not build area, pyresample missing...")

CASES = {}
