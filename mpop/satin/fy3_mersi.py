#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c20671.ad.smhi.se>

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

"""A reader for the FY3 Mersi-1 
"""


import numpy as np
import os
import logging
from datetime import datetime
import glob
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import h5py
import pdb

LOGGER = logging.getLogger('mersi-1')


def load(satscene, *args, **kwargs):
    """Read data from file and load it into *satscene*.
    A possible *calibrate* keyword argument is passed to the AAPP reader. 
    Should be 0 for off (counts), 1 for default (brightness temperatures and
    reflectances), and 2 for radiances only.

    If *use_extern_calib* keyword argument is set True, use external
    calibration data.

    """
    del args

    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))
    options = {}
    for option, value in conf.items(satscene.instrument_name + "-level2",
                                    raw=True):
        options[option] = value

    if kwargs.get("filename") is not None:
        options["full_filename"] = kwargs["filename"]
    if kwargs.get("calibrate") is not None:
        options["calibrate"] = kwargs["calibrate"]
    else:
        options["calibrate"] = True

    LOGGER.info("Loading instrument '%s'", satscene.instrument_name)

    try:
        CASES[satscene.instrument_name](satscene, options)
    except KeyError:
        raise KeyError("Unknown instrument '%s'" % satscene.instrument_name)


def load_mersi(satscene, options):
    """Read the Mersi-1 hdf file"""

    if "filename_1000m" not in options:
        raise IOError("No 1km mersi-1 filename given, cannot load.")

    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname
              }

    filename_1000m = \
        os.path.join(satscene.time_slot.strftime(options["dir"]) % values,
                     satscene.time_slot.strftime(
                         options["filename_1000m"])
                     % values)

    LOGGER.debug("Filename= %s", filename_1000m)

    datasets = ['EV_250_Aggr.1KM_RefSB',
                'EV_250_Aggr.1KM_Emissive',
                'EV_1KM_RefSB']

    calibrate = options['calibrate']

    # Get the calibration information:
    h5f = h5py.File(filename_1000m)
    # The K0, K1 and K2 coefficients:
    vis_cal_coeff = h5f['Calibration']['VIS_Cal_Coeff'][:]
    # See also "Update of Calibration for Reflective Solar Bands of MERSI / FY-3C"
    # http://satellite.cma.gov.cn/PortalSite/Download/FY3C/CalibrationCoefficient/Update%20of%20Calibration%20for%20Reflective%20Solar%20Bands%20of%20MERSI_20140618.doc

    sv_dn_average = h5f['Calibration']['SV_DN_average'][:]
    # Expand array over all lines (10 lines per scan):
    sv_dn_average = np.repeat(sv_dn_average, 10, axis=1)

    date_orig = h5f.attrs['DN_Normalized_LUT_UpdateDate']
    dtobj_orig = datetime.strptime(date_orig, '%Y-%m-%d')
    obs_beg_date = h5f.attrs["Observing Beginning Date"]
    obs_beg_time = h5f.attrs["Observing Beginning Time"]
    dtobj_obs = datetime.strptime(
        obs_beg_date + obs_beg_time, '%Y-%m-%d%H:%M:%S.%f')
    h5f.close()

    # Get the days since 'launch' or since coefficients update:
    dsl = (dtobj_obs - dtobj_orig).days
    slopes = (vis_cal_coeff[:, 0] +
              vis_cal_coeff[:, 1] * dsl +
              vis_cal_coeff[:, 2] * dsl * dsl)
    # The slopes are available for band 1-4 and 6-20.
    # To keep consistency with the other cal-coefficients we add the IR band as
    # well, and set the slope to 1:
    slopes = np.concatenate((slopes[0:4], [1], slopes[4:]))

    mersi_band_index = 0
    with h5py.File(filename_1000m) as h5f:

        for dset in datasets:
            band_data = h5f['Data'][dset]
            valid_range = band_data.attrs['valid_range']
            LOGGER.debug("valid-range = " + str(valid_range))
            # FIXME! There seem to be useful data outside the valid range!
            valid_range = (0, 65535)
            fillvalue = band_data.attrs['FillValue']
            band_names = band_data.attrs['band_name'].split(',')
            slope = band_data.attrs['Slope']
            intercept = band_data.attrs['Intercept']

            LOGGER.debug('band names = ' + str(band_names))
            for (i, band) in enumerate(band_names):
                if band not in satscene.channels_to_load:
                    continue

                LOGGER.debug("Reading channel %s, i=%d", band, i)

                # Take care of the case when there is only one
                # single band (band 5: IR) in the dataset:
                if len(band_data.shape) == 2:
                    data = band_data
                else:
                    data = band_data[i]

                bandmask = np.logical_or(np.less(data, valid_range[0]),
                                         np.greater(data, valid_range[1]))

                if calibrate:
                    data = slopes[mersi_band_index] * (
                        data - np.array([sv_dn_average[mersi_band_index]]).transpose())

                satscene[band] = np.ma.masked_array(data,
                                                    mask=bandmask,
                                                    copy=False)

                satscene[band].info = {
                    'var_name': 'ch' + str(band),
                    'var_data': satscene[band].data,
                    'var_dim_names': ('x', 'y'),
                    '_FillValue': fillvalue,
                    'standard_name': '',
                    'short_name': band,
                    'scale_factor': slope,
                    'add_offset': intercept,
                }

                mersi_band_index = mersi_band_index + 1

    satscene.info = {
        'Antenna': 'None',
        'Receiver': 'Unknown',
        'Time': satscene.time_slot.strftime("%Y-%m-%d %H:%M:%S UTC"),
        'Area_Name': "swath",
        #'Projection': 'satproj',
        'Platform Name': satscene.fullname,
        'Service': '',
        #'Columns' : satscene.channels[0].shape[1],
        #'Lines' : satscene.channels[0].shape[0],
        'SampleX': 1.0,
        'SampleY': 1.0,
        'title': 'MERSI Level 1',
    }

    # Get geolocation information


CASES = {
    "mersi/1": load_mersi,
}
