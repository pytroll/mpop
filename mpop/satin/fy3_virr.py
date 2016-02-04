#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015, 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Katerina.Melnik <kmelnik@scanex.ru>

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

"""A VIRR reader for FY3-B and maybe A....
"""
import numpy as np
import os
import logging
from datetime import datetime
from ConfigParser import ConfigParser
from mpop import CONFIG_PATH
import h5py
from pyspectral.blackbody import blackbody_wn_rad2temp as rad2temp

LOGGER = logging.getLogger('virr')


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

    LOGGER.debug("Calibrate = " + str(options["calibrate"]))
    LOGGER.info("Loading instrument '%s'", satscene.instrument_name)

    try:
        CASES[satscene.instrument_name](satscene, options)
    except KeyError:
        raise KeyError("Unknown instrument '%s'" % satscene.instrument_name)


def load_virr(satscene, options):
    """Read the VIRR hdf5 file"""

    if "filename" not in options:
        raise IOError("No 1km virr filename given, cannot load")

    values = {"orbit": satscene.orbit,
              "satname": satscene.satname,
              "instrument": satscene.instrument_name,
              "satellite": satscene.fullname
              }

    filename = \
        os.path.join(satscene.time_slot.strftime(options["dir"]) % values,
                     satscene.time_slot.strftime(
                         options["filename"])
                     % values)

    LOGGER.debug("Filename= %s", filename)

    datasets = ['EV_Emissive',
                'EV_RefSB']

    calibrate = options['calibrate']
    LOGGER.debug("Calibrate = " + str(calibrate))

    h5f = h5py.File(filename, 'r')

    # Get geolocation information
    lons = h5f['Longitude'][:]
    lats = h5f['Latitude'][:]
    # Mask out unrealistic values:
    mask = np.logical_or(lats > 90., lons > 90.)
    lons = np.ma.masked_array(lons, mask=mask)
    lats = np.ma.masked_array(lats, mask=mask)
    sunz = h5f['SolarZenith'][:]
    slope = h5f['SolarZenith'].attrs['Slope'][0]
    intercept = h5f['SolarZenith'].attrs['Intercept'][0]
    sunz = sunz * slope + intercept
    sunz = np.where(np.greater(sunz, 85.0), 85.0, sunz)

    # Get the calibration information
    # Emissive radiance coefficients:
    emis_offs = h5f['Emissive_Radiance_Offsets'][:]
    emis_scales = h5f['Emissive_Radiance_Scales'][:]

    # Central wave number (unit =  cm-1) for the three IR bands
    # It is ordered according to decreasing wave number (increasing wavelength):
    # 3.7 micron, 10.8 micron, 12 micron
    emiss_centroid_wn = h5f.attrs['Emmisive_Centroid_Wave_Number']

    # VIS/NIR calibration stuff:
    refsb_cal_coeff = h5f.attrs['RefSB_Cal_Coefficients']
    visnir_scales = refsb_cal_coeff[0::2]
    visnir_offs = refsb_cal_coeff[1::2]

    refsb_effective_wl = h5f.attrs['RefSB_Effective_Wavelength']

    # Read the band data:
    for dset in datasets:
        band_data = h5f[dset]
        valid_range = band_data.attrs['valid_range']
        LOGGER.debug("valid-range = " + str(valid_range))
        fillvalue = band_data.attrs['_FillValue']
        band_names = band_data.attrs['band_name'].split(',')
        slope = band_data.attrs['Slope']
        intercept = band_data.attrs['Intercept']
        units = band_data.attrs['units']
        long_name = band_data.attrs['long_name']

        LOGGER.debug('band names = ' + str(band_names))

        for (i, band) in enumerate(band_names):
            if band not in satscene.channels_to_load:
                continue

            LOGGER.debug("Reading channel %s, i=%d", band, i)
            data = band_data[i]

            bandmask = np.logical_or(np.less(data, valid_range[0]),
                                     np.greater(data, valid_range[1]))

            if calibrate:
                if dset in ['EV_Emissive']:
                    data = (np.array([emis_offs[:, i]]).transpose() +
                            data * np.array([emis_scales[:, i]]).transpose())
                    # Radiance to Tb conversion.
                    # Pyspectral wants SI units,
                    # but radiance data are in mW/m^2/str/cm^-1 and wavenumbers are in cm^-1
                    # Therefore multply wavenumber by 100 and radiances by
                    # 10^-5
                    data = rad2temp(emiss_centroid_wn[i] * 100., data * 1e-5)
                    LOGGER.debug("IR data calibrated")

                if dset in ['EV_RefSB']:
                    data = (visnir_offs[i] +
                            data * visnir_scales[i]) / np.cos(np.deg2rad(sunz))

            satscene[band] = np.ma.masked_array(data,
                                                mask=bandmask,
                                                copy=False)

    from pyresample import geometry
    satscene.area = geometry.SwathDefinition(lons=lons, lats=lats)

    h5f.close()


CASES = {
    "virr": load_virr,
}
