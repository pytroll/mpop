#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2015, 2016 Adam.Dybbroe

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

"""A reader for the OSISAF SST netCDF format
"""

import os.path
from ConfigParser import ConfigParser

import glob

import mpop.channel
from mpop import CONFIG_PATH
from mpop.plugin_base import Reader

from netCDF4 import Dataset
import numpy as np
import time
from datetime import datetime

import logging
LOG = logging.getLogger(__name__)


class InfoObject(object):

    """Simple data and metadata/header container.
    """

    def __init__(self):
        self.info = {}
        self.data = None


class OsisafNarSstProduct(mpop.channel.GenericChannel):

    def __init__(self, filename=None):
        mpop.channel.GenericChannel.__init__(self)
        self.name = "SST"
        self.mda = {}
        #self.header = {}
        self._projectables = []
        self._keys = []
        self._refs = {}

        self.time = None
        self.dtime = None
        self.sec_1981 = None
        self.sec_1970_1981 = time.mktime(
            (1981, 1, 1, 0, 0, 0, 0, 0, 0)) - time.timezone

        self.file = None
        self.sst = None
        self.l2pf = None  # L2P flags
        self.dt = None  # DT_analysis (K): (SST Deviation from previous day)
        self.stdv = None  # Standard deviation
        self.bias = None  # Bias
        self.lon = None
        self.lat = None

        self.shape = None
        if filename:
            self.read(filename)

    def read(self, filename, load_lonlat=True):
        """Read the OSISAF SST netCDF formatet data (from Ifremer)"""
        LOG.debug("OSISAF netCDF file format...")

        self.file = Dataset(filename, 'r')

        self.fillheader()

        # SST (K):
        sst_data = self.file.variables['sea_surface_temperature']

        sstdata = sst_data[0]
        self.sst = InfoObject()
        # For some strange reason the array seem to start in the lower left!?
        self.sst.data = sstdata[::-1]
        self.sst.info = self.get_data_header(self.sst.info, sst_data)
        self._projectables.append('sst')

        # dtime:
        dtime = self.file.variables['sst_dtime']
        dtime_data = dtime[0] * dtime.scale_factor + dtime.add_offset
        dtime_obj = InfoObject()
        dtime_obj.data = dtime_data[::-1]
        dtime_obj.info = self.get_data_header(dtime_obj.info, dtime)
        self.sec_1981 = dtime_obj.data + self.file.variables['time'][0]
        self.dtime = dtime_obj
        self._projectables.append('dtime')

        # DT_analysis (K): (SST Deviation from previous day)
        dta = self.file.variables['dt_analysis']
        gain = 0.1
        nodata = 255
        offset = -12.7
        data = dta[0] * dta.scale_factor + dta.add_offset
        valid_min = dta.valid_min
        valid_max = dta.valid_max

        dt_data = np.where(np.logical_and(np.greater(dta[0], valid_min),
                                          np.less(dta[0], valid_max)),
                           (data - offset) / gain, nodata).astype('B')
        dt = InfoObject()
        dt.data = dt_data[::-1]
        dt.info = self.get_data_header(dt.info, dta)
        dt.info["nodata"] = nodata
        dt.info["gain"] = gain
        dt.info["offset"] = offset
        self.dt = dt
        self._projectables.append('dt')

        # Bias:
        bias = self.file.variables['sses_bias']
        gain = 0.01
        offset = -1.27
        nodata = 255
        x = bias[0] * bias.scale_factor + bias.add_offset
        valid_min = bias.valid_min
        valid_max = bias.valid_max
        bias_data = np.where(np.logical_and(np.greater(bias[0], valid_min),
                                            np.less(bias[0], valid_max)),
                             (x - offset) / gain, nodata).astype('B')

        bias_obj = InfoObject()
        bias_obj.data = bias_data[::-1]
        bias_obj.info = self.get_data_header(bias_obj.info, bias)
        bias_obj.info["nodata"] = nodata
        bias_obj.info["gain"] = gain
        bias_obj.info["offset"] = offset
        self.bias = bias_obj
        self._projectables.append('bias')

        # Standard deviation:
        stdv = self.file.variables['sses_standard_deviation']
        gain = 0.01
        offset = 0.0
        nodata = 255
        x = stdv[0] * stdv.scale_factor + stdv.add_offset
        valid_min = stdv.valid_min
        valid_max = stdv.valid_max
        stdv_data = np.where(np.logical_and(np.greater(stdv[0], valid_min),
                                            np.less(stdv[0], valid_max)),
                             (x - offset) / gain, nodata).astype('B')

        stdv_obj = InfoObject()
        stdv_obj.data = stdv_data[::-1]
        stdv_obj.info = self.get_data_header(stdv_obj.info, stdv)
        stdv_obj.info["nodata"] = nodata
        stdv_obj.info["gain"] = gain
        stdv_obj.info["offset"] = offset
        self.stdv = stdv_obj
        self._projectables.append('stdv')

        # L2P flags:
        l2pf = self.file.variables['l2p_flags'][0]
        l2pf_obj = InfoObject()
        l2pf_obj.data = l2pf[::-1]
        l2pf_obj.info = self.get_data_header(l2pf_obj.info, l2pf)
        self.l2pf = l2pf_obj
        self._projectables.append('l2pf')

        # Longitudes:
        lon = self.file.variables['lon']
        self.lon = InfoObject()
        self.lon.data = lon[::-1].astype('f')
        self.lon.info = self.get_data_header(self.lon.info, lon)

        # Latitudes:
        lat = self.file.variables['lat']
        self.lat = InfoObject()
        self.lat.data = lat[::-1].astype('f')
        self.lat.info = self.get_data_header(self.lat.info, lat)

        return

    def project(self, coverage):
        """Project the data"""
        LOG.debug("Projecting channel %s..." % (self.name))
        import copy
        res = copy.copy(self)

        # Project the data
        for var in self._projectables:
            LOG.info("Projecting " + str(var))
            res.__dict__[var] = copy.copy(self.__dict__[var])
            res.__dict__[var].data = coverage.project_array(
                self.__dict__[var].data)

        res.name = self.name
        res.resolution = self.resolution
        res.filled = True
        res.area = coverage.out_area
        return res

    def is_loaded(self):
        """Tells if the channel contains loaded data.
        """
        return True

    def get_data_header(self, header, datasetObj):
        """Get the attribute names in the datasets that are common"""
        import copy
        retv = copy.copy(header)
        try:
            retv["longt_name"] = datasetObj.long_name
        except AttributeError:
            print "WARNING: No attribute 'long_name'"
        try:
            retv["standard_name"] = datasetObj.standard_name
        except AttributeError:
            print "WARNING: No attribute 'standard_name'"
        try:
            retv["comment"] = datasetObj.comment
        except AttributeError:
            pass
        try:
            retv["units"] = datasetObj.units
        except AttributeError:
            pass
        try:
            retv["valid_min"] = datasetObj.valid_min
        except AttributeError:
            pass
        try:
            retv["valid_max"] = datasetObj.valid_max
        except AttributeError:
            pass
        try:
            retv["scale_factor"] = datasetObj.scale_factor
        except AttributeError:
            pass
        try:
            retv["add_offset"] = datasetObj.add_offset
        except AttributeError:
            pass

        return retv

    def fillheader(self):
        """Get global header info and store in object"""
        pltname = self.file.platform
        self.mda["platform"] = ''.join(pltname.split('_'))
        self.mda["sensor"] = self.file.sensor
        self.mda["spatial_resolution"] = self.file.spatial_resolution
        self.mda["title"] = self.file.title
        self.mda["comment"] = self.file.comment

        self.mda["contact"] = self.file.creator_email
        self.mda["institution"] = self.file.institution
        self.mda["start_time"] = datetime.strptime(
            self.file.start_time, '%Y%m%dT%H%M%SZ')
        self.mda["stop_time"] = datetime.strptime(
            self.file.stop_time, '%Y%m%dT%H%M%SZ')
        # Reference time of sst file in seconds since 1981-01-01 00:00:00:
        self.mda["time"] = self.file.variables['time'][0]
        self.mda['westernmost_longitude'] = self.file.westernmost_longitude
        self.mda['easternmost_longitude'] = self.file.easternmost_longitude
        self.mda['southernmost_latitude'] = self.file.southernmost_latitude
        self.mda['northernmost_latitude'] = self.file.northernmost_latitude


def get_filename(scene, area_name):

    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, scene.fullname + ".cfg"))

    filename = conf.get(scene.instrument_name + "-level4",
                        "sst_product_filename",
                        raw=True,
                        vars=os.environ)
    directory = conf.get(scene.instrument_name + "-level4",
                         "sst_product_dir",
                         vars=os.environ)
    pathname_tmpl = os.path.join(directory, filename)
    LOG.debug("Path = " + str(pathname_tmpl))

    if not scene.orbit:
        orbit = ""
    else:
        orbit = scene.orbit

    filename_tmpl = (scene.time_slot.strftime(pathname_tmpl)
                     % {"area": area_name,
                        "satellite": scene.satname + scene.number})

    product = 'sst'
    file_list = glob.glob(filename_tmpl)
    if len(file_list) > 1:
        LOG.warning("More than 1 file matching for " + product + "! "
                    + str(file_list))
        return None
    elif len(file_list) == 0:
        LOG.warning(
            "No " + product + " matching!: " + filename_tmpl)
        return
    else:
        filename = file_list[0]

    return filename


class OSISAF_SST_Reader(Reader):

    pformat = "nc_osisaf_l2"

    def load(self, satscene, *args, **kwargs):
        """Read data from file and load it into *satscene*.
        """
        lonlat_is_loaded = False

        prodfilename = kwargs.get('filename')

        if "SST" not in satscene.channels_to_load:
            LOG.warning("No SST product requested. Nothing to be done...")
            return

        try:
            area_name = satscene.area_id or satscene.area.area_id
        except AttributeError:
            area_name = "satproj_?????_?????"

        conf = ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

        # Reading the products
        product = "sst"
        classes = {product: OsisafNarSstProduct}

        LOG.debug("Loading " + product)

        if isinstance(prodfilename, (list, tuple, set)):
            for fname in prodfilename:
                kwargs['filename'] = fname
                self.load(satscene, *args, **kwargs)
            return
        elif (prodfilename and
              'OSISAF-L3C' in os.path.basename(prodfilename)):
            if os.path.basename(prodfilename).split("_")[2] == 'SST':
                filename = prodfilename
            else:
                LOG.warning(
                    "Product file name is not as expected: " + str(prodfilename))
                return
        else:
            filename = get_filename(satscene, area_name)

        LOG.info("Filename = " + str(filename))

        chn = classes[product]()
        chn.read(filename, lonlat_is_loaded == False)
        satscene.channels.append(chn)

        lons, lats = chn.lon.data, chn.lat.data
        lonlat_is_loaded = True

        nodata_mask = False
        try:
            from pyresample import geometry
            lons = np.ma.masked_array(lons, nodata_mask)
            lats = np.ma.masked_array(lats, nodata_mask)
            area = geometry.SwathDefinition(lons=lons,
                                            lats=lats)
        except ImportError:
            area = None

        if area:
            chn.area = area
        else:
            chn.lon = lons
            chn.lat = lats

        LOG.info("Loading OSISAF SST parameters done")

        return
