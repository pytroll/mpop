#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014, 2015 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c14526.ad.smhi.se>
#   Panu Lahtinen <panu.lahtinen@fmi.fi>

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

"""PPS netcdf cloud product reader
"""

import os.path
from ConfigParser import ConfigParser
from ConfigParser import NoOptionError

from datetime import datetime
import glob
import numpy as np

import mpop.channel
from mpop import CONFIG_PATH
from mpop.plugin_base import Reader

import logging
LOG = logging.getLogger(__name__)

NEW_PRODNAMES = {'cloudtype': 'CT',
                 'cloudmask': 'CMA',
                 'precipclouds': 'PC',
                 'cpp': 'CPP',
                 'ctth': 'CTTH'}

PPS_DATASETS = ['Cloud Type',
                'Multilayer Cloud Detection',
                "SAFNWC PPS PC likelihood of intense precipitation",
                "SAFNWC PPS PC likelihood of moderate precipitation",
                "SAFNWC PPS PC likelihood of light precipitation",
                ]

class InfoObject(object):

    """Simple data and info container.
    """

    def __init__(self):
        self.info = {}
        self.data = None


class NwcSafPpsChannel(mpop.channel.GenericChannel):

    def __init__(self, filename=None):
        mpop.channel.GenericChannel.__init__(self)
        self.mda = {}
        self._projectables = []
        self._keys = []
        self._refs = {}
        self.shape = None
        if filename:
            self.read(filename)

    def read(self, filename, load_lonlat=True):
        """Read the PPS v2014 formatet data"""
        LOG.debug("New netCDF CF file format!")

        import h5py

        # Try see if it is bzipped:
        import bz2
        if filename.endswith('bz2'):
            bz2file = bz2.BZ2File(filename)
            import tempfile
            tmpfilename = tempfile.mktemp()
            try:
                ofpt = open(tmpfilename, 'wb')
                ofpt.write(bz2file.read())
                ofpt.close()
            except IOError:
                import traceback
                traceback.print_exc()
                LOG.info("Failed to read bzipped file %s", str(filename))

            filename = tmpfilename

        h5f = h5py.File(filename, 'r')
        self.mda.update(h5f.attrs.items())

        self.mda["satellite"] = h5f.attrs['platform']
        self.mda["orbit"] = h5f.attrs['orbit_number']
        try:
            self.mda["time_slot"] = \
                datetime.strptime(h5f.attrs['time_coverage_start'][:-2],
                                  "%Y%m%dT%H%M%S")
        except AttributeError:
            LOG.debug("No time information in product file!")

        variables = {}

        for key, item in h5f.items():
            if item.attrs.get("CLASS") != 'DIMENSION_SCALE':
                variables[key] = item

        # processed variables
        processed = set()

        non_processed = set(variables.keys()) - processed

        for var_name in non_processed:
            if var_name in ['lon', 'lat', 'lon_reduced', 'lat_reduced']:
                continue

            var = variables[var_name]
            if ("standard_name" not in var.attrs.keys() and
                "long_name" not in var.attrs.keys()):
                LOG.info("Delayed processing of %s", var_name)
                continue

            # Don't know how to unambiguously decide if the array is really a
            # data array or a palette or something else!
            # FIXME!
            if "standard_name" in var.attrs.keys():
                self._projectables.append(var_name)
            elif "long_name" in var.attrs.keys():
                dset_found = False
                for item in PPS_DATASETS:
                    if var.attrs['long_name'].find(item) >= 0:
                        self._projectables.append(var_name)
                        dset_found = True
                        break
                if not dset_found:
                    self.mda[var_name] = var[:]
                    # try:
                    #     self.mda[var_name] = var[:].filled(0)
                    # except AttributeError:
                    #     self.mda[var_name] = var[:]
                    continue

            setattr(self, var_name, InfoObject())
            # for key, item in var.attrs.items():
            for key in var.attrs.keys():
                if key != "DIMENSION_LIST":
                    item = var.attrs[key]
                    getattr(self, var_name).info[key] = item

            data = var[:]
            if len(data.shape) == 3 and data.shape[0] == 1:
                LOG.debug("Rip of the first dimension of length 1")
                data = data[0]

            if 'valid_range' in var.attrs.keys():
                data = np.ma.masked_outside(data, *var.attrs['valid_range'])
            elif '_FillValue' in var.attrs.keys():
                data = np.ma.masked_where(data, var.attrs['_FillValue'])
            dataset = (data * var.attrs.get("scale_factor", 1)
                       + var.attrs.get("add_offset", 0))

            getattr(self, var_name).data = dataset

            LOG.debug("long_name: %s", str(var.attrs['long_name']))
            LOG.debug("Var = %s, shape = %s",
                      str(var_name), str(dataset.shape))

            if self.shape is None:
                self.shape = dataset.shape
            elif self.shape != dataset.shape:
                LOG.debug("Shape = %s. Not the same shape as previous field.",
                          str(dataset.shape))
                #raise ValueError("Different variable shapes !")

            #dims = var.dimensions
            #dim = dims[0]

            processed |= set([var_name])

        non_processed = set(variables.keys()) - processed
        if len(non_processed) > 0:
            LOG.warning("Remaining non-processed variables: %s",
                        str(non_processed))

        # Get lon,lat:
        # from pyresample import geometry
        # area = geometry.SwathDefinition(lons=lon, lats=lat)

        return

    def project(self, coverage):
        """Project the data"""
        LOG.debug("Projecting channel %s...", self.name)
        import copy
        res = copy.copy(self)

        # Project the data
        for var in self._projectables:
            LOG.info("Projecting %s", str(var))
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
        # return len(self._projectables) > 0

    def save(self, filename, old=True, **kwargs):
        """Save to old format"""
        del kwargs
        if old:
            try:
                from nwcsaf_formats.ppsv2014_to_oldformat import write_product
                write_product(self, filename)
            except ImportError:
                LOG.error("Could not save to old format")
                raise
        else:
            raise NotImplementedError("Can't save to new pps format yet.")


class PPSReader(Reader):
    """Reader class for PPS files"""
    pformat = "nc_pps_l2"

    def load(self, satscene, *args, **kwargs):
        """Read data from file and load it into *satscene*.
        """
        lonlat_is_loaded = False

        geofilename = kwargs.get('geofilename')
        prodfilename = kwargs.get('filename')

        products = []
        if "CTTH" in satscene.channels_to_load:
            products.append("ctth")
        if "CT" in satscene.channels_to_load:
            products.append("cloudtype")
        if "CMA" in satscene.channels_to_load:
            products.append("cloudmask")
        if "PC" in satscene.channels_to_load:
            products.append("precipclouds")
        if "CPP" in satscene.channels_to_load:
            products.append("cpp")

        if len(products) == 0:
            return

        try:
            area_name = satscene.area_id or satscene.area.area_id
        except AttributeError:
            area_name = "satproj_?????_?????"

        # Looking for geolocation file

        conf = ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

        try:
            geodir = conf.get(satscene.instrument_name + "-level3",
                              "cloud_product_geodir",
                              vars=os.environ)
        except NoOptionError:
            LOG.warning("No option 'geodir' in level3 section")
            geodir = None

        if not geofilename and geodir:
            # Load geo file from config file:
            try:
                if not satscene.orbit:
                    orbit = ""
                else:
                    orbit = satscene.orbit
                geoname_tmpl = conf.get(satscene.instrument_name + "-level3",
                                        "cloud_product_geofilename",
                                        raw=True,
                                        vars=os.environ)
                filename_tmpl = satscene.time_slot.strftime(geoname_tmpl) % \
                                {"orbit": str(orbit).zfill(5) or "*",
                                 "area": area_name,
                                 "satellite": satscene.satname + \
                                 satscene.number}

                file_list = glob.glob(os.path.join(geodir, filename_tmpl))
                if len(file_list) > 1:
                    LOG.warning("More than 1 file matching for geoloaction: %s",
                                str(file_list))
                elif len(file_list) == 0:
                    LOG.warning("No geolocation file matching!: %s",
                                filename_tmpl)
                else:
                    geofilename = file_list[0]
            except NoOptionError:
                geofilename = None

        # Reading the products

        classes = {"ctth": CloudTopTemperatureHeight,
                   "cloudtype": CloudType,
                   "cloudmask": CloudMask,
                   "precipclouds": PrecipitationClouds,
                   "cpp": CloudPhysicalProperties
                   }

        nodata_mask = False

        area = None
        lons = None
        lats = None
        chn = None
        shape = None
        read_external_geo = {}
        for product in products:
            LOG.debug("Loading %s", product)

            if isinstance(prodfilename, (list, tuple, set)):
                for fname in prodfilename:
                    kwargs['filename'] = fname
                    self.load(satscene, *args, **kwargs)
                return
            elif (prodfilename and
                  (os.path.basename(prodfilename).startswith('S_NWC') or
                   os.path.basename(prodfilename).startswith("W_XX-EUMETSAT"))):
                if NEW_PRODNAMES[product] in os.path.basename(prodfilename):
                    filename = prodfilename
                else:
                    continue
            else:
                filename = conf.get(satscene.instrument_name + "-level3",
                                    "cloud_product_filename",
                                    raw=True,
                                    vars=os.environ)
                directory = conf.get(satscene.instrument_name + "-level3",
                                     "cloud_product_dir",
                                     vars=os.environ)
                pathname_tmpl = os.path.join(directory, filename)
                LOG.debug("Path = %s", str(pathname_tmpl))

                if not satscene.orbit:
                    orbit = ""
                else:
                    orbit = satscene.orbit

                filename_tmpl = \
                    satscene.time_slot.strftime(pathname_tmpl) % \
                    {"orbit": str(orbit).zfill(5) or "*",
                     "area": area_name,
                     "satellite": satscene.satname + satscene.number,
                     "product": product}

                file_list = glob.glob(filename_tmpl)
                if len(file_list) == 0:
                    product_name = NEW_PRODNAMES.get(product, product)
                    LOG.info("No %s product in old format matching",
                             str(product))
                    filename_tmpl = \
                        satscene.time_slot.strftime(pathname_tmpl) % \
                        {"orbit": str(orbit).zfill(5) or "*",
                         "area": area_name,
                         "satellite": satscene.satname + satscene.number,
                         "product": product_name}

                    file_list = glob.glob(filename_tmpl)

                if len(file_list) > 1:
                    LOG.warning("More than 1 file matching for %s: %s",
                                str(product), str(file_list))
                    continue
                elif len(file_list) == 0:
                    LOG.warning("No %s matching: %s",
                                str(product), filename_tmpl)
                    continue
                else:
                    filename = file_list[0]

            chn = classes[product]()
            chn.read(filename, lonlat_is_loaded == False)

            try:
                chn.lat_reduced = None
                chn.lon_reduced = None
            except AttributeError:
                pass

            # concatenate if there's already a channel of the same type (name)
            # NOTE! There must be a better way to do this
            if chn.name in satscene:
                LOG.info("Concatenating data")
                for atr in dir(satscene[chn.name]):
                    try:
                        old_data = getattr(satscene[chn.name], atr)
                        new_data = getattr(chn, atr)
                        old_data.data = np.concatenate((old_data.data,
                                                        new_data.data))
                    except AttributeError:
                        pass
            else:
                LOG.info("Adding new channel %s", chn.name)
                satscene.channels.append(chn)

            # Check if geolocation is loaded:
            if not chn.area:
                read_external_geo[product] = satscene.channels[-1].name
                shape = chn.shape

        # Check if some 'channel'/product needs geolocation. If some
        # product does not have geolocation, get it from the
        # geofilename:
        if not read_external_geo:
            LOG.info("Loading PPS parameters done.")
            return

        # Load geolocation
        interpolate = False
        if geofilename:
            geodict = get_lonlat(geofilename)
        else:
            LOG.info("No Geo file specified: "
                     "Geolocation will be loaded from product")
            geodict = get_lonlat(filename)

        lons, lats = geodict['lon'], geodict['lat']
        if lons.shape != shape or lats.shape != shape:
            interpolate = True
            row_indices = geodict['row_indices']
            column_indices = geodict['col_indices']

        if interpolate:
            from geotiepoints import SatelliteInterpolator

            cols_full = np.arange(shape[1])
            rows_full = np.arange(shape[0])

            satint = SatelliteInterpolator((lons, lats),
                                           (row_indices,
                                            column_indices),
                                           (rows_full, cols_full))
            #satint.fill_borders("y", "x")
            lons, lats = satint.interpolate()

        try:
            from pyresample import geometry
            lons = np.ma.masked_array(lons, nodata_mask)
            lats = np.ma.masked_array(lats, nodata_mask)
            area = geometry.SwathDefinition(lons=lons,
                                            lats=lats)
        except ImportError:
            area = None

        for chn_name in read_external_geo.values():
            if area:
                try:
                    lats = np.concatenate((satscene[chn_name].area.lats, lats))
                    lons = np.concatenate((satscene[chn_name].area.lons, lons))
                    satscene[chn_name].area = \
                        geometry.SwathDefinition(lons=lons, lats=lats)
                except AttributeError:
                    satscene[chn_name].area = area
            else:
                satscene[chn_name].lat = \
                    np.concatenate((satscene[chn_name].lat, lats))
                satscene[chn_name].lon = \
                    np.concatenate((satscene[chn_name].lon, lons))

        LOG.info("Loading PPS parameters done.")

        return


class CloudType(NwcSafPpsChannel):
    """CloudType PPS channel object"""
    def __init__(self, filename=None):
        NwcSafPpsChannel.__init__(self, filename)
        self.name = "CT"


class CloudTopTemperatureHeight(NwcSafPpsChannel):
    """Cloud top temperature and height PPS channel object"""
    def __init__(self, filename=None):
        NwcSafPpsChannel.__init__(self, filename)
        self.name = "CTTH"


class CloudMask(NwcSafPpsChannel):
    """Cloud mask PPS channel object"""
    def __init__(self, filename=None):
        NwcSafPpsChannel.__init__(self, filename)
        self.name = "CMA"


class PrecipitationClouds(NwcSafPpsChannel):
    """Precipitation clouds PPS channel object"""
    def __init__(self, filename=None):
        NwcSafPpsChannel.__init__(self, filename)
        self.name = "PC"


class CloudPhysicalProperties(NwcSafPpsChannel):
    """Cloud physical proeperties PPS channel"""
    def __init__(self, filename=None):
        NwcSafPpsChannel.__init__(self, filename)
        self.name = "CPP"


def get_lonlat(filename):
    """Read lon,lat from the given netCDF4 CF file"""
    import h5py

    col_indices = None
    row_indices = None

    LOG.debug("Geolocations read from %s", filename)

    h5f = h5py.File(filename, 'r')

    try:
        lon = h5f['lon']
        lons = (lon[:] * lon.attrs.get("scale_factor", 1)
                + lon.attrs.get("add_offset", 0))
        lat = h5f['lat']
        lats = (lat[:] * lat.attrs.get("scale_factor", 1)
                + lat.attrs.get("add_offset", 0))
    except KeyError:
        lon = h5f['lon_reduced']
        lons = lon[:]
        lat = h5f['lat_reduced']
        lats = lat[:]

    lons = np.ma.masked_equal(lons, lon.attrs["_FillValue"])
    lats = np.ma.masked_equal(lats, lat.attrs["_FillValue"])

    # FIXME: this is to mask out the npp bowtie deleted pixels...
    if "NPP" in h5f.attrs['platform']:

        new_mask = np.zeros((16, 3200), dtype=bool)
        new_mask[0, :1008] = True
        new_mask[1, :640] = True
        new_mask[14, :640] = True
        new_mask[15, :1008] = True
        new_mask[14, 2560:] = True
        new_mask[1, 2560:] = True
        new_mask[0, 2192:] = True
        new_mask[15, 2192:] = True
        new_mask = np.tile(new_mask, (lons.shape[0] / 16, 1))
        lons = np.ma.masked_where(new_mask, lons)
        lats = np.ma.masked_where(new_mask, lats)

    if "column_indices" in h5f.keys():
        col_indices = h5f["column_indices"][:]
    if "row_indices" in h5f.keys():
        row_indices = h5f["row_indices"][:]
    if "nx_reduced" in h5f:
        col_indices = h5f["nx_reduced"][:]
    if "ny_reduced" in h5f:
        row_indices = h5f["ny_reduced"][:]

    return {'lon': lons,
            'lat': lats,
            'col_indices': col_indices, 'row_indices': row_indices}
