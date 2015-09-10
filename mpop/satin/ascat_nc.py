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
              "satellite": satscene.fullname
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
        data = np.ma.masked_array(ds.variables[chn_name][:],np.isnan(ds.variables[chn_name][:]))
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
