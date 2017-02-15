#!/usr/bin/python
"""Loader for s2, jpeg2000 format.
"""
#Matias Takala FMI 2016

import glob
import os
import pickle
import re
from ConfigParser import ConfigParser

import numpy.ma as ma
from pyresample import utils

import glymur
from mpop import CONFIG_PATH
from mpop.satellites import GenericFactory

#in this version Q_V is hardcoded but could be read from metadata
QUANTIFICATION_VALUE = 10000


def parse_tile(file):
    tile = re.findall('T(\d{2}\w{3})_', file)
    f = open('s2tiles.pickle', 'r')
    s2tiles = pickle.load(f)
    f.close()
    return [tile[0], s2tiles[tile[0]]]


def read_jp2_data(file):
    jp2 = glymur.Jp2k(file)
    data = jp2[:] / (QUANTIFICATION_VALUE + 0.0)
    return data


def open_s2_tile(fname):
    data = read_jp2_data(fname)
    size = data.shape
    params = parse_tile(fname)
    areadef = utils.get_area_def(
        params[0], "Sentinel 2 tile " + params[0], 'PROJ EPSG:' + params[1][0],
        'init=epsg:' + params[1][0], size[0], size[1], params[1][1])
    return ([data, areadef])


def load(satscene):
    """Load jpeg2000 data.
    """

    # Read config file content
    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

    for chn_name in satscene.channels_to_load:
        values = {"orbit": satscene.orbit,
                  "satname": satscene.satname.upper(),
                  "number": satscene.number,
                  "instrument": satscene.instrument_name.upper(),
                  "satellite": satscene.fullname.upper(),
                  "band": chn_name}
        filename = os.path.join(
            conf.get("msi-level2", "dir"),
            satscene.time_slot.strftime(conf.get(
                "msi-level2", "filename", raw=True)) % values)
        filelist = glob.glob(filename)
        data_area = open_s2_tile(filelist[0])
        satscene[chn_name] = ma.masked_array(data_area[0])
        satscene[chn_name].area = data_area[1]
