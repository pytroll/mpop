#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Adam.Dybbroe

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

"""A reader for the EUMETSAT MPEF OCA cloud parameters. Data are segmented LRIT
encoded Grib messages

"""

import tempfile
import os.path
from ConfigParser import ConfigParser
from mpop.satin.gribformat import Grib
import mpop.channel
from mpop import CONFIG_PATH
from mpop.plugin_base import Reader
from trollsift import parser
from glob import glob
import pyresample as pr
import numpy as np

import logging
LOG = logging.getLogger(__name__)

CFG_DIR = os.environ.get('PPP_CONFIG_DIR', './')
AREA_DEF_FILE = os.path.join(CFG_DIR, "areas.def")
if not os.path.exists(AREA_DEF_FILE):
    raise IOError('Config file %s does not exist!' % AREA_DEF_FILE)


LRIT_PATTERN = "L-000-{platform_name:_<5s}_-MPEF________-OCAE_____-{segment:_<9s}-{nominal_time:%Y%m%d%H%M}-{compressed:_<2s}"


SCENE_TYPE_LAYERS = {111: 'Single Layer Water Cloud',
                     112: 'Single Layer Ice Cloud',
                     113: 'Multi Layer Cloud'}

OCA_FIELDS = [{'Pixel scene type': 'Scene type'},
              {'24': 'Measurement Cost',
               'abbrev': 'JM', 'units': ''},
              {'25': 'Upper Layer Cloud Optical Thickness', 'units': '',
               'abbrev': 'ULCOT'},
              {'26': 'Upper Layer Cloud Top Pressure', 'units': 'Pa',
               'abbrev': 'ULCTP'},
              {'27': 'Upper Layer Cloud Effective Radius', 'units': 'm',
               'abbrev': 'ULCRE'},
              {'28': 'Error in Upper Layer Cloud Optical Thickness', 'units': '',
               'abbrev': 'ERR-ULCOT'},
              {'29': 'Error in Upper Layer Cloud Top Pressure', 'units': 'Pa',
               'abbrev': 'ERR-ULCTP'},
              {'30': 'Error in Upper Layer Cloud Effective Radius', 'units': 'm',
               'abbrev': 'ERR-ULCRE'},
              {'31': 'Lower Layer Cloud Optical Thickness',
                  'units': '', 'abbrev': 'LLCOT'},
              {'32': 'Lower Layer Cloud Top Pressure',
                  'units': 'Pa', 'abbrev': 'LLCTP'},
              {'33': 'Error in Lower Layer Cloud Optical Thickness',
                  'units': '', 'abbrev': 'ERR-LLCOT'},
              {'34': 'Error in Lower Layer Cloud Top Pressure', 'units': 'Pa',
               'abbrev': 'ERR-LLCTP'}]

FIELDNAMES = {'scenetype': ('Pixel scene type', None),
              'cost': ('24', None),
              'ul_cot': ('25', '28'),
              'ul_ctp': ('26', '29'),
              'reff': ('27', '30'),
              'll_cot': ('31', '33'),
              'll_ctp': ('32', '34')}


class OCAField(object):

    """One OCA data field with metadata"""

    def __init__(self, units=None, long_name='', standard_name=''):
        self.units = units
        self.data = None
        self.error = None
        self.long_name = long_name
        self.standard_name = standard_name
        self.info = {}


class OCAData(mpop.channel.GenericChannel):

    """The OCA scene data"""

    def __init__(self):
        mpop.channel.GenericChannel.__init__(self)
        self.name = "OCA"
        self.mda = {}
        self._keys = []
        self._refs = {}

        self._lritfiles = None
        self._gribfilename = None
        self._store_grib = False

        self.resolution = 3000

        self.scenetype = OCAField()
        self.cost = OCAField()
        self.ul_cot = OCAField()
        self.ll_cot = OCAField()
        self.ul_ctp = OCAField()
        self.ll_ctp = OCAField()
        self.reff = OCAField()
        self._projectables = []
        for field in FIELDNAMES.keys():
            self._projectables.append(field)

        self.timeslot = None
        self.area_def = pr.utils.load_area(AREA_DEF_FILE, 'met09globeFull')
        self.shape = None

    def readgrib(self):
        """Read the data"""

        oca = Grib(self._gribfilename)
        self.scenetype.data = oca.get('Pixel scene type')[::-1, ::-1]
        self.scenetype.long_name = OCA_FIELDS[0]['Pixel scene type']

        for field in FIELDNAMES.keys():

            setattr(getattr(self, field), 'data', oca.get(
                FIELDNAMES[field][0])[::-1, ::-1])
            param = [s for s in OCA_FIELDS if FIELDNAMES[field][0] in s][0]
            if 'units' in param:
                setattr(getattr(self, field), 'units', param['units'])
            if 'abbrev' in param:
                setattr(getattr(self, field), 'standard_name', param['abbrev'])
            setattr(getattr(self, field), 'long_name',
                    param[FIELDNAMES[field][0]])
            param_name = FIELDNAMES[field][1]
            if param_name:
                setattr(
                    getattr(self, field), 'error', oca.get(param_name)[::-1, ::-1])

        if not self._store_grib:
            os.remove(self._gribfilename)

    def read_from_lrit(self, filenames, gribfilename=None):
        """Read and concatenate the LRIT segments"""

        self._lritfiles = filenames

        if len(filenames) == 0:
            print("No files provided!")
            return

        if gribfilename:
            self._store_grib = True
            self._gribfilename = gribfilename
        else:
            self._store_grib = False
            self._gribfilename = tempfile.mktemp(suffix='.grb')

        p__ = parser.Parser(LRIT_PATTERN)

        bstr = {}
        nsegments = 0
        for lritfile in self._lritfiles:
            if os.path.basename(lritfile).find('PRO') > 0:
                print("PRO file... %s: Skip it..." % lritfile)
                continue

            res = p__.parse(os.path.basename(lritfile))
            segm = int(res['segment'].strip('_'))
            if not self.timeslot:
                self.timeslot = res['nominal_time']
            LOG.debug("Segment = %d", segm)
            nsegments = nsegments + 1

            with open(lritfile) as fpt:
                fpt.seek(103)
                bstr[segm] = fpt.read()

        fstr = bstr[1]
        for idx in range(2, nsegments + 1):
            fstr = fstr + bstr[idx]

        with open(self._gribfilename, 'wb') as fpt:
            fpt.write(fstr)

        self.readgrib()

    def project(self, coverage):
        """Project the data"""
        LOG.debug("Projecting channel %s...", (self.name))
        import copy
        res = copy.copy(self)

        res.name = self.name
        res.resolution = self.resolution
        res.filled = True
        res.area = coverage.out_area
        resolution_str_x = str(int(res.area.pixel_size_x)) + 'm'
        resolution_str_y = str(int(res.area.pixel_size_y)) + 'm'

        time_axis = 0

        # Project the data
        for var in self._projectables:
            LOG.info("Projecting " + str(var))
            res.__dict__[var] = copy.copy(self.__dict__[var])
            data = coverage.project_array(self.__dict__[var].data)
            valid_min = np.min(data)
            valid_max = np.max(data)

            res.__dict__[var].data = data
            res.__dict__[var].info['var_name'] = var
            res.__dict__[var].info[
                'var_data'] = np.ma.expand_dims(data, time_axis)

            dim_names = ['y' + resolution_str_y,
                         'x' + resolution_str_x]
            dim_names.insert(0, 'time')
            res.__dict__[var].info['var_dim_names'] = dim_names
            res.__dict__[var].info['long_name'] = res.__dict__[var].long_name
            res.__dict__[var].info[
                'standard_name'] = res.__dict__[var].standard_name
            res.__dict__[var].info['valid_range'] = np.array(
                [valid_min, valid_max])
            # res.__dict__[var].info['resolution'] = res.resolution

        return res

    def is_loaded(self):
        """Tells if the channel contains loaded data.
        """
        return True


class OCAReader(Reader):

    pformat = "mpef_oca"

    def load(self, satscene, *args, **kwargs):
        """Read data from files and load it into *satscene*.
        """
        lonlat_is_loaded = False

        lritfiles = kwargs.get('filenames')

        if "OCA" not in satscene.channels_to_load:
            LOG.warning("No OCA product requested. Nothing to be done...")
            return

        area_name = satscene.area_id or satscene.area.area_id
        # platform_name = satscene.satname

        conf = ConfigParser()
        conf.read(os.path.join(CONFIG_PATH, satscene.fullname + ".cfg"))

        # Reading the products
        product = "oca"
        classes = {product: OCAData}

        LOG.debug("Loading " + product)

        if not lritfiles:
            dummy, lritfiles = get_lrit_filenames(satscene, area_name)

        LOG.info("Filenames = " + str(lritfiles))

        chn = classes[product]()
        chn.read_from_lrit(lritfiles)

        # Prepare info object for netCDF writer:
        resolution_str = str(int(chn.resolution)) + 'm'
        for field in chn._projectables:

            getattr(chn, field).info['var_name'] = field
            getattr(chn, field).info['var_data'] = getattr(chn, field).data
            getattr(chn, field).info[
                'var_dir_names'] = getattr(chn, field).data

            getattr(chn, field).info['var_dim_names'] = ('y' + resolution_str,
                                                         'x' + resolution_str)
            getattr(chn, field).info['long_name'] = getattr(
                chn, field).long_name
            getattr(chn, field).info['standard_name'] = getattr(
                chn, field).standard_name
            valid_min = np.min(getattr(chn, field).data)
            valid_max = np.max(getattr(chn, field).data)
            getattr(chn, field).info['valid_range'] = np.array(
                [valid_min, valid_max])
            getattr(chn, field).info['resolution'] = chn.resolution

        satscene.channels.append(chn)

        LOG.info("Loading MPEF OCA cloud parameters done")

        return


def get_lrit_filenames(scene, area_name):
    """Get the set of lrit filenames for the given scene
    """

    conf = ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, scene.fullname + ".cfg"))

    filename = conf.get(scene.instrument_name + "-level4",
                        "filename",
                        raw=True,
                        vars=os.environ)
    directory = conf.get(scene.instrument_name + "-level4",
                         "dir",
                         vars=os.environ)
    pathname_tmpl = os.path.join(directory, filename)
    LOG.debug("Path = " + str(pathname_tmpl))

    fparser = parser.Parser(pathname_tmpl)

    lrit_files = glob(
        parser.globify(pathname_tmpl, {'nominal_time': scene.time_slot}))

    prologue = None
    segmfiles = []
    segm_numbers = []
    for item in lrit_files:
        p__ = fparser.parse(item)
        segm = p__['segment'].strip('_')
        if segm == 'PRO':
            prologue = item
        else:
            segm_numbers.append(int(segm))
            segmfiles.append(item)

    if not prologue:
        LOG.warning("No prologue file found for timeslot")

    segm_numbers.sort()
    if range(1, 11) == segm_numbers:
        LOG.info("All ten segment files found")
    else:
        LOG.warning("Less than 10 segments found: %s", str(segm_numbers))

    return prologue, segmfiles
