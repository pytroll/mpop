#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2009, 2010, 2011, 2012, 2013, 2014.

# SMHI,
# Folkborgsvägen 1,
# Norrköping,
# Sweden

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This file is part of mpop.

# mpop is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# mpop is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with mpop.  If not, see <http://www.gnu.org/licenses/>.

"""This module handles coverage objects. Such objects are used to
transform area projected data by changing either the area or the
projection or both. A typical usage is to transform one large area in
satellite projection to an area of interest in polar projection for
example.
"""
import os
import ConfigParser
import logging

import numpy as np
from pyresample import image, utils, geometry, kd_tree
from pyresample.bilinear import get_sample_from_bil_info, get_bil_info
try:
    from pyresample.ewa import ll2cr, fornav
except ImportError:
    ll2cr, fornav = None, None

from mpop import CONFIG_PATH

logger = logging.getLogger(__name__)

area_file = None


def get_area_file():
    global area_file
    if area_file:
        return area_file

    conf = ConfigParser.ConfigParser()
    conf.read(os.path.join(CONFIG_PATH, "mpop.cfg"))

    try:
        area_file = os.path.join(conf.get("projector",
                                          "area_directory") or
                                 CONFIG_PATH,
                                 conf.get("projector", "area_file"))
    except ConfigParser.NoSectionError:
        area_file = ""
        logger.warning("Couldn't find the mpop.cfg file. "
                       "Do you have one ? is it in $PPP_CONFIG_DIR ?")
    return area_file


def get_area_def(area_name):
    """Get the definition of *area_name* from file. The file is defined to use
    is to be placed in the $PPP_CONFIG_DIR directory, and its name is defined
    in mpop's configuration file.
    """
    return utils.parse_area_file(get_area_file(), area_name)[0]


def _get_area_hash(area):
    """Calculate a (close to) unique hash value for a given area.
    """
    try:
        return hash(area.lons.tostring() + area.lats.tostring())
    except AttributeError:
        try:
            return hash(area.tostring())
        except AttributeError:
            return hash(str(area))


class Projector(object):

    """This class define projector objects. They contain the mapping
    information necessary for projection purposes. For efficiency reasons,
    generated projectors can be saved to disk for later reuse. Use the
    :meth:`save` method for this.

    To define a projector object, on has to specify *in_area* and
    *out_area*, and can also input the *in_lonlats* or the *mode*.
    Available modes area:
    - 'quick' (works only if both in- and out-areas are AreaDefinitions)
    - 'bilinear' (out-area needs to be AreaDefinition with proj4_string)
    - 'ewa'
    - 'nearest'.
    *radius* defines the radius of influence for nearest neighbour
    search in 'nearest' and 'bilinear' modes.
    """

    def __init__(self, in_area, out_area,
                 in_latlons=None, mode=None,
                 radius=10000, nprocs=1):

        if (mode is not None and
                mode not in ["quick", "nearest", "ewa", "bilinear"]):
            raise ValueError("Projector mode must be one of 'nearest', "
                             "'quick', 'ewa', 'bilinear'")

        self.area_file = get_area_file()

        self.in_area = None
        self.out_area = None
        self._cache = None
        self._filename = None
        self.mode = "quick"
        self.radius = radius
        self.conf = ConfigParser.ConfigParser()
        self.conf.read(os.path.join(CONFIG_PATH, "mpop.cfg"))

        # TODO:
        # - Rework so that in_area and out_area can be lonlats.
        # - Add a recompute flag ?

        # Setting up the input area
        in_id = self._setup_input_area(in_area, latlons=in_latlons)

        # Setting up the output area
        out_id = self._setup_output_area(out_area, latlons=None)

        # if self.in_area == self.out_area:
        #    return

        # choosing the right mode if necessary
        if mode is None:
            try:
                dicts = in_area.proj_dict, out_area.proj_dict
                del dicts
                self.mode = "quick"
            except AttributeError:
                self.mode = "nearest"
        else:
            self.mode = mode

        projections_directory = "/var/tmp"
        try:
            projections_directory = self.conf.get("projector",
                                                  "projections_directory")
        except ConfigParser.NoSectionError:
            pass

        self._filename = get_precompute_cache_fname(in_id, out_id,
                                                    in_area, out_area,
                                                    self.mode,
                                                    projections_directory)

        try:
            self._cache = {}
            self._file_cache = np.load(self._filename)
        except:
            logger.info("Computing projection from %s to %s...",
                        in_id, out_id)

            if self.mode == "nearest":
                self._cache = calc_nearest_params(in_area, out_area,
                                                  radius, nprocs=nprocs)

            elif self.mode == "quick":
                self._cache = calc_quick_params(in_area, out_area)

            elif self.mode == "ewa":
                if ll2cr is not None:
                    self._cache = calc_ewa_params(in_area, out_area)
                else:
                    raise ImportError("Can't import pyresample.ewa")

            elif self.mode == "bilinear":
                self._cache = calc_bilinear_params(in_area, out_area,
                                                   radius, nprocs=nprocs)

    def _setup_input_area(self, area, latlons=None):
        """Setup self.in_area and return area id"""
        try:
            self.in_area, in_id = get_area_and_id(area, latlons=latlons)
        except TypeError:
            raise utils.AreaNotFound("Input area " +
                                     str(area) +
                                     " must be defined in " +
                                     self.area_file +
                                     ", be an area object"
                                     " or longitudes/latitudes must be "
                                     "provided.")

        return in_id

    def _setup_output_area(self, area, latlons=None):
        """Setup output area"""
        try:
            self.out_area, out_id = get_area_and_id(area, latlons=latlons)
        except AttributeError:
            raise utils.AreaNotFound("Output area " +
                                     str(area) +
                                     " must be defined in " +
                                     self.area_file + " or "
                                     "be an area object.")
        return out_id

    def save(self, resave=False):
        """Save the precomputation to disk, and overwrite existing file in case
        *resave* is true.
        """
        if (not os.path.exists(self._filename)) or resave:
            logger.info("Saving projection to " +
                        self._filename)
            np.savez(self._filename, **self._cache)

    def _project_array_nearest(self, data):
        """Project array *data* using nearest neighbour resampling"""
        if 'valid_index' not in self._cache:
            self._cache['valid_index'] = self._file_cache['valid_index']
            self._cache['valid_output_index'] = \
                self._file_cache['valid_output_index']
            self._cache['index_array'] = self._file_cache['index_array']

        valid_index, valid_output_index, index_array = \
            (self._cache['valid_index'],
             self._cache['valid_output_index'],
             self._cache['index_array'])

        res = kd_tree.get_sample_from_neighbour_info('nn',
                                                     self.out_area.shape,
                                                     data,
                                                     valid_index,
                                                     valid_output_index,
                                                     index_array,
                                                     fill_value=None)
        return res

    def _project_array_quick(self, data):
        """Project array *data* using quick interpolation"""
        if 'row_idx' not in self._cache:
            self._cache['row_idx'] = self._file_cache['row_idx']
            self._cache['col_idx'] = self._file_cache['col_idx']
        row_idx, col_idx = self._cache['row_idx'], self._cache['col_idx']
        img = image.ImageContainer(data, self.in_area, fill_value=None)
        res = np.ma.array(img.get_array_from_linesample(row_idx, col_idx),
                          dtype=data.dtype)

        return res

    def _project_array_ewa(self, data):
        """Project array *data* using EWA interpolation"""
        # TODO: should be user configurable?
        rows_per_scan = None

        if 'ewa_cols' in self. not_cache:
            self._cache['ewa_cols'] = self._file_cache['ewa_cols']
            self._cache['ewa_rows'] = self._file_cache['ewa_rows']
        num_valid_points, res = fornav(self._cache['ewa_cols'],
                                       self._cache['ewa_rows'],
                                       self.out_area, data,
                                       rows_per_scan=rows_per_scan)
        del num_valid_points

        return res

    def _project_array_bilinear(self, data):
        """Project array *data* using bilinear interpolation"""
        if 'bilinear_t' not in self._cache:
            self._cache['bilinear_t'] = self._file_cache['bilinear_t']
            self._cache['bilinear_s'] = self._file_cache['bilinear_s']
            self._cache['input_idxs'] = self._file_cache['input_idxs']
            self._cache['idx_arr'] = self._file_cache['idx_arr']

        res = get_sample_from_bil_info(data.ravel(),
                                       self._cache['bilinear_t'],
                                       self._cache['bilinear_s'],
                                       self._cache['input_idxs'],
                                       self._cache['idx_arr'],
                                       output_shape=self.out_area.shape)
        res = np.ma.masked_invalid(res)

        return res

    def project_array(self, data):
        """Project an array *data* along the given Projector object.
        """

        if self.mode == "nearest":
            res = self._project_array_nearest(data)

        elif self.mode == "quick":
            res = self._project_array_quick(data)

        elif self.mode == "ewa":
            if fornav is not None:
                res = self._project_array_ewa(data)
            else:
                raise ImportError("Can't import pyresample.ewa")

        elif self.mode == "bilinear":
            res = self._project_array_bilinear(data)

        return res


def calc_nearest_params(in_area, out_area, radius, nprocs=1):
    """Calculate projection parameters for nearest neighbour
    interpolation"""
    valid_index, valid_output_index, index_array, distance_array = \
        kd_tree.get_neighbour_info(in_area,
                                   out_area,
                                   radius,
                                   neighbours=1,
                                   nprocs=nprocs)
    del distance_array
    cache = {}
    cache['valid_index'] = valid_index
    cache['valid_output_index'] = valid_output_index
    cache['index_array'] = index_array

    return cache


def calc_quick_params(in_area, out_area):
    """Calculate projection parameters for quick interpolation mode"""
    ridx, cidx = utils.generate_quick_linesample_arrays(in_area,
                                                        out_area)
    cache = {}
    cache['row_idx'] = ridx
    cache['col_idx'] = cidx

    return cache


def calc_bilinear_params(in_area, out_area, radius, nprocs=1):
    """Calculate projection parameters for bilinear interpolation"""
    bilinear_t, bilinear_s, input_idxs, idx_arr = \
        get_bil_info(in_area, out_area, radius, neighbours=32,
                     nprocs=nprocs, masked=False)

    cache = {}
    cache['bilinear_s'] = bilinear_s
    cache['bilinear_t'] = bilinear_t
    cache['input_idxs'] = input_idxs
    cache['idx_arr'] = idx_arr

    return cache


def calc_ewa_params(in_area, out_area):
    """Calculate projection parameters for EWA interpolation"""
    swath_points_in_grid, cols, rows = ll2cr(in_area, out_area)
    del swath_points_in_grid
    cache = {}
    # cache['ewa_swath_points_in_grid'] = \
    #     swath_points_in_grid
    cache['ewa_cols'] = cols
    cache['ewa_rows'] = rows

    return cache


def get_precompute_cache_fname(in_id, out_id, in_area, out_area, mode,
                               proj_dir):
    """Form filename for precompute cache"""

    filename = (in_id + "2" + out_id + "_" +
                str(_get_area_hash(in_area)) + "to" +
                str(_get_area_hash(out_area)) + "_" +
                mode + ".npz")

    return os.path.join(proj_dir, filename)


def get_area_and_id(area, latlons=None):
    try:
        area_def = get_area_def(area)
        area_id = area
    except (utils.AreaNotFound, AttributeError):
        try:
            area_id = area.area_id
            area_def = area
        except AttributeError:
            if latlons is None:
                raise
            else:
                # TODO: Note that latlons are in order (lons, lats)
                area_def = geometry.SwathDefinition(lons=latlons[0],
                                                    lats=latlons[1])
                area_id = area

    return area_def, area_id
