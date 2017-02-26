#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2009, 2012, 2013, 2014.

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

"""Test module for mpop.projector.
"""
import unittest
import numpy as np

from mock import MagicMock, patch
import sys
sys.modules['pyresample'] = MagicMock()
sys.modules['pyresample.bilinear'] = MagicMock()

from pyresample import geometry, utils

from mpop.projector import Projector
import mpop.projector


class TestProjector(unittest.TestCase):

    """Class for testing the Projector class.
    """

    proj = None

    @patch('mpop.projector.get_bil_info')
    @patch.object(utils, 'generate_quick_linesample_arrays')
    @patch.object(mpop.projector.kd_tree, 'get_neighbour_info')
    @patch.object(mpop.projector, '_get_area_hash')
    def test_init(self, gah, gni, gqla, bil_info):
        """Creation of coverage.
        """

        # in case of wrong number of arguments

        self.assertRaises(TypeError, Projector)
        self.assertRaises(TypeError, Projector, random_string(20))

        # in case of string arguments

        in_area_id = random_string(20)
        out_area_id = random_string(20)

        area_type = utils.parse_area_file.return_value.__getitem__.return_value

        gni.side_effect = [("a", "b", "c", "d")] * 10

        self.proj = Projector(in_area_id, out_area_id)
        self.assertEquals(utils.parse_area_file.call_count, 2)
        area_file = mpop.projector.get_area_file()
        utils.parse_area_file.assert_any_call(area_file, in_area_id)
        utils.parse_area_file.assert_any_call(area_file, out_area_id)

        self.assertEquals(self.proj.in_area, area_type)
        self.assertEquals(self.proj.out_area, area_type)

        # in case of undefined areas

        mock = MagicMock(side_effect=Exception("raise"))
        with patch.object(utils, 'parse_area_file', mock):
            self.assertRaises(Exception,
                              Projector,
                              "raise",
                              random_string(20))
            self.assertRaises(Exception,
                              Projector,
                              random_string(20),
                              "raise")

        # in case of geometry objects as input

        with patch.object(utils, 'AreaNotFound', Exception):
            mock = MagicMock(side_effect=[utils.AreaNotFound("raise"),
                                          MagicMock()])
            with patch.object(utils, 'parse_area_file', mock):
                in_area = geometry.AreaDefinition()
                self.proj = Projector(in_area, out_area_id)
                self.assertEquals(self.proj.in_area, in_area)

        in_area = geometry.SwathDefinition()
        utils.parse_area_file.return_value.__getitem__.side_effect = [
            AttributeError, out_area_id]
        self.proj = Projector(in_area, out_area_id)
        self.assertEquals(self.proj.in_area, in_area)

        out_area = geometry.AreaDefinition()
        utils.parse_area_file.return_value.__getitem__.side_effect = [
            in_area_id, AttributeError]
        self.proj = Projector(in_area_id, out_area)
        self.assertEquals(self.proj.out_area, out_area)

        # in case of lon/lat is input

        utils.parse_area_file.return_value.__getitem__.side_effect = [
            AttributeError, out_area_id]
        lonlats = ("great_lons", "even_greater_lats")

        self.proj = Projector("raise", out_area_id, lonlats)
        geometry.SwathDefinition.assert_called_with(lons=lonlats[0],
                                                    lats=lonlats[1])

        utils.parse_area_file.return_value.__getitem__.side_effect = None
        # in case of wrong mode

        self.assertRaises(ValueError,
                          Projector,
                          random_string(20),
                          random_string(20),
                          mode=random_string(20))

        utils.parse_area_file.return_value.__getitem__.side_effect = ["a", "b",
                                                                      "c", "d"]
        gqla.side_effect = [("ridx", "cidx")]
        # quick mode cache
        self.proj = Projector(in_area_id, out_area_id, mode="quick")
        cache = getattr(self.proj, "_cache")
        self.assertTrue(cache['row_idx'] is not None)
        self.assertTrue(cache['col_idx'] is not None)

        # nearest mode cache
        self.proj = Projector(in_area_id, out_area_id, mode="nearest")
        cache = getattr(self.proj, "_cache")
        self.assertTrue(cache['valid_index'] is not None)
        self.assertTrue(cache['valid_output_index'] is not None)
        self.assertTrue(cache['index_array'] is not None)

        # bilinear mode cache
        bil_info.return_value = (1, 2, 3, 4)

        def spam(val):
            return 'adef'

        with patch.object(mpop.projector, 'get_area_def', spam):
            self.proj = Projector(in_area_id, out_area_id, mode="bilinear")
        cache = getattr(self.proj, "_cache")
        self.assertTrue(cache['bilinear_t'] is not None)
        self.assertTrue(cache['bilinear_s'] is not None)
        self.assertTrue(cache['input_idxs'] is not None)
        self.assertTrue(cache['idx_arr'] is not None)

    @patch.object(np.ma, "array")
    @patch.object(mpop.projector.kd_tree, 'get_sample_from_neighbour_info')
    @patch.object(np, "load")
    def test_project_array(self, npload, gsfni, marray):
        """Test the project_array function.
        """
        in_area_id = random_string(20)
        out_area_id = random_string(20)
        data = np.random.standard_normal((3, 1))

        utils.parse_area_file.return_value.__getitem__.side_effect = [
            "a", "b", "c", "d"]
        # test quick
        self.proj = Projector(in_area_id, out_area_id, mode="quick")
        self.proj.project_array(data)
        mpop.projector.image.ImageContainer.assert_called_with(
            data, "a", fill_value=None)
        mpop.projector.image.ImageContainer.return_value.\
            get_array_from_linesample.assert_called_with(
                self.proj._cache["row_idx"], self.proj._cache["col_idx"])
        marray.assert_called_once_with(
            mpop.projector.image.ImageContainer.return_value.
            get_array_from_linesample.return_value,
            dtype=np.dtype('float64'))

        # test nearest
        in_area = MagicMock()
        out_area = MagicMock()
        utils.parse_area_file.return_value.__getitem__.side_effect = [
            in_area, out_area]
        self.proj = Projector(in_area_id, out_area_id, mode="nearest")
        self.proj.project_array(data)
        mpop.projector.kd_tree.get_sample_from_neighbour_info.\
            assert_called_with('nn',
                               out_area.shape,
                               data,
                               npload.return_value.__getitem__.return_value,
                               npload.return_value.__getitem__.return_value,
                               npload.return_value.__getitem__.return_value,
                               fill_value=None)

    @patch.object(mpop.projector.kd_tree, 'get_neighbour_info')
    def test_calc_nearest_params(self, gni):
        gni.return_value = (1, 2, 3, 4)
        res = mpop.projector.calc_nearest_params('in_area', 'out_area',
                                                 'radius', nprocs='nprocs')
        self.assertTrue(isinstance(res, dict))
        self.assertTrue('valid_index' in res)
        self.assertEqual(res['valid_index'], 1)
        self.assertTrue('valid_output_index' in res)
        self.assertEqual(res['valid_output_index'], 2)
        self.assertTrue('index_array' in res)
        self.assertEqual(res['index_array'], 3)

    @patch.object(mpop.projector.utils, 'generate_quick_linesample_arrays')
    def test_calc_quick_params(self, gqla):
        gqla.return_value = (1, 2)
        res = mpop.projector.calc_quick_params('in_area', 'out_area')

        self.assertTrue(isinstance(res, dict))
        self.assertTrue('row_idx' in res)
        self.assertEqual(res['row_idx'], 1)
        self.assertTrue('col_idx' in res)
        self.assertEqual(res['col_idx'], 2)

    @patch.object(mpop.projector, 'get_bil_info')
    def test_calc_bilinear_params(self, gbi):
        gbi.return_value = (1, 2, 3, 4)
        res = mpop.projector.calc_bilinear_params('in_area', 'out_area',
                                                  'radius', nprocs='nprocs')
        self.assertTrue(isinstance(res, dict))
        self.assertTrue('bilinear_t' in res)
        self.assertEqual(res['bilinear_t'], 1)
        self.assertTrue('bilinear_s' in res)
        self.assertEqual(res['bilinear_s'], 2)
        self.assertTrue('input_idxs' in res)
        self.assertEqual(res['input_idxs'], 3)
        self.assertTrue('idx_arr' in res)
        self.assertEqual(res['idx_arr'], 4)

    @patch.object(mpop.projector, 'll2cr')
    def test_calc_ewa_params(self, ll2):
        ll2.return_value = (0, 1, 2)
        res = mpop.projector.calc_ewa_params('in_area', 'out_area')
        self.assertTrue(isinstance(res, dict))
        self.assertTrue('ewa_cols' in res)
        self.assertEqual(res['ewa_cols'], 1)
        self.assertTrue('ewa_rows' in res)
        self.assertEqual(res['ewa_rows'], 2)

    def test_get_precompute_cache_fname(self):
        res = mpop.projector.get_precompute_cache_fname('in_id', 'out_id',
                                                        'in_area', 'out_area',
                                                        'mode', 'proj_dir')
        cor_res = "proj_dir/in_id2out_id_-" + \
                  "6296787761359943868to8984161303220364208_mode.npz"
        self.assertTrue(res == cor_res)

    @patch.object(mpop.projector, 'get_area_def')
    @patch.object(mpop.projector.geometry, 'SwathDefinition')
    def test_get_area_and_id(self, swath_def, gad):
        # Case when get_area_def works
        swath_def.return_value = 1
        gad.return_value = 'adef'
        res = mpop.projector.get_area_and_id('area')
        self.assertTrue(res[0] == 'adef')
        self.assertTrue(res[1] == 'area')
        # Case when AttributeError is raised
        with self.assertRaises(AttributeError):
            gad.side_effect = AttributeError
            res = mpop.projector.get_area_and_id('area')
        # Case when AttributeError is raised and latlons are given
        gad.side_effect = AttributeError
        res = mpop.projector.get_area_and_id('area', latlons=[1, 2])
        self.assertEqual(res[0], 1)
        self.assertTrue(res[1], 'area')


def random_string(length,
                  choices="abcdefghijklmnopqrstuvwxyz"
                  "ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
    """Generates a random string with elements from *set* of the specified
    *length*.
    """
    import random
    return "".join([random.choice(choices)
                    for dummy in range(length)])


def suite():
    """The test suite for test_projector.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestProjector))

    return mysuite
