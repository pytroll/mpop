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

"""Unit testing the pps level-2 netCDF reader
"""
import unittest
from mpop.satin.nc_pps_l2 import PPSReader

CTTH_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_CTTH_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'
CT_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_CT_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'
CPP_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_CPP_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'
CMA_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_CMA_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'
PC_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_PC_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'
GEO_TESTFILE_LOCAL_1 = '/path/to/my/products/S_NWC_CMA_noaa19_33897_20150906T1240015Z_20150906T1240598Z.nc'


CTTH_TESTFILE_EARS_1 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CTTH_C_EUMS_20150914093100_15510.nc.bz2'
CT_TESTFILE_EARS_1 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CT_C_EUMS_20150914093100_15510.nc.bz2'
CMA_TESTFILE_EARS_1 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CMA_C_EUMS_20150914093100_15510.nc.bz2'
CTTH_TESTFILE_EARS_2 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CTTH_C_EUMS_20150914093200_15510.nc.bz2'
CT_TESTFILE_EARS_2 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CT_C_EUMS_20150914093200_15510.nc.bz2'
CMA_TESTFILE_EARS_2 = '/path/to/my/products/W_XX-EUMETSAT-Darmstadt,SING+LEV+SAT,METOPB+CMA_C_EUMS_20150914093200_15510.nc.bz2'


class TestPPSReader(unittest.TestCase):

    """Class for testing the PPSReader reader class.
    """

    def setUp(self):
        self.reader = PPSReader(object)

    def test_determine_prod_and_geo_files_local(self):
        """Test the private method _determine_prod_and_geo_files """
        self.reader._source = 'local'
        self.reader._cloud_product_geodir = None
        self.reader._geolocation_product_name = 'CMA'
        pro, geo = self.reader._determine_prod_and_geo_files(
            CTTH_TESTFILE_LOCAL_1)
        self.assertEqual(pro.keys(), ['CTTH'])
        self.assertEqual(pro['CTTH'], [CTTH_TESTFILE_LOCAL_1])
        self.assertEqual(geo.keys(), ['CTTH'])
        self.assertEqual(geo['CTTH'], [GEO_TESTFILE_LOCAL_1])

    def test_determine_prod_and_geo_files_ears(self):
        """Test the private method _determine_prod_and_geo_files """
        self.reader._source = 'ears'
        self.reader._cloud_product_geodir = None
        self.reader._geolocation_product_name = None
        pro, geo = self.reader._determine_prod_and_geo_files(
            CTTH_TESTFILE_EARS_1)
        self.assertEqual(pro.keys(), ['CTTH'])
        self.assertEqual(pro['CTTH'], [CTTH_TESTFILE_EARS_1])

        pro, geo = self.reader._determine_prod_and_geo_files(
            [CTTH_TESTFILE_EARS_1, CT_TESTFILE_EARS_1])
        nkeys = len(pro.keys())
        self.assertEqual(nkeys, 2)
        for key in pro.keys():
            self.assertTrue(key in ["CTTH", "CT"])

        for key in pro.keys():
            self.assertEqual(geo[key], pro[key])

    def tearDown(self):
        pass


def suite():
    """The test suite for test_viirs_sdr.
    """
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestPPSReader))

    return mysuite
