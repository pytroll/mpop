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

"""Utility functions to read Grib messages
"""

import os
import pygrib
import os.path


class Grib(object):

    def __init__(self, fname):

        self._abspath = os.path.abspath(fname)

    @property
    def nmsgs(self):
        '''Number of GRIB messages in file.
        '''

        prop = 'nmsgs'
        attr = '_{}'.format(prop)

        if not hasattr(self, attr):
            grbs = pygrib.open(self._abspath)
            nmsgs = grbs.messages
            grbs.close()

            setattr(self, attr, nmsgs)

        return getattr(self, attr)

    def get(self, gmessage, key='values'):
        '''
        Returns the value for the 'key' for a given message number 'gmessage' or
        message field name 'gmessage'.
        '''

        grbs = pygrib.open(self._abspath)

        if type(gmessage) == int:
            mnbr = gmessage
        elif type(gmessage) == str:
            msg_found = False
            msgnum = 1
            while msgnum < self.nmsgs + 1:
                if grbs[msgnum]['parameterName'] == gmessage:
                    msg_found = True
                    break
                msgnum = msgnum + 1

            if msg_found:
                mnbr = msgnum
            else:
                print("No Grib message found with parameter name = %s" %
                      gmessage)
                return None

        if grbs[mnbr].valid_key(key):

            arr = grbs[mnbr][key]
            grbs.close()
            return arr
        else:
            grbs.close()
            return
