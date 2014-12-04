#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2014 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <adam.dybbroe@smhi.se>
#   Northaholic <northaholic@icloud.com>

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

"""Interface to Eumetcast level 1.5 Nativeformat. Uses the MIPP reader.
"""

import logging

LOG = logging.getLogger(__name__)

from mpop.plugin_base import Reader


class NativeReader(Reader):

    pformat = "mipp_native"

    def load(self, satscene, calibrate=1, filename=None, **kwargs):
        """Load the data"""

        from mipp import native
        image = native.MSG.NativeImage('meteosat10', filename=filename)

        for channel in satscene.channels_to_load:
            print channel
            chobj = getattr(image, channel.lower())
            satscene[channel] = chobj.data
