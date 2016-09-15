#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2009, 2013, 2015, 2016.

# SMHI,
# Folkborgsvägen 1,
# Norrköping,
# Sweden

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This file is part of the mpop.

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

"""Palette holder module.
"""

import numpy as np


def tv_legend():
    """Palette for TV.
    """
    legend = []
    legend.append((0,   0,   0))  # Unprocessed: Black
    legend.append((0, 120,   0))  # Land
    legend.append((0,   0, 215))  # Sea: Blue
    legend.append((0, 120,   0))  # Land (Snow on land)
    legend.append((0,   0, 215))  # Sea: Blue (Snow/Ice on sea)

    for i in range(5, 256):
        # All other pixel values are grey according to IR temp.
        legend.append((i, i, i))

    return convert_palette(legend)


def vv_legend():
    """Palette for Swedish road authorities (Vägverket).
    """
    legend = []
    legend.append((0,   0,   0))  # Unprocessed: Black
    legend.append((0, 120,   0))  # Land
    legend.append((0,   0, 215))  # Sea: Blue
    # Cloud type values 5 to 8:
    legend.append((255, 150,   0))  # Very low cumuliform
    legend.append((255, 100,   0))  # Very low
    legend.append((255, 220,   0))  # Low cumuliform
    legend.append((255, 180,   0))  # Low

    for i in range(7, 256):
        # All other pixel values are grey according to IR temp.
        legend.append((i, i, i))

    return convert_palette(legend)


def cloud_phase():
    """Palette for cloud phase.
    """
    legend = []
    legend.append((0,    0,   0))  # Unprocessed: Black
    legend.append((0,    0, 215))  # Water Clouds: Blue
    legend.append((240, 240, 240))  # Ice Clouds: Almost White
    legend.append((120, 120,   0))  # Uncertain Phase: ?

    return convert_palette(legend)


def cms_modified():
    """Palette for regular cloud classification.
    """
    return nwcsaf_cloudtype()


def nwcsaf_cloudtype():
    """Palette for regular cloud classification.
    """
    legend = []
    legend.append((100, 100, 100))  # Unprocessed: Grey
    legend.append((0, 120,   0))
    legend.append((0,   0,   0))  # Sea: Black
    legend.append((250, 190, 250))  # Snow
    legend.append((220, 160, 220))  # Sea-ice

    legend.append((255, 150,   0))  # Very low cumuliform
    legend.append((255, 100,   0))  # Very low
    legend.append((255, 220,   0))  # Low cumuliform
    legend.append((255, 180,   0))  # Low
    legend.append((255, 255, 140))  # Medium cumuliform
    legend.append((240, 240,   0))  # Medium
    legend.append((250, 240, 200))  # High cumiliform
    legend.append((215, 215, 150))  # High
    legend.append((255, 255, 255))  # Very high cumuliform
    legend.append((230, 230, 230))  # Very high

    legend.append((0,  80, 215))  # Semi-transparent thin
    legend.append((0, 180, 230))  # Semi-transparent medium
    legend.append((0, 240, 240))  # Semi-transparent thick
    legend.append((90, 200, 160))  # Semi-transparent above
    legend.append((200,   0, 200))  # Broken
    legend.append((95,  60,  30))  # Undefined: Brown

    return convert_palette(legend)


def ctth_height():
    """CTTH height palette.
    """
    legend = []
    legend.append((0,   0,   0))
    legend.append((255,   0, 216))  # 0 meters
    legend.append((126,   0,  43))
    legend.append((153,  20,  47))
    legend.append((178,  51,   0))
    legend.append((255,  76,   0))
    legend.append((255, 102,   0))
    legend.append((255, 164,   0))
    legend.append((255, 216,   0))
    legend.append((216, 255,   0))
    legend.append((178, 255,   0))
    legend.append((153, 255,   0))
    legend.append((0, 255,   0))
    legend.append((0, 140,  48))
    legend.append((0, 178, 255))
    legend.append((0, 216, 255))
    legend.append((0, 255, 255))
    legend.append((238, 214, 210))
    legend.append((239, 239, 223))
    legend.append((255, 255, 255))  # 10,000 meters
    for i in range(79):
        legend.append((255, 255, 255))
    legend.append((224, 224, 224))

    return convert_palette(legend)


def ctth_height_pps():
    """CTTH height palette for NWCSAF/PPS.
    Identical to the one found in the hdf5 files.
    """
    legend = []
    legend.append((255,   0, 216))  # 0 meters
    legend.append((255,   0, 216))  # 0 meters
    legend.append((255,   0, 216))  # 0 meters
    legend.append((126,   0,  43))
    legend.append((126,   0,  43))
    legend.append((153,  20,  47))
    legend.append((153,  20,  47))
    legend.append((153,  20,  47))
    legend.append((178,  51,   0))
    legend.append((178,  51,   0))
    legend.append((255,  76,   0))
    legend.append((255,  76,   0))
    legend.append((255,  76,   0))
    legend.append((255, 102,   0))
    legend.append((255, 102,   0))
    legend.append((255, 164,   0))
    legend.append((255, 164,   0))
    legend.append((255, 164,   0))
    legend.append((255, 216,   0))
    legend.append((255, 216,   0))
    legend.append((216, 255,   0))
    legend.append((216, 255,   0))
    legend.append((178, 255,   0))
    legend.append((178, 255,   0))
    legend.append((178, 255,   0))
    legend.append((153, 255,   0))
    legend.append((153, 255,   0))
    legend.append((0, 255,   0))
    legend.append((0, 255,   0))
    legend.append((0, 255,   0))
    legend.append((0, 140,  48))
    legend.append((0, 140,  48))
    legend.append((0, 178, 255))
    legend.append((0, 178, 255))
    legend.append((0, 178, 255))
    legend.append((0, 216, 255))
    legend.append((0, 216, 255))
    legend.append((0, 255, 255))
    legend.append((0, 255, 255))
    legend.append((0, 255, 255))
    legend.append((238, 214, 210))
    legend.append((238, 214, 210))
    legend.append((239, 239, 223))
    legend.append((239, 239, 223))
    for idx in range(47, 150):
        legend.append((255, 255, 255))  # 10,000 meters

    for idx in range(150, 256):
        legend.append((0, 0, 0))

    return convert_palette(legend)


def chlorophyll_a():
    """Chlorophyll-A legend for MERIS"""
    raise NotImplementedError("This palette is not implemented - " +
                              "it was earlier though...")


# --------------------------------------------
#   Define colour palette LUT for SST palette image
#   colour shading; blue(cold) -> green -> yellow -> red(warm)
#


def sstlut_osisaf_metno():

    legend = []
    legend.append((0, 0, 0))  # /*  Black (out of image)  */
    legend.append((82, 82, 82))  # /*  Dark grey (land)  */
    legend.append((187, 187, 187))  # /* Light grey (cloud contaminated)  */
    legend.append((255, 255, 255))  # /* White (sea ice and snow)  */

    legend.append((255, 0, 255))  # /* Starting at 4 = pink */
    legend.append((195, 0, 129))  # /* Dark pink */
    legend.append((129, 0, 47))  # /* Dark red */
    legend.append((195, 0, 0))  # /* Medium dark red */

    r = 255  # /* Red */
    g = 0
    b = 0
    for i in range(8, 11, 1):
        legend.append((r, g, b))
        r = r - 19
        g = g + 43

    r = 200  # /* Brown */
    g = 128
    b = 0
    for i in range(11, 16, 1):
        legend.append((r, g, b))
        r = r + 11
        g = g + 26
        b = b + 13

    r = 65535 / 256  # /* Yellow */
    g = 65535 / 256
    b = 16185 / 256
    legend.append((r, g, b))
    legend.append((52000 / 256, 65535 / 256, 13500 / 256))
    legend.append((35000 / 256, 65535 / 256, 7000 / 256))

    r = 0  # /* Green */
    g = 65535 / 256
    b = 0
    for i in range(19, 22, 1):
        legend.append((r, g, b))
        g = g - 12422 / 256

    r = 0  # /* Dark Green */
    g = 28269 / 256
    b = 0
    for i in range(22, 26, 1):
        legend.append((r, g, b))
        g = g - 7067 / 256
        b = b + 16383 / 256

    r = 0
    g = 0
    b = 65535 / 256
    legend.append((r, g, b))  # Blue
    legend.append((25700 / 256, 25700 / 256, 65535 / 256))  # Dark purple
    legend.append((38550 / 256, 38550 / 256, 65535 / 256))  # Light purple

    for i in range(29, 256):
        legend.append((0, 0, 0))

    return convert_palette(legend)


def convert_palette(palette):
    """Convert palette from [0,255] range to [0,1].
    """
    new_palette = []
    for i in palette:
        new_palette.append((i[0] / 255.0,
                            i[1] / 255.0,
                            i[2] / 255.0))
    return new_palette


def convert_palette2colormap(palette):
    """Convert palette from [0,255] range to [0,1].
    """
    from trollimage.colormap import Colormap
    j = 0
    n_pal = len(palette) - 1
    values = []
    colors = []

    red = [r for (r, g, b) in palette]
    green = [g for (r, g, b) in palette]
    blue = [b for (r, g, b) in palette]

    max_pal = max(max(red), max(blue), max(green))
    if max_pal <= 1.0:
        # print "palette already normalized"
        denom = 1.0
    else:
        # print "palette normalized to 255"
        denom = 255.0

    for i in palette:
        values.append((n_pal - j) / float(n_pal))
        colors.append((i[0] / denom, i[1] / denom, i[2] / denom))
        j = j + 1
    # reverse order to the entries
    values = values[::-1]
    colors = colors[::-1]

    # for i in palette:
    #    values.append( j /  float(n_pal))
    #    colors.append((i[0] / 255.0, i[1] / 255.0, i[2] / 255.0))
    #    j=j+1

    # attention:
    # Colormap(values, colors) uses the second input option of Colormap
    # values has to be a list (not a tuple) and
    # colors has to be the corresponding list of color tuples

    return Colormap(values, colors)


class LogColors(object):

    """
    Defines colors to use with `logdata2image`

    """

    def __init__(self, nodata, zeros, over, breaks):
        self.nodata = nodata
        self.zeros = zeros
        self.over = over
        self.breaks = breaks

    def palette(self, N=256):
        """
        Build a palette for logarithmic data images.

        """

        max_value = self.breaks[-1][0]

        pal = np.zeros((N, 3), dtype=np.uint8)

        b_last, rgb_last = self.breaks[0]
        for b, rgb in self.breaks[1:]:
            # Get a slice of the palette array for the current interval
            p = pal[
                np.log(b_last + 1) * N / np.log(max_value):np.log(b + 1) * N / np.log(max_value)]
            for i in range(3):  # red, green, blue
                p[:, i] = np.linspace(rgb_last[i], rgb[i], p.shape[0])
            b_last = b
            rgb_last = rgb

        pal[0] = self.nodata
        pal[1] = self.zeros
        pal[-1] = self.over

        return pal


class TriColors(LogColors):

    """
    Use three color tones in the intervals between the elements of *breaks*.

    """
    color_tones = [((0, 0, 200), (150, 150, 255)),  # dark to light blue
                   ((150, 150, 0), (255, 255, 8)),  # greyish to bright yellow
                   ((230, 150, 100), (230, 0, 0))]  # green to red

    nodata = (0, 0, 0)  # black
    # zeros = (20, 0, 20) # dark purple
    # black  #There is no need to mark zeros with another col
    zeros = (0, 0, 0)
    over = (255, 0, 0)  # bright red

    def __init__(self, breaks):
        breaks = [(breaks[0], TriColors.color_tones[0][0]),
                  (breaks[1], TriColors.color_tones[0][1]),

                  (breaks[1], TriColors.color_tones[1][0]),
                  (breaks[2], TriColors.color_tones[1][1]),

                  (breaks[2], TriColors.color_tones[2][0]),
                  (breaks[3], TriColors.color_tones[2][1])]

        LogColors.__init__(self, TriColors.nodata, TriColors.zeros,
                           TriColors.over, breaks)

CPP_COLORS = {'cpp_cot': TriColors([0, 3.6, 23, 700]),  # ISCCP intervals
              'cpp_reff': TriColors([0, 10, 20, 1000])}

CPP_COLORS['cot'] = CPP_COLORS['cpp_cot']
CPP_COLORS['reff'] = CPP_COLORS['cpp_reff']


def get_ctp_legend():
    """
    Get the Cloud Top Pressure color palette
    """

    legend = []
    legend.append((0, 0, 0))     # No data
    legend.append((255, 0, 216))  # 0: 1000-1050 hPa (=100000-105000 Pa)
    legend.append((126, 0, 43))  # 1: 950-1000 hPa
    legend.append((153, 20, 47))  # 2: 900-950 hPa
    legend.append((178, 51, 0))  # 3: 850-900 hPa
    legend.append((255, 76, 0))  # 4: 800-850 hPa
    legend.append((255, 102, 0))  # 5: 750-800 hPa
    legend.append((255, 164, 0))  # 6: 700-750 hPa
    legend.append((255, 216, 0))  # 7: 650-700 hPa
    legend.append((216, 255, 0))  # 8: 600-650 hPa
    legend.append((178, 255, 0))  # 9: 550-600 hPa
    legend.append((153, 255, 0))  # 10: 500-550 hPa
    legend.append((0, 255, 0))   # 11: 450-500 hPa
    legend.append((0, 140, 48))  # 12: 400-450 hPa
    legend.append((0, 178, 255))  # 13: 350-400 hPa
    legend.append((0, 216, 255))  # 14: 300-350 hPa
    legend.append((0, 255, 255))  # 15: 250-300 hPa
    legend.append((238, 214, 210))  # 16: 200-250 hPa
    legend.append((239, 239, 223))  # 17: 150-200 hPa
    legend.append((255, 255, 255))  # 18: 100-150 hPa
    legend.append((255, 255, 255))  # 19: 50-100 hPa
    legend.append((255, 255, 255))  # 20: 0-50 hPa  (=0-5000 Pa)

    palette = convert_palette(legend)
    return palette


def get_reff_legend():
    return get_log_legend('reff')


def get_cot_legend():
    return get_log_legend('cot')


def get_log_legend(product_name):
    # This is the same data as is used in logdata2image (when indata as for
    # the calls from cppimage)
    legend = CPP_COLORS[product_name].palette()
    palette = convert_palette(legend)
    return palette


def oca_get_scenetype_legend():

    # Colorize using PPS/CPP palette
    legend = np.array([[170, 130, 255],  # purple/blue for liquid (cph == 1)
                       [220, 200, 255],  # almost white for ice (cph == 2)
                       [255, 200, 200]   # Redish for multi layer clouds
                       ])
    legend = np.vstack([np.zeros((111, 3)), legend])
    palette = convert_palette(legend)
    return palette
