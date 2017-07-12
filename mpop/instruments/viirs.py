#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017.

# Author(s):

#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This file is part of mpop.

# mpop is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# mpop is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# mpop.  If not, see <http://www.gnu.org/licenses/>.

"""This modules describes the viirs instrument.
It provides VIIRS specific methods for RGB-compositing.
"""
import logging

import numpy as np

import mpop.imageo.geo_image as geo_image
from mpop.instruments.visir import VisirCompositer

LOG = logging.getLogger(__name__)

try:
    from pyorbital.astronomy import sun_zenith_angle as sza
except ImportError:
    LOG.warning("Sun zenith angle correction not possible! " +
                "Check the availability of the pyorbital module in your environment")
    sza = None


# VIIRS
# Since there is overlap between I-bands and M-bands we need to
# specifically re-define some of the RGB composites already defined
# in the standard visir.py module. So, the same names, like "overview"
# can be invoked and based on M-bands only.
# In addition we define new composite names for the I-bands,
# like e.g. hr_overview, hr_night_fog, etc
#


class ViirsCompositer(VisirCompositer):

    """This class sets up the VIIRS instrument channel list.
    """

    instrument_name = "viirs"

    def overview(self, stretch='linear', gamma=1.6):
        """Make an Overview RGB image composite from VIIRS
        channels.
        """
        self.check_channels('M05', 'M07', 'M15')

        ch1 = self['M05'].check_range()
        ch2 = self['M07'].check_range()
        ch3 = -self['M15'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")
        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    overview.prerequisites = set(['M05', 'M07', 'M15'])

    def overview_sun(self, stretch='linear', gamma=1.6, fill_value=(0, 0, 0)):
        """Make an overview RGB image composite normalising with cosine to the
        sun zenith angle.
        """
        self.check_channels('M05', 'M07', 'M15')

        lonlats = self['M15'].area.get_lonlats()

        red = self['M05'].sunzen_corr(self.time_slot, lonlats, limit=88.,
                                      sunmask=95).data
        green = self['M07'].sunzen_corr(self.time_slot, lonlats, limit=88.,
                                        sunmask=95).data
        blue = -self['M15'].data

        img = geo_image.GeoImage((red, green, blue),
                                 self.area,
                                 self.time_slot,
                                 fill_value=fill_value,
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    overview_sun.prerequisites = set(['M05', 'M07', 'M15'])

    def hr_overview(self):
        """Make a high resolution Overview RGB image composite 
        from the VIIRS I-bands only - 375 meter resolution.
        """
        self.check_channels('I01', 'I02', 'I05')

        ch1 = self['I01'].check_range()
        ch2 = self['I02'].check_range()
        ch3 = -self['I05'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch="crude")
        img.enhance(gamma=1.6)

        return img

    hr_overview.prerequisites = set(['I01', 'I02', 'I05'])

    def truecolor(self, stretch='linear', gamma=2.0):
        """Make a True Color RGB image composite from
        M-bands only.
        """
        #elf.check_channels('M02', 'M04', 'M05')
        self.check_channels('M03', 'M04', 'M05')

        ch1 = self['M05'].check_range()
        ch2 = self['M04'].check_range()
        #h3 = self['M02'].check_range()
        ch3 = self['M03'].check_range()

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=None,
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        if gamma:
            img.enhance(gamma=gamma)

        return img

    truecolor.prerequisites = set(['M03', 'M04', 'M05'])
#   truecolor.prerequisites = set(['M02', 'M04', 'M05'])

    def natural(self):
        """Make a Natural Colors RGB image composite from
        M-bands only.
        """
        self.check_channels('M05', 'M06', 'M07', 'M10')

        ch1 = self['M10'].check_range()
        ch2 = self['M07'].check_range()
        ch3 = self['M05'].check_range()

        ch2b = self['M06'].check_range()
        ch2 = np.ma.where(ch2.mask, ch2b, ch2)

        common_mask = np.logical_or(ch1.mask, ch2.mask)
        common_mask = np.logical_or(common_mask, ch3.mask)
        ch1.mask = common_mask
        ch2.mask = common_mask
        ch3.mask = common_mask

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((0, 90),
                                         (0, 90),
                                         (0, 90)))

        img.enhance(gamma=1.8)

        return img

    natural.prerequisites = set(['M05', 'M06', 'M07', 'M10'])

    def hr_natural(self):
        """Make a high resolution Day Natural Colors RGB image 
        composite from I-bands only - 375 meter resolution.
        """
        self.check_channels('I01', 'I02', 'I03')

        ch1 = self['I03'].check_range()
        ch2 = self['I02'].check_range()
        ch3 = self['I01'].check_range()

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((0, 90),
                                         (0, 90),
                                         (0, 90)))

        img.enhance(gamma=1.8)

        return img

    hr_natural.prerequisites = set(['I01', 'I02', 'I03'])

    def vis06(self):
        """Make a black and white image of the VIS 0.635um channel.
        """
        return self.channel_image("M05")

    vis06.prerequisites = set(['M05'])

    def hr_vis06(self):
        """Make a high res black and white image of the 
        'visible' (VIS) I-band at 0.640um.
        """
        return self.channel_image('I01')

    hr_vis06.prerequisites = set(['I01'])

    def green_snow(self):
        """Make a Green Snow RGB image composite.
        """
        self.check_channels('M05', 'M10', 'M15')

        ch1 = self['M10'].check_range()
        ch2 = self['M05'].check_range()
        ch3 = -self['M15'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch="crude")
        img.enhance(gamma=1.6)

        return img

    green_snow.prerequisites = set(['M05', 'M10', 'M15'])

    def hr_green_snow(self):
        """Make a Green Snow RGB image composite.
        """
        self.check_channels('I01', 'I03', 'I05')

        ch1 = self['I03'].check_range()
        ch2 = self['I01'].check_range()
        ch3 = -self['I05'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch="crude")
        img.enhance(gamma=1.6)

        return img

    hr_green_snow.prerequisites = set(['I01', 'I03', 'I05'])

    def red_snow(self):
        """Make a Red Snow RGB image composite.
        """
        self.check_channels('M05', 'M10', 'M15')

        ch1 = self['M05'].check_range()
        ch2 = self['M10'].check_range()
        ch3 = -self['M15'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch="crude")

        return img

    red_snow.prerequisites = set(['M05', 'M10', 'M15'])

    def hr_red_snow(self):
        """Make a high resolution Red Snow RGB image composite
        from the I-bands only.
        """
        self.check_channels('I01', 'I03', 'I05')

        ch1 = self['I01'].check_range()
        ch2 = self['I03'].check_range()
        ch3 = -self['I05'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch="crude")

        return img

    hr_red_snow.prerequisites = set(['I01', 'I03', 'I05'])

    def dnb_overview(self, stretch='linear'):
        """Make a nighttime overview RGB image composite from VIIRS
        DNB and M bands.
        """
        self.check_channels('DNB', 'M15')

        lonlats = self['M15'].area.get_lonlats()

        if sza:
            sunz = sza(self.time_slot, lonlats[0], lonlats[1])
            sunz = np.ma.masked_outside(sunz, 103, 180)
            sunzmask = sunz.mask

            red = np.ma.masked_where(sunzmask, self['DNB'].data)
            green = np.ma.masked_where(sunzmask, self['DNB'].data)
            blue = np.ma.masked_where(sunzmask, -self['M15'].data)
        else:
            LOG.warning("No masking of solar contaminated pixels performed!")
            red = self['DNB'].data
            green = self['DNB'].data
            blue = -self['M15'].data

        img = geo_image.GeoImage((red, green, blue),
                                 self.area,
                                 self.time_slot,
                                 fill_value=None,
                                 mode="RGB")

        img.enhance(stretch=stretch)

        return img

    dnb_overview.prerequisites = set(['DNB', 'M15'])

    # def dnb_overview(self):
    #     """Make an Overview RGB image composite from VIIRS
    #     channels.
    #     """
    #     self.check_channels('DNB', 'M15')

    #     ch1 = self['DNB'].data
    #     ch2 = self['DNB'].data
    #     ch3 = -self['M15'].data

    #     img = geo_image.GeoImage((ch1, ch2, ch3),
    #                              self.area,
    #                              self.time_slot,
    #                              fill_value=None,
    #                              mode="RGB")

    #     img.enhance(stretch="linear")

    #     return img

    # dnb_overview.prerequisites = set(['DNB', 'M15'])

    def night_color(self):
        """Make a Night Overview RGB image composite.
        Same as cloudtop ... just different.
        """
        return self.cloudtop(stretch="histogram")

    night_color.prerequisites = set(['M12', 'M15', 'M16'])

    def night_fog(self):
        """Make a Night Fog RGB image composite.
        """
        self.check_channels('M12', 'M15', 'M16')

        ch1 = self['M16'].data - self['M15'].data
        ch2 = self['M15'].data - self['M12'].data
        ch3 = self['M15'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 6),
                                         (243, 293)))

        img.enhance(gamma=(1.0, 2.0, 1.0))

        return img

    night_fog.prerequisites = set(['M12', 'M15', 'M16'])

    def dust(self):
        """Make a Dust RGB image composite.
        """
        self.check_channels('M14', 'M15', 'M16')

        ch1 = self['M16'].data - self['M15'].data
        ch2 = self['M15'].data - self['M14'].data
        ch3 = self['M15'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 15),
                                         (261, 289)))

        img.enhance(gamma=(1.0, 2.5, 1.0))

        return img

    dust.prerequisites = set(['M14', 'M15', 'M16'])

    def ash(self):
        """Make a Ash RGB image composite.
        """
        self.check_channels('M14', 'M15', 'M16')

        ch1 = self['M16'].data - self['M15'].data
        ch2 = self['M15'].data - self['M14'].data
        ch3 = self['M15'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (-4, 5),
                                         (243, 303)))

        return img

    ash.prerequisites = set(['M14', 'M15', 'M16'])

    def fog(self):
        """Make a Fog RGB image composite.
        """
        self.check_channels('M14', 'M15', 'M16')

        ch1 = self['M16'].data - self['M15'].data
        ch2 = self['M15'].data - self['M14'].data
        ch3 = self['M15'].data
        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB",
                                 crange=((-4, 2),
                                         (0, 6),
                                         (243, 283)))

        img.enhance(gamma=(1.0, 2.0, 1.0))

        return img

    fog.prerequisites = set(['M14', 'M15', 'M16'])

    def cloudtop(self, stretch=None):
        """Make a Cloudtop RGB image composite.
        """
        self.check_channels('M12', 'M15', 'M16')

        ch1 = -self['M12'].data
        ch2 = -self['M15'].data
        ch3 = -self['M16'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        if stretch:
            img.enhance(stretch=stretch)
        else:
            img.enhance(stretch=(0.005, 0.005))

        return img

    cloudtop.prerequisites = set(['M12', 'M15', 'M16'])

    def dnb(self, stretch="histogram"):
        """Make a black and white image of the Day-Night band."""
        self.check_channels('DNB')

        img = geo_image.GeoImage(self['DNB'].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L")
        if stretch:
            img.enhance(stretch=stretch)
        return img

    dnb.prerequisites = set(['DNB'])

    def dnb_rgb(self, stretch="linear"):
        """Make a RGB Day-Night band using M15 as blue."""
        self.check_channels('DNB', 'M15')
        ch1 = self['DNB'].data
        ch2 = self['DNB'].data
        ch3 = -self['M15'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")
        if stretch:
            img.enhance(stretch=stretch)
        return img

    dnb_rgb.prerequisites = set(['DNB', 'M15'])

    def ir108(self):
        """Make a black and white image of the IR 10.8um channel.
        """
        self.check_channels("M15")

        img = geo_image.GeoImage(self["M15"].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L",
                                 crange=(-70 + 273.15, 57.5 + 273.15))
        img.enhance(inverse=True)
        return img

    ir108.prerequisites = set(["M15"])

    def hr_ir108(self):
        """Make a black and white image of the IR 10.8um channel (320m).
        """
        self.check_channels("I05")

        img = geo_image.GeoImage(self["I05"].data,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L",
                                 crange=(-70 + 273.15, 57.5 + 273.15))
        img.enhance(inverse=True)
        return img

    hr_ir108.prerequisites = set(["I05"])

    def chlorophyll(self, stretch=None):
        """ From http://oceancolor.gsfc.nasa.gov/REPROCESSING/R2009/ocv6/

        * Rrs1 = blue wavelength Rrs (e.g., 443, 490, or 510-nm)
        * Rrs2 = green wavelength Rrs (e.g., 547, 555, or 565-nm)
        * X = log10(Rrs1 / Rrs2)
        * chlor_a = 10^(a0 + a1*X + a2*X^2 + a3*X^3 + a4*X^4)

        sensor  default *      blue     green   a0       a1      a2       a3       a4
        OC3V    VIIRS   Y      443>486  550     0.2228  -2.4683  1.5867  -0.4275  -0.7768

        blue: M02(445)>M03(488)
        green: M04(555)

        * X = log10(max(M2, M3)/M4)
        """
        self.check_channels("M02", "M03", "M04")

        a0, a1, a2, a3, a4 = (0.2228, -2.4683, 1.5867, -0.4275, -0.7768)

        # X = np.maximum(self["M02"].data, self["M03"].data)/self["M04"].data
        X = self["M02"].data / self["M04"].data
        X = np.log10(X)
        chlor_a = 10 ** (a0 + a1 * X + a2 * (X ** 2) +
                         a3 * (X ** 3) + a4 * (X ** 4))
        print 'chlor_a:', chlor_a.min(), chlor_a.mean(), chlor_a.max()

        img = geo_image.GeoImage(chlor_a,
                                 self.area,
                                 self.time_slot,
                                 fill_value=0,
                                 mode="L")

        if stretch:
            img.enhance(stretch=stretch)
        return img

    chlorophyll.prerequisites = set(["M02", "M03", "M04"])

    def hr_cloudtop(self):
        """Make a Night Fog RGB image composite.
        """
        self.check_channels('I04', 'I05')

        ch1 = -self['I04'].data
        ch2 = self['I05'].data
        ch3 = self['I05'].data

        img = geo_image.GeoImage((ch1, ch2, ch3),
                                 self.area,
                                 self.time_slot,
                                 fill_value=(0, 0, 0),
                                 mode="RGB")

        img.enhance(stretch=(0.005, 0.005))

        return img

    hr_cloudtop.prerequisites = set(['I04', 'I05'])

    def snow_age(self):
        """Make a Snow age RGB image composite.
        """
        self.check_channels('M07', 'M08', 'M09', 'M10', 'M11')

        coeff = 255. / 160.

        lonlats = self['M11'].area.get_lonlats()

        m07 = self['M07'].sunzen_corr(
            self.time_slot, lonlats, limit=88., sunmask=95).data * coeff
        m08 = self['M08'].sunzen_corr(
            self.time_slot, lonlats, limit=88., sunmask=95).data * coeff
        m09 = self['M09'].sunzen_corr(
            self.time_slot, lonlats, limit=88., sunmask=95).data * coeff
        m10 = self['M10'].sunzen_corr(
            self.time_slot, lonlats, limit=88., sunmask=95).data * coeff
        m11 = self['M11'].sunzen_corr(
            self.time_slot, lonlats, limit=88., sunmask=95).data * coeff

        refcu = m11 - m10
        refcu[refcu < 0] = 0

        ch1 = m07 - refcu / 2. - m09 / 4.
        ch2 = m08 + refcu / 4. + m09 / 4.
        ch3 = m11 + m09

        # Bernard Bellec snow Look-Up Tables V 1.0 (c) Meteo-France
        # These Look-up Tables allow you to create the RGB snow product
        # for SUOMI-NPP VIIRS Imager according to the algorithm
        # presented at the second CSPP/IMAPP users' meeting at Eumetsat
        # in Darmstadt on 14-16 April 2015
        # The algorithm and the product are described in this
        # presentation :
        # http://www.ssec.wisc.edu/meetings/cspp/2015/Agenda%20PDF/Wednesday/Roquet_snow_product_cspp2015.pdf
        #
        # For further information you may contact
        # Bernard Bellec at Bernard.Bellec@meteo.fr
        # or
        # Pascale Roquet at Pascale.Roquet@meteo.fr

        luts = np.array([[0, 0, 0], [1, 2, 2], [3, 8, 5], [4, 12, 8], [6, 15, 10], [8, 18, 13], [9, 21, 16],
                         [11, 24, 19], [13, 26, 21], [14, 28, 24], [
                             16, 30, 27], [18, 32, 30], [19, 34, 32],
                         [21, 36, 35], [22, 38, 38], [24, 40, 40], [
                             26, 42, 43], [27, 43, 46], [29, 45, 49],
                         [31, 47, 51], [32, 49, 54], [34, 50, 57], [
                             36, 52, 60], [37, 54, 62], [39, 55, 65],
                         [40, 57, 68], [42, 59, 70], [44, 60, 73], [
                             45, 62, 76], [47, 64, 79], [49, 66, 81],
                         [50, 67, 84], [52, 69, 87], [53, 71, 90], [
                             55, 73, 92], [56, 75, 95], [58, 77, 98],
                         [59, 79, 100], [61, 81, 103], [62, 83, 106], [
                             64, 85, 109], [65, 86, 111], [67, 88, 114],
                         [68, 90, 117], [70, 92, 119], [71, 94, 121], [
                             73, 96, 124], [74, 98, 126], [76, 100, 129],
                         [77, 102, 131], [79, 104, 134], [80, 106, 136], [
                             82, 107, 139], [83, 109, 141], [85, 111, 144],
                         [86, 113, 146], [88, 115, 149], [89, 117, 151], [
                             91, 118, 154], [92, 120, 156], [94, 122, 159],
                         [95, 124, 161], [97, 126, 162], [98, 128, 164], [
                             100, 129, 166], [101, 131, 168],
                         [103, 133, 170], [104, 135, 172], [106, 137, 173], [
                             107, 138, 175], [109, 140, 177],
                         [110, 142, 179], [112, 144, 181], [113, 145, 183], [
                             114, 147, 184], [116, 149, 186],
                         [117, 151, 188], [118, 152, 190], [120, 154, 192], [
                             121, 156, 193], [123, 158, 194],
                         [124, 159, 196], [125, 161, 197], [127, 163, 199], [
                             128, 165, 200], [130, 166, 202],
                         [131, 168, 203], [132, 170, 205], [134, 172, 206], [
                             135, 173, 206], [136, 175, 207],
                         [138, 177, 208], [139, 178, 209], [141, 180, 210], [
                             142, 182, 211], [143, 184, 212],
                         [145, 185, 213], [146, 187, 214], [148, 189, 215], [
                             149, 191, 216], [150, 192, 217],
                         [152, 194, 218], [153, 196, 219], [154, 198, 220], [
                             156, 200, 220], [157, 201, 221],
                         [159, 203, 221], [160, 205, 222], [161, 207, 223], [
                             162, 209, 223], [163, 210, 224],
                         [164, 212, 225], [166, 213, 225], [167, 214, 226], [
                             168, 216, 227], [169, 217, 227],
                         [171, 218, 228], [173, 220, 228], [174, 221, 228], [
                             175, 222, 229], [176, 224, 229],
                         [177, 225, 229], [178, 226, 230], [179, 227, 230], [
                             181, 228, 230], [182, 229, 231],
                         [183, 230, 231], [184, 231, 232], [185, 232, 232], [
                             186, 233, 232], [187, 234, 233],
                         [188, 235, 233], [190, 236, 233], [191, 237, 234], [
                             192, 237, 234], [193, 238, 234],
                         [194, 239, 235], [195, 240, 235], [196, 240, 236], [
                             196, 241, 236], [197, 242, 236],
                         [198, 243, 237], [199, 243, 237], [200, 244, 237], [
                             201, 245, 238], [202, 245, 238],
                         [203, 245, 238], [204, 246, 239], [205, 246, 239], [
                             206, 246, 239], [207, 247, 239],
                         [208, 247, 239], [209, 247, 239], [209, 248, 240], [
                             210, 248, 240], [210, 248, 240],
                         [211, 248, 240], [212, 248, 240], [212, 248, 241], [
                             213, 248, 241], [214, 248, 241],
                         [215, 248, 241], [216, 248, 241], [217, 248, 242], [
                             217, 248, 242], [218, 248, 242],
                         [219, 248, 242], [219, 248, 242], [220, 248, 243], [
                             221, 248, 243], [221, 249, 243],
                         [222, 249, 243], [223, 249, 243], [223, 249, 244], [
                             223, 249, 244], [224, 249, 244],
                         [224, 249, 244], [225, 249, 245], [225, 249, 245], [
                             226, 249, 245], [226, 249, 245],
                         [227, 249, 245], [227, 249, 246], [228, 249, 246], [
                             228, 250, 246], [229, 250, 246],
                         [229, 250, 246], [230, 250, 247], [230, 250, 247], [
                             231, 250, 247], [231, 250, 247],
                         [232, 250, 247], [233, 250, 248], [233, 250, 248], [
                             233, 250, 248], [234, 250, 248],
                         [234, 250, 248], [234, 250, 249], [235, 251, 249], [
                             235, 251, 249], [235, 251, 249],
                         [236, 251, 249], [236, 251, 250], [237, 251, 250], [
                             237, 251, 250], [237, 251, 250],
                         [238, 251, 250], [238, 251, 250], [238, 251, 250], [
                             239, 251, 250], [239, 251, 250],
                         [240, 251, 250], [240, 251, 250], [240, 252, 250], [
                             241, 252, 250], [241, 252, 251],
                         [241, 252, 251], [242, 252, 251], [242, 252, 251], [
                             242, 252, 251], [243, 252, 251],
                         [243, 252, 251], [244, 252, 251], [244, 252, 251], [
                             244, 252, 251], [245, 252, 252],
                         [245, 252, 252], [245, 253, 252], [246, 253, 252], [
                             246, 253, 252], [247, 253, 252],
                         [248, 253, 252], [248, 253, 252], [248, 253, 252], [
                             249, 253, 252], [249, 253, 253],
                         [249, 253, 253], [250, 253, 253], [250, 253, 253], [
                             250, 253, 253], [250, 253, 253],
                         [251, 254, 253], [251, 254, 253], [251, 254, 253], [
                             252, 254, 253], [252, 254, 254],
                         [252, 254, 254], [253, 254, 254], [253, 254, 254], [
                             253, 254, 254], [253, 254, 254],
                         [254, 254, 254], [254, 254, 254], [254, 254, 254], [254, 254, 254], [255, 255, 255]]) / 255.0
        np.ma.clip(ch1, 0, 255, ch1)
        np.ma.clip(ch2, 0, 255, ch2)
        np.ma.clip(ch3, 0, 255, ch3)
        ch1 = np.ma.array(
            luts[:, 0][ch1.astype(np.uint8)], copy=False, mask=ch1.mask)
        ch2 = np.ma.array(
            luts[:, 1][ch2.astype(np.uint8)], copy=False, mask=ch2.mask)
        ch3 = np.ma.array(
            luts[:, 2][ch3.astype(np.uint8)], copy=False, mask=ch3.mask)

        img = geo_image.GeoImage(
            (ch1, ch2, ch3), self.area, self.time_slot, mode="RGB")

        return img

    snow_age.prerequisites = set(['M07', 'M08', 'M09', 'M10', 'M11'])
