#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2017, 2018.

# SMHI,
# Folkborgsvägen 1,
# Norrköping,
# Sweden

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>
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

"""This module defines satellite instrument channels as a generic class, to be
inherited when needed.
"""
import copy

import numpy as np
import logging

LOG = logging.getLogger(__name__)

try:
    from pyorbital.astronomy import sun_zenith_angle as sza
except ImportError:
    sza = None

from mpop.tools import viewzen_corr as vz_corr


class GeolocationIncompleteError(Exception):

    """Exception to try catch cases where the original data have not been read or
    expanded properly so that each pixel has a geo-location"""

    pass


class NotLoadedError(Exception):

    """Exception to be raised when attempting to use a non-loaded channel.
    """
    pass


class GenericChannel(object):

    """This is an abstract channel class. It can be a super class for
    calibrated channels data or more elaborate channels such as cloudtype or
    CTTH.
    """

    def __init__(self, name=None):
        object.__init__(self)

        # Channel name
        if name is not None and not isinstance(name, str):
            raise TypeError("Channel name must be a string, or None")
        self.name = name

        # Channel resolution, in meters.
        self.resolution = None

        # ID of the area on which the channel is defined.
        self.area_id = None

        # Area on which the channel is defined.
        self.area_def = None
        self.info = {}

    def __cmp__(self, ch2):
        if(isinstance(ch2, str)):
            return cmp(self.name, ch2)
        elif(ch2.name is not None and
             self.name is not None and
             ch2.name[0] == "_" and
             self.name[0] != "_"):
            return -1
        elif(ch2.name is not None and
             self.name is not None and
             ch2.name[0] != "_" and
             self.name[0] == "_"):
            return 1
        else:
            return cmp(self.name, ch2.name)

    def _get_area(self):
        """Getter for area.
        """
        return self.area_def or self.area_id

    def _set_area(self, area):
        """Setter for area.
        """
        if (area is None):
            self.area_def = None
            self.area_id = None
        elif(isinstance(area, str)):
            self.area_id = area
        else:
            try:
                dummy = area.area_extent
                dummy = area.x_size
                dummy = area.y_size
                dummy = area.proj_id
                dummy = area.proj_dict
                self.area_def = area
            except AttributeError:
                try:
                    dummy = area.lons
                    dummy = area.lats
                    self.area_def = area
                    self.area_id = None
                except AttributeError:
                    raise TypeError("Malformed area argument. "
                                    "Should be a string or an area object.")

    area = property(_get_area, _set_area)


class Channel(GenericChannel):

    """This is the satellite channel class. It defines satellite channels as a
    container for calibrated channel data.

    The *resolution* sets the resolution of the channel, in meters. The
    *wavelength_range* is a triplet, containing the lowest-, center-, and
    highest-wavelength values of the channel. *name* is simply the given name
    of the channel, and *data* is the data it should hold.
    """

    def __init__(self,
                 name=None,
                 resolution=0,
                 wavelength_range=[-np.inf, -np.inf, -np.inf],
                 data=None,
                 calibration_unit=None):

        GenericChannel.__init__(self, name)

        self._data = None
        self.wavelength_range = None

        if(name is None and
           wavelength_range == [-np.inf, -np.inf, -np.inf]):
            raise ValueError("Cannot define a channel with neither name "
                             "nor wavelength range.")

        if not isinstance(resolution, (int, float)):
            raise TypeError("Resolution must be an integer number of meters.")

        self.resolution = resolution

        if(not isinstance(wavelength_range, (tuple, list, set)) or
           len(wavelength_range) != 3 or
           not isinstance(wavelength_range[0], float) or
           not isinstance(wavelength_range[1], float) or
           not isinstance(wavelength_range[2], float)):
            raise TypeError("Wavelength_range should be a triplet of floats.")
        elif(not (wavelength_range[0] <= wavelength_range[1]) or
             not (wavelength_range[1] <= wavelength_range[2])):
            raise ValueError("Wavelength_range should be a sorted triplet.")

        self.wavelength_range = list(wavelength_range)
        self.unit = calibration_unit
        self.data = data

    def get_reflectance(self, tb11, sun_zenith=None, tb13_4=None):
        """Get the reflectance part of an NIR channel"""

        try:
            from pyspectral.near_infrared_reflectance import Calculator
        except ImportError:
            LOG.info("Couldn't load pyspectral")

        # Check the wavelength, and if outside 3-4 microns this functionality
        # doesn't give any meaning and should not be supported
        if (self.wavelength_range[1] < 3.0 or self.wavelength_range[1] > 4.0):
            LOG.warning("Deriving the near infrared reflectance" +
                        " of a band that is outside the 3-4 micron range" +
                        " is not supported!\n\tWill do nothing...")
            return

        # Check if the sun-zenith angle was provided:
        if sun_zenith is None:
            lonlats = self.area.get_lonlats()
            sun_zenith = sza(self.info['time'], lonlats[0], lonlats[1])

        try:
            refl39 = Calculator(self.info['satname'] + self.info['satnumber'],
                                self.info['instrument_name'], self.name)
        except NameError:
            LOG.warning("pyspectral missing!")
            return

        return refl39.reflectance_from_tbs(sun_zenith, self.data, tb11, tb_ir_co2=tb13_4)

    def __cmp__(self, ch2, key=0):
        if(isinstance(ch2, str)):
            return cmp(self.name, ch2)
        elif(ch2.name is not None and
             self.name is not None and
             ch2.name[0] == "_" and
             self.name[0] != "_"):
            return -1
        elif(ch2.name is not None and
             self.name is not None and
             ch2.name[0] != "_" and
             self.name[0] == "_"):
            return 1
        else:
            res = cmp(abs(self.wavelength_range[1] - key),
                      abs(ch2.wavelength_range[1] - key))
            if res == 0:
                return cmp(self.name, ch2.name)
            else:
                return res

    def __str__(self):
        if self.shape is not None:
            return ("'%s: (%.3f,%.3f,%.3f)μm, shape %s, resolution %sm'" %
                    (self.name,
                     self.wavelength_range[0],
                     self.wavelength_range[1],
                     self.wavelength_range[2],
                     self.shape,
                     self.resolution))
        else:
            return ("'%s: (%.3f,%.3f,%.3f)μm, resolution %sm, not loaded'" %
                    (self.name,
                     self.wavelength_range[0],
                     self.wavelength_range[1],
                     self.wavelength_range[2],
                     self.resolution))

    def is_loaded(self):
        """Tells if the channel contains loaded data.
        """
        return self._data is not None

    def check_range(self, min_range=1.0):
        """Check that the data of the channels has a definition domain broader
        than *min_range* and return the data, otherwise return zeros.
        """
        if not self.is_loaded():
            raise ValueError("Cannot check range of an non-loaded channel")

        if not isinstance(min_range, (float, int)):
            raise TypeError("Min_range must be a single number.")

        if isinstance(self._data, np.ma.core.MaskedArray):
            if self._data.mask.all():
                return self._data

        if((self._data.max() - self._data.min()) < min_range):
            return np.ma.zeros(self.shape)
        else:
            return self._data

    def show(self):
        """Display the channel as an image.
        """
        if not self.is_loaded():
            raise ValueError("Channel not loaded, cannot display.")

        from PIL import Image as pil

        data = ((self._data - self._data.min()) * 255.0 /
                (self._data.max() - self._data.min()))
        if isinstance(data, np.ma.core.MaskedArray):
            img = pil.fromarray(np.array(data.filled(0), np.uint8))
        else:
            img = pil.fromarray(np.array(data, np.uint8))
        img.show()

    def as_image(self, stretched=True):
        """Return the channel as a :class:`mpop.imageo.geo_image.GeoImage`
        object. The *stretched* argument set to False allows the data to remain
        untouched (as opposed to crude stretched by default to obtain the same
        output as :meth:`show`).
        """
        from mpop.imageo.geo_image import GeoImage

        img = GeoImage(self._data, self.area, None)
        if stretched:
            img.stretch("crude")
        return img

    def project(self, coverage_instance):
        """Make a projected copy of the current channel using the given
        *coverage_instance*.

        See also the :mod:`mpop.projector` module.
        """
        res = Channel(name=self.name,
                      resolution=self.resolution,
                      wavelength_range=self.wavelength_range,
                      data=None,
                      calibration_unit=self.unit)
        res.area = coverage_instance.out_area
        res.info = self.info
        if hasattr(self, 'palette'):      # UH, new
            res.palette = self.palette    # UH, new
        if self.is_loaded():
            LOG.info("Projecting channel %s (%fμm)..."
                     % (self.name, self.wavelength_range[1]))
            import pyresample
            if (hasattr(coverage_instance, 'in_area') and
                isinstance(coverage_instance.in_area, pyresample.geometry.SwathDefinition) and
                    hasattr(coverage_instance.in_area.lats, 'shape') and
                    coverage_instance.in_area.lats.shape != self._data.shape):
                raise GeolocationIncompleteError("Lons and lats doesn't match data! " +
                                                 "Data can't be re-projected unless " +
                                                 "each pixel of the swath has a " +
                                                 "geo-location atached to it.")
            data = coverage_instance.project_array(self._data)
            res.data = data
            return res
        else:
            raise NotLoadedError("Can't project, channel %s (%fμm) not loaded."
                                 % (self.name, self.wavelength_range[1]))

    def get_data(self):
        """Getter for channel data.
        """
        return self._data

    def set_data(self, data):
        """Setter for channel data.
        """
        if data is None:
            del self._data
            self._data = None
        elif isinstance(data, (np.ndarray, np.ma.core.MaskedArray)):
            self._data = data
        else:
            raise TypeError("Data must be a numpy (masked) array.")

    data = property(get_data, set_data)

    @property
    def shape(self):
        """Shape of the channel.
        """
        if self.data is None:
            return None
        else:
            return self.data.shape

    def sunzen_corr(self, time_slot, lonlats=None, limit=80., mode='cos',
                    sunmask=False):
        '''Perform Sun zenith angle correction for the channel at
        *time_slot* (datetime.datetime() object) and return the
        corrected channel.  The parameter *limit* can be used to set
        the maximum zenith angle for which the correction is
        calculated.  For larger angles, the correction is the same as
        at the *limit* (default: 80.0 degrees).  Coordinate values can
        be given as a 2-tuple or a two-element list *lonlats* of numpy
        arrays; if None, the coordinates will be read from the channel
        data.  Parameter *mode* is a placeholder for other possible
        illumination corrections. The name of the new channel will be
        *original_chan.name+'_SZC'*, eg. "VIS006_SZC".  This name is
        also stored to the info dictionary of the originating channel.
        '''

        if self.info.get('sun_zen_correction_applied'):
            LOG.debug("Sun zenith correction already applied, skipping")
            return self

        import mpop.tools

        try:
            from pyorbital import astronomy
        except ImportError:
            LOG.warning("Could not load pyorbital.astronomy")
            return None

        if lonlats is None or len(lonlats) != 2:
            # Read coordinates
            LOG.debug("No valid coordinates given, reading from the "
                      "channel data")
            lons, lats = self.area.get_lonlats()
        else:
            lons, lats = lonlats

        # Calculate Sun zenith angles and the cosine
        cos_zen = astronomy.cos_zen(time_slot, lons, lats)

        # Copy the channel
        new_ch = copy.deepcopy(self)

        # Set the name
        new_ch.name += '_SZC'

        if mode == 'cos':
            new_ch.data = mpop.tools.sunzen_corr_cos(new_ch.data,
                                                     cos_zen, limit=limit)
        else:
            # Placeholder for other correction methods
            pass

        # Add information about the corrected version to original
        # channel
        self.info["sun_zen_corrected"] = self.name + '_SZC'

        if sunmask:
            if isinstance(sunmask, (float, int)):
                sunmask = sunmask
            else:
                sunmask = 90.
            cos_limit = np.cos(np.radians(sunmask))
            LOG.debug("Masking out data where sun-zenith " +
                      "is greater than %f deg", sunmask)
            LOG.debug("cos_limit = %f", cos_limit)
            # Mask out data where the sun elevation is below a threshold:
            new_ch.data = np.ma.masked_where(
                cos_zen < cos_limit, new_ch.data, copy=False)

        new_ch.info["sun_zen_correction_applied"] = True

        return new_ch

    def get_viewing_geometry(self, orbital, time_slot, altitude=None):
        '''Calculates the azimuth and elevation angle as seen by the observer 
           at the position of the current area pixel. 
           inputs:
             orbital   an orbital object define by the tle file 
                       (see pyorbital.orbital import Orbital or mpop/scene.py get_oribtal)
             time_slot time object specifying the observation time
             altitude  optinal: altitude of the observer above the earth ellipsoid
           outputs:
             azi       azimuth viewing angle in degree (south is 0, counting clockwise)
             ele       elevation viewing angle in degree (zenith is 90, horizon is 0)
        '''

        try:
            from pyorbital.orbital import Orbital
        except ImportError:
            LOG.warning("Could not load pyorbital.orbial.Orbital")
            return None

        try:
            from pyorbital import tlefile
        except ImportError:
            LOG.warning("Could not load pyorbital.tlefile")
            return None

        (lons, lats) = self.area.get_lonlats()
        # Calculate observer azimuth and elevation
        if altitude == None:
            altitude = np.zeros(lons.shape)
        azi, ele = orbital.get_observer_look(time_slot, lons, lats, altitude)

        return (azi, ele)

    def vinc_vect(phi, lembda, alpha, s, f=None, a=None, degree=True):
        """ Vincenty's Direct formular

        Returns the lat and long of projected point and reverse azimuth
        given a reference point and a distance and azimuth to project.
        lats, longs and azimuths are passed in radians.

        Keyword arguments:
            phi    Latitude in degree/radians
            lembda Longitude in degree/radians
            alpha    Geodetic azimuth in degree/radians
            s    Ellipsoidal distance in meters
            f    WGS84 parameter
            a    WGS84 parameter
            degree Boolean if in/out values are in degree or radians.
                   Default is in degree

        Returns:
            (phiout,  lembdaout,  alphaout ) as a tuple

        """
        if degree:
            phi = np.deg2rad(phi)
            lembda = np.deg2rad(lembda)
            alpha = np.deg2rad(alpha)

        if f is None:
            f = 1 / 298.257223563
        if a is None:
            a = 6378137

        two_pi = 2.0 * np.pi

        if isinstance(alpha, np.ndarray):
            alpha[alpha < 0.0] += two_pi
            alpha[alpha > two_pi] -= two_pi

        else:
            if alpha < 0.0:
                alpha = alpha + two_pi
            if (alpha > two_pi):
                alpha = alpha - two_pi
        """
        alphama = np.ma.masked_less_equal(alphama, two_pi)
        alpha = alphama - two_pi
        alpha.mask = np.ma.nomask
        logger.debug(alpha)
        """
        b = a * (1.0 - f)

        tan_u1 = (1 - f) * np.tan(phi)
        u_1 = np.arctan(tan_u1)
        sigma1 = np.arctan2(tan_u1, np.cos(alpha))

        sinalpha = np.cos(u_1) * np.sin(alpha)
        cosalpha_sq = 1.0 - sinalpha * sinalpha

        u_2 = cosalpha_sq * (a * a - b * b) / (b * b)
        aa_ = 1.0 + (u_2 / 16384) * (4096 + u_2 * (-768 + u_2 *
                                                   (320 - 175 * u_2)))
        bb_ = (u_2 / 1024) * (256 + u_2 * (-128 + u_2 * (74 - 47 * u_2)))

        # Starting with the approximation
        sigma = (s / (b * aa_))
        last_sigma = 2.0 * sigma + 2.0  # something impossible

        # Iterate the following three equations
        # until there is no significant change in sigma

        # two_sigma_m , delta_sigma

        def iter_sigma(sigma, last_sigma, sigma1, s, b, aa_, bb_):
            while (abs((last_sigma - sigma) / sigma) > 1.0e-9):
                two_sigma_m = 2 * sigma1 + sigma

                delta_sigma = (bb_ * np.sin(sigma) *
                               (np.cos(two_sigma_m) + (bb_ / 4) *
                                (np.cos(sigma) *
                                 (-1 + 2 * np.power(np.cos(two_sigma_m), 2) -
                                  (bb_ / 6) * np.cos(two_sigma_m) *
                                  (-3 + 4 * np.power(np.sin(sigma), 2)) *
                                  (-3 + 4 * np.power(np.cos(two_sigma_m), 2))))))
                last_sigma = sigma
                sigma = (s / (b * aa_)) + delta_sigma

            return(sigma, two_sigma_m)

        # Check for array inputs
        arraybool = [isinstance(ele, np.ndarray) for ele in (sigma, last_sigma,
                                                             sigma1)]
        logger.debug("Sigma Arrays?: " + str(arraybool))
        if all(arraybool):
            viter_sigma = np.vectorize(iter_sigma)
            sigma, two_sigma_m = viter_sigma(sigma, last_sigma, sigma1, s, b, aa_,
                                             bb_)

        else:
            sigma, two_sigma_m = iter_sigma(sigma, last_sigma, sigma1, s, b, aa_,
                                            bb_)

        phiout = np.arctan2((np.sin(u_1) * np.cos(sigma) +
                             np.cos(u_1) * np.sin(sigma) * np.cos(alpha)),
                            ((1 - f) * np.sqrt(np.power(sinalpha, 2) +
                                               pow(np.sin(u_1) *
                                                   np.sin(sigma) -
                                                   np.cos(u_1) *
                                                   np.cos(sigma) *
                                                   np.cos(alpha), 2))))

        deltalembda = np.arctan2((np.sin(sigma) * np.sin(alpha)),
                                 (np.cos(u_1) * np.cos(sigma) -
                                  np.sin(u_1) * np.sin(sigma) * np.cos(alpha)))

        cc_ = (f / 16) * cosalpha_sq * (4 + f * (4 - 3 * cosalpha_sq))

        omega = (deltalembda - (1 - cc_) * f * sinalpha *
                 (sigma + cc_ * np.sin(sigma) * (np.cos(two_sigma_m) + cc_ *
                                                 np.cos(sigma) *
                                                 (-1 + 2 *
                                                  np.power(np.cos(two_sigma_m),
                                                           2)))))

        lembdaout = lembda + omega

        alphaout = np.arctan2(sinalpha, (-np.sin(u_1) * np.sin(sigma) +
                                         np.cos(u_1) * np.cos(sigma) *
                                         np.cos(alpha)))

        alphaout = alphaout + two_pi / 2.0

        if isinstance(alphaout, np.ndarray):
            alphaout[alphaout < 0.0] += two_pi
            alphaout[alphaout > two_pi] -= two_pi

        else:
            if alphaout < 0.0:
                alphaout = alphaout + two_pi
            if (alphaout > two_pi):
                alphaout = alphaout - two_pi

        if degree:
            phiout = np.rad2deg(phiout)
            lembdaout = np.rad2deg(lembdaout)
            alphaout = np.rad2deg(alphaout)

        return(phiout, lembdaout, alphaout)

    def parallax_corr(self, cth=None, time_slot=None, orbital=None, azi=None, ele=None, fill="False"):
        '''Perform the parallax correction for channel at
        *time_slot* (datetime.datetime() object), assuming the cloud top height cth
        and the viewing geometry given by the satellite orbital "orbital" and return the
        corrected channel. 
        Authors: Ulrich Hamann (MeteoSwiss), Thomas Leppelt (DWD)
        Example calls:
            * calling this function (using orbital and time_slot)
                 orbital = data.get_oribtal()
                 data["VIS006"].parallax_corr(cth=data["CTTH"].height, time_slot=data.time_slot, orbital=orbital)
            * calling this function (using viewing geometry)
                 orbital = data.get_oribtal()
                 (azi, ele) = get_viewing_geometry(self, orbital, time_slot)
                 data["VIS006"].parallax_corr(cth=data["CTTH"].height, azi=azi, ele=ele)
        Optional input:
          cth        The parameter cth is the cloud top height 
                     (or  the altitude of the object that should be shifted).
                     cth must have the same size and projection as the channel

          orbital    an orbital object define by the tle file 
                     (see pyorbital.orbital import Orbital or mpop/scene.py get_oribtal)
          azi        azimuth viewing angle in degree (south is 0, counting clockwise)
                     e.g. as given by self.get_viewing_geometry
          ele        elevation viewing angle in degree (zenith is 90, horizon is 0)
                     e.g. as given by self.get_viewing_geometry
          fill       specifies the interpolation method to fill the gaps
                     (basically areas behind the cloud that can't be observed by the satellite instrument)
                     "False" (default): no interpolation, gaps are np.nan values and mask is set accordingly
                     "nearest": fill gaps with nearest neighbour
                     "bilinear": use scipy.interpolate.griddata with linear interpolation 
                                 to fill the gaps

        output: 
          parallax corrected channel
                     the content of the channel will be parallax corrected.
                     The name of the new channel will be
                     *original_chan.name+'_PC'*, eg. "IR_108_PC". This name is
                     also stored to the info dictionary of the originating channel.
        '''

        # get time_slot from info, if present
        if time_slot == None:
            if "time" in self.info.keys():
                time_slot = self.info["time"]

        if azi == None or ele == None:
            if time_slot == None or orbital == None:
                print "*** Error in parallax_corr (mpop/channel.py)"
                print "    parallax_corr needs either time_slot and orbital"
                print "    data[\"IR_108\"].parallax_corr(data[\"CTTH\"].height, time_slot=data.time_slot, orbital=orbital)"
                print "    or the azimuth and elevation angle"
                print "    data[\"IR_108\"].parallax_corr(data[\"CTTH\"].height, azi=azi, ele=ele)"
                quit()
            else:
                print (
                    "... calculate viewing geometry (orbit and time are given)")
                (azi, ele) = self.get_viewing_geometry(orbital, time_slot)
        else:
            print ("... azimuth and elevation angle given")

        # mask the cloud top height
        cth_ = np.ma.masked_where(cth < 0, cth, copy=False)

        # Elevation displacement
        dz = cth_ / np.tan(np.deg2rad(ele))

        # Create the new channel (by copying) and initialize the data with None
        # values
        new_ch = copy.deepcopy(self)
        new_ch.data[:, :] = np.nan

        # Set the name
        new_ch.name += '_PC'

        # Add information about the corrected version to original channel
        self.info["parallax_corrected"] = self.name + '_PC'

        # get projection coordinates in meter
        (proj_x, proj_y) = self.area.get_proj_coords()

        print "... calculate parallax shift"
        # shifting pixels according to parallax corretion
        # shift West-East   in m  # ??? sign correct ???
        proj_x_pc = proj_x - np.sin(np.deg2rad(azi)) * dz
        # shift North-South in m
        proj_y_pc = proj_y + np.cos(np.deg2rad(azi)) * dz

        # get indices for the pixels for the original position
        (y, x) = self.area.get_xy_from_proj_coords(proj_x, proj_y)
        # comment: might be done more efficient with meshgrid
        # >>> x = np.arange(-5.01, 5.01, 0.25)
        # >>> y = np.arange(-5.01, 5.01, 0.25)
        # >>> xx, yy = np.meshgrid(x, y)
        # get indices for the pixels at the parallax corrected position
        (y_pc, x_pc) = self.area.get_xy_from_proj_coords(proj_x_pc, proj_y_pc)

        # copy cloud free satellite pixels (surface observations)
        ind = np.where(cth_.mask == True)
        new_ch.data[x[ind], y[ind]] = self.data[x[ind], y[ind]]

        print "... copy data to parallax corrected position"
        # copy cloudy pixel with new position modified with parallax shift
        ind = np.where(x_pc.mask == False)
        new_ch.data[x_pc[ind], y_pc[ind]] = self.data[x[ind], y[ind]]

        # Mask out data gaps (areas behind the clouds)
        new_ch.data = np.ma.masked_where(
            np.isnan(new_ch.data), new_ch.data, copy=False)

        if fill.lower() == "false":
            return new_ch
        elif fill == "nearest":
            print "*** fill missing values with nearest neighbour"
            from scipy.ndimage import distance_transform_edt
            invalid = np.isnan(new_ch.data)
            ind = distance_transform_edt(
                invalid, return_distances=False, return_indices=True)
            new_ch.data = new_ch.data[tuple(ind)]
        elif fill == "bilinear":
            # this function does not interpolate at the outer boundaries
            from scipy.interpolate import griddata
            ind = np.where(new_ch.data.mask == False)
            points = np.transpose(np.append([y[ind]], [x[ind]], axis=0))
            values = new_ch.data[ind]
            new_ch.data = griddata(points, values, (y, x), method='linear')

            # fill the remaining pixels with nearest neighbour
            from scipy.ndimage import distance_transform_edt
            invalid = np.isnan(new_ch.data)
            ind = distance_transform_edt(
                invalid, return_distances=False, return_indices=True)
            new_ch.data = new_ch.data[tuple(ind)]
        else:
            print "*** Error in parallax_corr (channel.py)"
            print "    unknown gap fill method ", fill
            quit()

        return new_ch

    def viewzen_corr(self, view_zen_angle_data):
        """Apply atmospheric correction on a copy of this channel data
        using the given satellite zenith angle data of the same shape.
        Returns a new channel containing the corrected data.
        The name of the new channel will be *original_chan.name+'_VZC'*,
        eg. "IR108_VZC".  This name is also stored to the info dictionary of
        the originating channel.
        """

        # copy channel data which will be corrected in place
        chn_data = self.data.copy()
        CHUNK_SZ = 500
        for start in xrange(0, chn_data.shape[1], CHUNK_SZ):
            # apply correction on channel data
            vz_corr(chn_data[:, start:start + CHUNK_SZ],
                    view_zen_angle_data[:, start:start + CHUNK_SZ])

        new_ch = Channel(name=self.name + "_VZC",
                         resolution=self.resolution,
                         wavelength_range=self.wavelength_range,
                         data=chn_data,
                         calibration_unit=self.unit)

        # Add information about the corrected version to original channel
        self.info["view_zen_corrected"] = self.name + '_VZC'

        return new_ch

    # Arithmetic operations on channels.

    def __pow__(self, other):
        return Channel(name="new", data=self.data ** other)

    def __rpow__(self, other):
        return Channel(name="new", data=self.data ** other)

    def __mul__(self, other):
        return Channel(name="new", data=self.data * other)

    def __rmul__(self, other):
        return Channel(name="new", data=self.data * other)

    def __add__(self, other):
        return Channel(name="new", data=self.data + other)

    def __radd__(self, other):
        return Channel(name="new", data=self.data + other)

    def __sub__(self, other):
        return Channel(name="new", data=self.data - other)

    def __rsub__(self, other):
        return Channel(name="new", data=self.data - other)

    def __div__(self, other):
        return Channel(name="new", data=self.data / other)

    def __rdiv__(self, other):
        return Channel(name="new", data=self.data / other)

    def __neg__(self):
        return Channel(name="new", data=-self.data)

    def __abs__(self):
        return Channel(name="new", data=abs(self.data))
