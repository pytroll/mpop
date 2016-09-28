# -*- coding: utf-8 -*-
# Copyright (c) 2014, 2015
#
# Author(s):
#
#   Panu Lahtinen <pnuu+git@iki.fi>
#
# This file is part of mpop.
#
# mpop is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# mpop is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# mpop.  If not, see <http://www.gnu.org/licenses/>.

'''Helper functions for eg. performing Sun zenith angle correction.
'''

import numpy as np


def sunzen_corr_cos(data, cos_zen, limit=80.):
    '''Perform Sun zenith angle correction to the given *data* using
    cosine of the zenith angle (*cos_zen*).  The correction is limited
    to *limit* degrees (default: 80.0 degrees).  For larger zenith
    angles, the correction is the same as at the *limit*.  Both *data*
    and *cos_zen* are given as 2-dimensional Numpy arrays or Numpy
    MaskedArrays, and they should have equal shapes.
    '''

    # Convert the zenith angle limit to cosine of zenith angle
    cos_limit = np.cos(np.radians(limit))

    # Cosine correction
    lim_y, lim_x = np.where(cos_zen > cos_limit)
    data[lim_y, lim_x] /= cos_zen[lim_y, lim_x]
    # Use constant value (the limit) for larger zenith
    # angles
    lim_y, lim_x = np.where(cos_zen <= cos_limit)
    data[lim_y, lim_x] /= cos_limit

    return data


def estimate_cth(IR_108, cth_atm="standard"):

    '''
    Estimation of the cloud top height using the 10.8 micron channel
    limitations: this is the most simple approach
                 a simple fit of the ir108 to the temperature profile
                 * no correction for water vapour or any other trace gas
                 * no viewing angle dependency
                 * no correction for semi-transparent clouds

    optional input:
      cth_atm    * "standard", "tropics", "midlatitude summer", "midlatitude winter", "subarctic summer", "subarctic winter"
                  Matching the 10.8 micron temperature with atmosphere profile
                  (s)  AFGL atmospheric constituent profile. U.S. standard atmosphere 1976. (AFGL-TR-86-0110) 
                  (t)  AFGL atmospheric constituent profile. tropical.                      (AFGL-TR-86-0110)
                  (mw) AFGL atmospheric constituent profile. midlatitude summer.            (AFGL-TR-86-0110) 
                  (ms) AFGL atmospheric constituent profile. midlatitude winter.            (AFGL-TR-86-0110)
                  (ss) AFGL atmospheric constituent profile. subarctic summer.              (AFGL-TR-86-0110) 
                  (sw) AFGL atmospheric constituent profile. subarctic winter.              (AFGL-TR-86-0110)
                  Ulrich Hamann (MeteoSwiss)
                * "tropopause"
                  Assuming a fixed tropopause height and a fixed temperature gradient
                  Richard Mueller (DWD)
    output: 
      parallax corrected channel
                 the content of the channel will be parallax corrected.
                 The name of the new channel will be
                 *original_chan.name+'_PC'*, eg. "IR_108_PC". This name is
                 also stored to the info dictionary of the originating channel.

    Versions: 05.07.2016 initial version
              Ulrich Hamann (MeteoSwiss), Richard Mueller (DWD)
    '''

    print "*** estimating CTH using the 10.8 micro meter brightness temperature "

    if cth_atm.lower() != "tropopause":

        # define atmospheric temperature profile    
        import os
        from numpy import loadtxt, zeros, where, logical_and
        import mpop 

        mpop_dir = os.path.dirname(mpop.__file__)
        afgl_file = mpop_dir+"/afgl.dat"
        print "... assume ", cth_atm, " atmosphere for temperature profile"

        if cth_atm.lower()=="standard" or cth_atm.lower()=="s":
            z, T = loadtxt(afgl_file, usecols=(0, 1), unpack=True, comments="#")
        elif cth_atm.lower()=="tropics" or cth_atm.lower()=="t":
            z, T = loadtxt(afgl_file, usecols=(0, 2), unpack=True, comments="#")
        elif cth_atm.lower()=="midlatitude summer" or cth_atm.lower()=="ms":
            z, T = loadtxt(afgl_file, usecols=(0, 3), unpack=True, comments="#")
        elif cth_atm.lower()=="midlatitude winter" or cth_atm.lower()=="ws":
            z, T = loadtxt(afgl_file, usecols=(0, 4), unpack=True, comments="#")
        elif cth_atm.lower()=="subarctic summer" or cth_atm.lower()=="ss":
            z, T = loadtxt(afgl_file, usecols=(0, 5), unpack=True, comments="#")
        elif cth_atm.lower()=="subarctic winter" or cth_atm.lower()=="ss":
            z, T = loadtxt(afgl_file, usecols=(0, 6), unpack=True, comments="#")
        else:
            print "*** Error in estimate_cth (mpop/tools.py)"
            print "unknown temperature profiel for CTH estimation: cth_atm = ", cth_atm
            quit()

        height = zeros(IR_108.shape)
        # warmer than lowest level -> clear sky 
        height[where(IR_108 > T[-1])] = -1.
        print "     z0(km)   z1(km)   T0(K)   T1(K)  number of pixels"
        print "------------------------------------------------------"
        for i in range(z.size)[::-1]:

            # search for temperatures between layer i-1 and i
            ind =  np.where( logical_and( T[i-1]< IR_108, IR_108 < T[i]) )
            # interpolate CTH according to ir108 temperature
            height[ind] = z[i] + (IR_108[ind]-T[i])/(T[i-1]-T[i]) * (z[i-1]-z[i])
            # verbose output
            print " {0:8.1f} {1:8.1f} {2:8.1f} {3:8.1f} {4:8d}".format(z[i], z[i-1], T[i], T[i-1], len(ind[0]))

            # if temperature increases above 8km -> tropopause detected
            if z[i]>=8. and T[i] <= T[i-1]:
                # no cloud above tropopose
                break
            # no cloud heights above 20km
            if z[i]>=20.:
                break

        # if height is still 0 -> cloud colder than tropopause -> cth == tropopause height
        height[np.where( height == 0 )] = z[i]
        
    else:

        Htropo=11.0 # km
        # this is an assumption it should be optimized 
        # by making it dependent on region and season. 
        # It might be good to include the ITC in the  
        # region of interest, that would make a fixed Htropo 
        # value more reliable. 
        Tmin = np.amin(IR_108) 
        # for Tmin it might be better to use the 5th or 10th percentile 
        # else overshoting tops induces further uncertainties  
        # in the calculation of the cloud height. 
        # However numpy provides weird results for 5th percentile. 
        # Hence, for the working version the minima is used 

        print "... assume tropopause height ", Htropo, ", tropopause temperature ", Tmin, "K (", Tmin-273.16, "deg C)"
        print "    and constant temperature gradient 6.5 K/km"

        height = -(IR_108 - Tmin)/6.5 + Htropo 
        # calculation of the height, the temperature gradient 
        # 6.5 K/km is an assumption  
        # derived from USS and MPI standard profiles. It 
        # has to be improved as well 

    # convert to masked array
    # convert form km to meter
    height = np.ma.masked_where(height <= 0, height, copy=False) * 1000.

    if False:
        from trollimage.image import Image as trollimage
        from trollimage.colormap import rainbow
        from copy import deepcopy 
        # cloud top height
        prop = height
        min_data = prop.min()
        max_data = prop.max()
        print " estimated CTH(meter) (min/max): ", min_data, max_data
        min_data =     0
        max_data = 12000    
        colormap = deepcopy(rainbow)
        colormap.set_range(min_data, max_data)
        img = trollimage(prop, mode="L") #, fill_value=[0,0,0]
        img.colorize(colormap)
        img.show()

    # return cloud top height in meter
    return height


def viewzen_corr(data, view_zen):
    """Apply atmospheric correction on the given *data* using the
    specified satellite zenith angles (*view_zen*). Both input data
    are given as 2-dimensional Numpy (masked) arrays, and they should
    have equal shapes.
    The *data* array will be changed in place and has to be copied before.
    """
    def ratio(value, v_null, v_ref):
        return (value - v_null) / (v_ref - v_null)

    def tau0(t):
        T_0 = 210.0
        T_REF = 320.0
        TAU_REF = 9.85
        return (1 + TAU_REF)**ratio(t, T_0, T_REF) - 1

    def tau(t):
        T_0 = 170.0
        T_REF = 295.0
        TAU_REF = 1.0
        M = 4
        return TAU_REF * ratio(t, T_0, T_REF)**M

    def delta(z):
        Z_0 = 0.0
        Z_REF = 70.0
        DELTA_REF = 6.2
        return (1 + DELTA_REF)**ratio(z, Z_0, Z_REF) - 1

    y0, x0 = np.ma.where(view_zen == 0)
    data[y0, x0] += tau0(data[y0, x0])

    y, x = np.ma.where((view_zen > 0) & (view_zen < 90) & (~data.mask))
    data[y, x] += tau(data[y, x]) * delta(view_zen[y, x])
