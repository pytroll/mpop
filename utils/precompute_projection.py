#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2017
#
# Author(s):
#
#   Panu Lahtinen <panu.lahtinen@fmi.fi>
#
# This file is part of mpop.
#
# mpop is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mpop is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mpop.  If not, see <http://www.gnu.org/licenses/>.

"""Module for calculating and projection mapping look-up tables, which
can be used with MPOP by setting 'precompute=True'"""

import numpy as np
import sys

from mpop.projector import (calc_nearest_params,
                            calc_bilinear_params,
                            calc_quick_params,
                            get_precompute_cache_fname,
                            get_area_and_id)


def save(fname, **data):
    """Save preprojection data to npz file."""
    np.savez(fname, **data)


def calc_preproj_params(out_dir, mode, in_area_name, out_area_name,
                        radius=None, nprocs=1):
    """Calculate preprojection parameters and save to *out_dir*"""
    in_area, in_id = get_area_and_id(in_area_name)
    out_area, out_id = get_area_and_id(out_area_name)

    fname = get_precompute_cache_fname(in_id, out_id,
                                       in_area_name, out_area_name,
                                       mode, out_dir)

    if mode == "nearest":
        data = calc_nearest_params(in_area, out_area,
                                   radius, nprocs=nprocs)
    elif mode == "bilinear":
        data = calc_bilinear_params(in_area, out_area, radius, nprocs=nprocs)
    elif mode == "quick":
        data = calc_quick_params(in_area, out_area)

    save(fname, **data)


def print_usage():
    """Print usage"""
    print("USAGE:")
    print("python precompute_projection.py <in_area_name> "
          "<out_area_name>")
    print("python precompute_projection.py <in_area_name> "
          "<out_area_name> <mode>")
    print("python precompute_projection.py <in_area_name> "
          "<out_area_name> <mode> <search_radius>")
    print("python precompute_projection.py <in_area_name> "
          "<out_area_name> <mode> <search_radius> <nprocs>")
    print("python precompute_projection.py <in_area_name> "
          "<out_area_name> <mode> <search_radius> <nprocs> <out_dir>")


def main():
    try:
        in_area_name = sys.argv[1]
        out_area_name = sys.argv[2]
    except IndexError:
        print_usage()
        return

    try:
        mode = sys.argv[3]
    except IndexError:
        mode = "nearest"
    try:
        radius = int(sys.argv[4])
    except IndexError:
        radius = 50000
    try:
        nprocs = int(sys.argv[5])
    except IndexError:
        nprocs = 1

    try:
        out_dir = sys.argv[6]
    except IndexError:
        out_dir = '.'

    calc_preproj_params(out_dir, mode, in_area_name, out_area_name,
                        radius=radius, nprocs=nprocs)


if __name__ == "__main__":
    main()
