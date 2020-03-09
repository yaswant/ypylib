"""
(c) Crown copyright 2018, the Met Office.

Demonstrate how to convert (to regular grid) and extract region from a cube in
rotated pole grid
"""
# -*- coding: utf-8 -*-
# author: Yaswant Pradhan
# date: 2018-06-11

from __future__ import print_function
import os
import iris
import numpy as np
import matplotlib.pyplot as plt
import iris.quickplot as qplt


def coord_rot2reg(rotated_lons, rotated_lats, pole_lon=177.5, pole_lat=37.5):
    """Convert coordinate values from rotated pole grid to regular grid.

    Parameters
    ----------
    rotated_lons : array_like
        Longitude values in rotated pole grid.
    rotated_lats : array_like
        Latitude values in rotated pole grid.
    pole_lon : real
        The true longitude of the rotated pole in degrees.
    pole_lat : real
        The true latitude of the rotated pole in degrees.
    """
    lons, lats = iris.analysis.cartography.unrotate_pole(
        np.asarray(rotated_lons, 'float'),
        np.asarray(rotated_lats, 'float'),
        pole_lon, pole_lat)
    return lons, lats


def coord_reg2rot(lons, lats, pole_lon=177.5, pole_lat=37.5):
    """Convert coordinate values from regular grid to rotated pole grid.

    Parameters
    ----------
    lons : array_like
        Longitude values in regular grid.
    lats : array_like
        Latitude values in regular grid.
    pole_lon : real
        The true longitude of the rotated pole in degrees.
    pole_lat : real
        The true latitude of the rotated pole in degrees.
    """
    rlon, rlat = iris.analysis.cartography.rotate_pole(
        np.asarray(lons, 'float'),
        np.asarray(lats, 'float'),
        pole_lon, pole_lat)
    return rlon, rlat


def reg_2dcube_template(latlim=None, lonlim=None, nlat=None, nlon=None):
    """Create a 2-dimensional cube template with regular grid coordinates.

    Parameters
    ----------
    latlim : sequence
        bottom-left and top-right corner latitude values
    lonlim : sequence
        bottom-left and top-right corner longitude values
    nlat : int
        number of linearly spaced latitudes
    nlon : int
        number of linearly spaced longitudes
    """
    latitude = iris.coords.DimCoord(
        np.linspace(latlim[0], latlim[1], nlat),
        standard_name='latitude', units='degrees',
        coord_system=iris.coord_systems.GeogCS(6371229))
    longitude = iris.coords.DimCoord(
        np.linspace(lonlim[0], lonlim[1], nlon),
        standard_name='longitude', units='degrees',
        coord_system=iris.coord_systems.GeogCS(6371229))
    cube = iris.cube.Cube(
        np.zeros((nlat, nlon), np.float32),
        dim_coords_and_dims=[(latitude, 0), (longitude, 1)])
    return cube


def reg_cube_template(rot_cube, latlim=None, lonlim=None,
                      nlat=None, nlon=None, nlev=None):
    """Create a cube template with regular grid coordinates.

    Parameters
    ----------
    latlim : sequence
        bottom-left and top-right corner latitude values
    lonlim : sequence
        bottom-left and top-right corner longitude values
    nlev : int
        number of levels (vertical/z coordinate)
    nlat : int
        number of linearly spaced latitudes (y coordinate)
    nlon : int
        number of linearly spaced longitudes (x coordinate)
    """

    # get boundary information from the rotated pole grid
    rlats = rot_cube.coords('grid_latitude')[0].points
    rlons = rot_cube.coords('grid_longitude')[0].points
    ny, nx = len(rlats), len(rlons)

    # translate rotated to regular lon lat values for
    # bl_: bottom-left, tl_: top-left, br_: bottom-right, tr_: top-right,
    # tc_lat: top-centre latitude, bc_lat: bottom-centre latitude
    bl_lon, bl_lat = coord_rot2reg(rlons.min(), rlats.min())
    tl_lon, tl_lat = coord_rot2reg(rlons.min(), rlats.max())
    br_lon, br_lat = coord_rot2reg(rlons.max(), rlats.min())
    tr_lon, tr_lat = coord_rot2reg(rlons.max(), rlats.max())
    _, tc_lat = coord_rot2reg(np.median(rlons), rlats.max())
    _, bc_lat = coord_rot2reg(np.median(rlons), rlats.min())

    # regular lat lon bounds that will include all data points in rotated grid
    mm_lat = np.min([bl_lat, br_lat, bc_lat]), np.max([tl_lat, tr_lat, tc_lat])
    mm_lon = np.min([bl_lon, tl_lon]), np.max([br_lon, tr_lon])

    # update regular grid specifications
    latlim = [latlim, mm_lat][latlim is None]
    lonlim = [lonlim, mm_lon][lonlim is None]
    nlat = [nlat, ny][nlat is None]
    nlon = [nlon, nx][nlon is None]

    cs = iris.coord_systems.GeogCS(6371229.0)

    # get dimension coordinates
    for coord in rot_cube.coords(dim_coords=True):
        if coord.name() in ('model_level_number'):
            nlev = [nlev, coord.shape[0]][nlev is None]
            lev = iris.coords.DimCoord(
                np.arange(nlev),
                standard_name=coord.standard_name,
                units=coord.units)
        if coord.name() in ('latitude', 'grid_latitude'):
            lat = iris.coords.DimCoord(
                np.linspace(latlim[0], latlim[1], nlat),
                # standard_name=coord.standard_name,
                standard_name='latitude',
                units=coord.units, coord_system=cs)
        if coord.name() in ('longitude', 'grid_longitude'):
            lon = iris.coords.DimCoord(
                np.linspace(lonlim[0], lonlim[1], nlon),
                standard_name='longitude',
                units=coord.units, coord_system=cs)

    dtype = rot_cube.data.dtype
    if rot_cube.ndim == 2:
        reg_cube = iris.cube.Cube(
            np.empty((nlat, nlon), dtype),
            # np.zeros((nlat, nlon), dtype),
            dim_coords_and_dims=[(lat, 0), (lon, 1)])
    elif rot_cube.ndim == 3:
        reg_cube = iris.cube.Cube(
            np.empty((nlev, nlat, nlon), dtype),
            # np.zeros((nlev, nlat, nlon), dtype),
            dim_coords_and_dims=[(lev, 0), (lat, 1), (lon, 2)])
    else:
        raise NotImplementedError(
            'Multi-dimensional cube (>3) not supported at this time')

    return reg_cube


def rotpol2regular(rot_cube, nearest=False,
                   latlim=None, lonlim=None, nlat=None, nlon=None):
    """convert a cube from rotated pole grid to regular grid.

    Parameters
    ----------
    rot_cube : instance
        An instance of iris.cube.Cube or CubeList in rotated pole grid

    nearest : bool
        Use interpolation using nearest neighbour interpolation. Default is
        Linear interpolation.
    latlim : sequence
        Minimum and Maximum regular latitude value of output cube. Read and
        translated from the input cube if not present
    lonlim : sequence
        Minimum and Maximum regular longitude value of output cube. Read and
        translated from the input cube if not present
    nlat : int
        Number of linearly spaced latitudes. If not present, then this equals
        to the number of elements in the rotated pole grid_latitude
    nlon : int
        Number of linearly spaced longitudes. If not present, then this equals
        to the number of elements in the rotated pole grid_longitude
    """
    # Create a template cube with regular grid
    reg_cube = reg_cube_template(
        rot_cube, latlim=latlim, lonlim=lonlim, nlat=nlat, nlon=nlon)

    if nearest is False:
        scheme = iris.analysis.Linear(extrapolation_mode='mask')
    else:
        scheme = iris.analysis.Nearest(extrapolation_mode='mask')

    # TODO: Remove derived coordinate altitude from cube or add orography cube
    # to the cube
    #
    # print(dir(rot_cube.derived_coords[0]))
    # rot_cube.remove_coord('altitude')
    # rot_cube.remove_aux_factory(rot_cube.derived_coords[0])
    # rot_cube.remove_aux_factory(rot_cube.aux_factories[0])

    reg_cube = rot_cube.regrid(reg_cube, scheme)

    print('ROT CUBE: ', rot_cube)
    print('REG CUBE: ', reg_cube)

    # Re-grid cube with specific scheme
    return reg_cube


def demo_cube_conversion():
    # -------------------------------------------------------------------------
    # Demo: conversion and extraction from a rotated pole grid cube
    # Save result as netcdf and show on map
    # -------------------------------------------------------------------------

    # Read a sample UKV surface air temperature data in rotated pole grid
    pp_file = os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528.pp')
    rot_cube = iris.load_cube(pp_file)

    iris.save(
        rot_cube,
        os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528_rot.nc'),
        zlib=True
    )

    coastlines_kw = dict(resolution='10m', linewidth=0.5)

    # Plot the original cube
    qplt.pcolormesh(rot_cube)
    plt.gca().coastlines(**coastlines_kw)
    plt.gca().gridlines()
    plt.show()

    # -------------------------------------------------------------------------
    # Example 1: Extract rotated-pole subset from full domain
    # Note: the cut out will be based on bottom-left and top-right corner
    # therefore, the top corner longitude values may be missing at
    # bottom corners
    # -------------------------------------------------------------------------
    rlon_lim, rlat_lim = coord_reg2rot([-6, 2], [50, 55])
    sub_rot_cube = rot_cube.intersection(
        grid_longitude=rlon_lim, grid_latitude=rlat_lim)
    iris.save(
        sub_rot_cube,
        os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528_sub_rot.nc'),
        zlib=True
    )

    qplt.pcolormesh(sub_rot_cube)
    plt.gca().coastlines(**coastlines_kw)
    plt.show()

    # -------------------------------------------------------------------------
    # Example 2: Convert the rotated-pole cube to regular cube (full domain)
    # -------------------------------------------------------------------------
    reg_cube = rotpol2regular(rot_cube)
    iris.save(
        reg_cube,
        os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528_reg.nc'),
        zlib=True
    )

    qplt.pcolormesh(reg_cube)
    plt.gca().coastlines(**coastlines_kw)
    plt.show()

    # -------------------------------------------------------------------------
    # Example 3: Convert the rotated-pole cube to regular cube
    # (subset on-the-fly)
    # -------------------------------------------------------------------------
    sub_reg_cube = rotpol2regular(
        rot_cube, latlim=[50, 55], lonlim=[-6, 2], nlat=320, nlon=512)
    iris.save(
        sub_reg_cube,
        os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528_sub_reg.nc'),
        zlib=True
    )

    qplt.pcolormesh(sub_reg_cube)
    plt.gca().coastlines(**coastlines_kw)
    plt.show()

    # -------------------------------------------------------------------------
    # Example 4: Convert the rotated-pole cube to regular cube
    # (larger than the original domain)
    # -------------------------------------------------------------------------
    sub_reg_cube = rotpol2regular(
        rot_cube, latlim=[48, 61], lonlim=[-14, 6], nlat=260, nlon=400)
    iris.save(
        sub_reg_cube,
        os.path.expandvars('$SCRATCH/test/UKV_OP_03_20180528_sup_reg.nc'),
        zlib=True
    )

    qplt.pcolormesh(sub_reg_cube)
    plt.gca().coastlines(**coastlines_kw)
    plt.gca().gridlines(linewidth=0.5)
    plt.show()


if __name__ == '__main__':
    # demo_cube_conversion()
    pp_file = os.path.expandvars('$SCRATCH/test/20180607T0300Z.pp')
    cube = iris.load(pp_file)
    # print(cube)

    reg_cub = []
    if type(cube) is iris.cube.CubeList:
        for cub in cube:
            print(cub.name())
            print(cub)
            iris.save(
                rotpol2regular(cub),
                os.path.expandvars('$SCRATCH/test/') + cub.name() + '.nc',
                zlib=True)
            # reg_cub.append(rotpol2regular(cub))
            # reg_cub.append(rotpol2regular(
            #     cub, latlim=(50, 55), lonlim=(-6, 2), nlat=333, nlon=533)
            # )

    # print(iris.cube.CubeList(reg_cub))

    # iris.save(reg_cub, '$SCRATCH/test/UKV_OP_03_20180528_sub_rot.nc',
    #           zlib=True)

    # reg_cube = reg_cube_template(
    #     latlim=(50, 55), lonlim=(-6, 2), nlev=70, nlat=333, nlon=533
    # )
    # scheme = iris.analysis.Linear(extrapolation_mode='mask')
    # out = cube[2].regrid(reg_cube, scheme)
    # iris.save(out, '$SCRATCH/test/20180607T0300Z_3d.nc', zlib=True)

    # rlon_lim, rlat_lim = coord_reg2rot([-6, 2], [50, 55])
    # sub_rot_cube = cube_list[0].intersection(
    #     grid_longitude=rlon_lim, grid_latitude=rlat_lim
    # )
    # iris.save(
    #     sub_rot_cube,
    #     '$SCRATCH/test/UKV_OP_03_20180528_sub_rot.nc', zlib=True
    # )

    # qplt.pcolormesh(sub_rot_cube)
    # # plt.gca().coastlines(**coastlines_kw)
    # plt.show()
