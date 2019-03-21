#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Geostationary satellite imagery manipulation tool.

"""
from __future__ import (division, print_function)
import os
import numpy as np
import cartopy.crs as ccrs
from ypylib.sdf import h5Parse
from ypylib.utils import log

__version__ = "1.0"
__author__ = "Yaswant Pradhan"
__copyright__ = "(c) Crown copyright 2019, the Met Office."


# Constants
RPOL_EARTH = 6356.7523         # radius from earth centre to pole (km)
REQ_EARTH = 6378.1370          # radius from earth centre to equator (km)
rpolsq_reqsq = RPOL_EARTH ** 2 / REQ_EARTH ** 2
reqsq_rpolsq = REQ_EARTH ** 2 / RPOL_EARTH ** 2
rsqdiff_reqsq = (REQ_EARTH ** 2 - RPOL_EARTH ** 2) / REQ_EARTH ** 2

RMDI = -32768.0 * 32768.0      # Real missing data indicator
IMDI = -32768                  # Integer missing data indicator
ISDI = 32767                   # Integer to indicate space
RMDItol = 1e-3                 # RMDI tolerance

# Earth model [2]
REQ = 6378169.0                # Earth semi-major axis (m)
RPOL = 6356583.8               # Earth semi-minor axis (m)


class GeoProjection(object):
    """
    A geostationary projection definition for SPS processing

    """
    def __init__(self, satellite, channel_resolution=None,
                 apply_SEVIRI_grid_correction=False):
        """
        Defines the parameters specific to a satellite

        Args:

        * satellite (string)
            Satellite (e.g. MSG)

        Kwargs:

        * channel_resolution (string)
            Used to indicate channels with atypical resolution

        * apply_SEVIRI_grid_correction (boolean)
            True if SEVIRI grid correction to be applied (MSG
            data prior to 6 Dec 2017)

        """
        self.satellite = satellite
        self.channel_resolution = channel_resolution
        self.apply_SEVIRI_grid_correction = apply_SEVIRI_grid_correction
        self.sweep = 'y'

        if self.satellite[:3].upper() in ('MSG', 'IOD'):
            if self.channel_resolution == 'HRV':
                self.cfac = -40927014
                self.lfac = 40927014
                self.nc = 11136
                self.nl = 11136
                self.coff = 5566
                self.loff = 5566
                correction_offset = 1.5
            else:
                self.cfac = -13642337
                self.lfac = 13642337
                self.nc = 3712
                self.nl = 3712
                self.coff = 1856
                self.loff = 1856
                correction_offset = 0.5

            if self.apply_SEVIRI_grid_correction is True:
                log.info("Applying SEVIRI grid correction")
                self.coff += correction_offset
                self.loff += correction_offset

            if satellite.upper() in ('MSG_IODC', 'IODC'):
                self.sub_lon = 41.5
            elif satellite == 'MSG_RSS':
                self.sub_lon = 9.5
            else:
                self.sub_lon = 0

            # Intermediate coords deltas
            self.dx = 2**16 / self.cfac
            self.dy = 2**16 / self.lfac

            self.req = REQ
            self.rpol = RPOL
            self.satellite_height = 35785831

        elif self.satellite.upper() == 'HIM8':
            self.cfac = 20466275
            self.lfac = -20466275
            self.nc = 5500
            self.nl = 5500
            self.coff = 2750.5
            self.loff = 2750.5

            self.sub_lon = 140.7

            # Intermediate coords deltas
            self.dx = 2**16 / self.cfac
            self.dy = 2**16 / self.lfac

            # WGS84 = different from MSG
            self.req = 6378137
            self.rpol = 6356752.3
            self.satellite_height = 35785831

        elif self.satellite.upper() in ('GOES16', 'GOES17'):
            self.nc = 5424
            self.nl = 5424
            self.coff = 2712.5
            self.loff = 2712.5

            if satellite.upper() == 'GOES16':
                self.sub_lon = -75.0
            elif satellite.upper() == 'GOES17':
                self.sub_lon = -137.0

            # Intermediate coords deltas
            self.dx = np.rad2deg(0.000056)
            self.dy = np.rad2deg(-0.000056)

            # WGS84 = different from MSG
            self.req = 6378137
            self.rpol = 6356752.3
            self.satellite_height = 35786023

            # 'sweep_angle_axis' (CF conventions 1.7)
            self.sweep = 'x'

        elif self.satellite.upper() == 'GOES15':
            self.cfac = 10216334
            self.lfac = -10216334
            self.nc = 2816
            self.nl = 2816
            self.coff = 1408.5
            self.loff = 1408.5

            # GOES West
            self.sub_lon = -135.0

            # Intermediate coords deltas
            self.dx = 2**16 / self.cfac
            self.dy = 2**16 / self.lfac

            # WGS84 = different from MSG
            self.req = 6378137
            self.rpol = 6356752.3
            self.satellite_height = 35785831

    def as_cartopy_crs(self):
        """Return this projection as a Cartopy CRS instance"""
        import iris.exceptions as iexcept
        import cartopy.crs as ccrs
        proj = ccrs.Geostationary(central_longitude=self.sub_lon,
                                  satellite_height=self.satellite_height,
                                  sweep_axis=self.sweep,
                                  globe=ccrs.Globe(semimajor_axis=self.req,
                                                   semiminor_axis=self.rpol,
                                                   ellipse=None))
        # Calculate extent
        proj.extent = None  # (x0, x1, y0, y1)
        try:
            proj.extent = (
                self.iris_cube().coord('projection_x_coordinate').points[-1],
                self.iris_cube().coord('projection_x_coordinate').points[0],
                self.iris_cube().coord('projection_y_coordinate').points[0],
                self.iris_cube().coord('projection_y_coordinate').points[-1]
            )
        except iexcept.CoordinateNotFoundError as err:
            log.warning(err)

        return proj

    def forward_projection(self, lons_and_lats):
        """
        Convert lons/lats to intermediate coords and col/line vals
        (4.4.4 of [1]). Assumes lons/lats are in degrees, and lons
        in range -180 to 180.

        Method copied from Fortran:
            SpsMod_Coordinates/Sps_GeostationaryProjection.inc

        Columns/lines are 1-based (not 0-based as in Python)

        Args:

        * lons_and_lats (numpy ndarray):
            2D numpy array containing lons and lats in degrees,
            stacked 1D arrays so size of 2nd array is 2, with
            lons = lons_and_lats[:, 0], lats = lons_and_lats[:, 1]

        Returns 4-element tuple of 1D arrays: x, y, columns, lines

        """
        errmsg = ("lons_and_lats should be numpy array, 2 dims, second"
                  " of which should be of length 2")
        if (type(lons_and_lats) is not np.ndarray or
                len(lons_and_lats.shape) != 2 or
                lons_and_lats.shape[1] != 2):
            log.warning(errmsg)
            return

        lons_r = np.deg2rad(lons_and_lats[:, 0])
        lats_r = np.deg2rad(lons_and_lats[:, 1])

        # Geocentric latitude clat
        req = 1.0e-3 * self.req
        rpol = 1.0e-3 * self.rpol
        clatf = rpol * rpol / (req * req)
        h = 1.0e-3 * self.satellite_height + req
        rf = (req * req - rpol * rpol) / (req * req)

        clat = np.arctan(clatf * np.tan(lats_r))
        cos_clat = np.cos(clat)
        cos2_clat = cos_clat * cos_clat

        # Distance rl of point P from the centre of the earth
        rl = rpol / np.sqrt(1.0 - rf * cos2_clat)

        # Cartesian components of the vector rs (r1,r2,r3)
        # which points from the satellite to P.
        r1 = h - rl * cos_clat * np.cos(lons_r - np.deg2rad(self.sub_lon))
        r2 = -rl * cos_clat * np.sin(lons_r - np.deg2rad(self.sub_lon))
        r3 = rl * np.sin(clat)

        # Check to see if P is visible (the maximum viewing extent
        # is where a line from the satellite is tangent to the earth).
        rs = r1 * r1 + r2 * r2 + r3 * r3
        mask = np.logical_or.reduce((np.abs(lons_r) > np.pi,
                                     np.abs(lats_r) > np.pi / 2.0,
                                     (rs + rl * rl) > (h * h)))

        # Convert to planar coordinates
        if self.sweep == 'x':
            rn = np.sqrt(r1 * r1 + r3 * r3)
            x = np.where(mask, RMDI, np.arctan(-r2 / rn))
            y = np.where(mask, RMDI, np.arctan(r3 / r1))
        else:
            rn = np.sqrt(rs)
            x = np.where(mask, RMDI, np.arctan(-r2 / r1))
            y = np.where(mask, RMDI, np.arcsin(r3 / rn))

        # Deduce column / lines
        cols = np.where(mask, IMDI,
                        np.round(self.coff +
                                 np.rad2deg(x / self.dx))).astype('int64')
        lines = np.where(mask, IMDI,
                         np.round(self.loff +
                                  np.rad2deg(y / self.dy))).astype('int64')

        return (x, y, cols, lines)

    def inverse_projection(self, cols_and_lines=None, maxlon=None):
        """
        Generate lat, lon coords (4.4.3.2 of [1])

        Method copied from Fortran:
            SpsMod_Coordinates/Sps_GeostationaryProjection.inc

        Returns a two-element tuple containing same-shape
        1D arrays for lons and lats (in degrees).

        """
        import warnings
        warnings.simplefilter('ignore', RuntimeWarning)

        if cols_and_lines is None:
            # Do full disk
            col_start = -self.coff + 1
            col_stop = self.nc + col_start
            cols = np.arange(col_start, col_stop)

            line_start = -self.loff + 1
            line_stop = self.nl + line_start
            lines = np.arange(line_start, line_stop)

            colsm, linesm = np.meshgrid(cols, lines)
            cols = colsm.flatten()
            lines = linesm.flatten()

        else:
            errmsg = ("cols_and_lines should be numpy array, 2 dims, second"
                      " of which should be of length 2")
            if (type(cols_and_lines) is not np.ndarray or
                    len(cols_and_lines.shape) != 2 or
                    cols_and_lines.shape[1] != 2):
                log.warning(errmsg)
                return

            cols = cols_and_lines[:, 0] - self.coff + 1
            lines = cols_and_lines[:, 1] - self.loff + 1

        # Generate intermediate coords
        x = np.deg2rad(cols * self.dx)
        y = np.deg2rad(lines * self.dy)

        sin_y = np.sin(y)
        sin2_y = np.sin(y) * np.sin(y)
        cos_y = np.cos(y)
        cos2_y = np.cos(y) * np.cos(y)

        sin_x = np.sin(x)
        cos_x = np.cos(x)

        # Solve quadratic to get the vector length, but check that the
        # discriminant is positive, else the point is off the earth disk.
        req = 1e-3 * self.req
        rpol = 1e-3 * self.rpol
        latf = req * req / (rpol * rpol)
        h = 1e-3 * self.satellite_height + req

        a = cos2_y + latf * sin2_y
        b = h * cos_x * cos_y
        sd = b * b - a * (h * h - req * req)

        space_mask = sd < 0
        sn = (b - np.sqrt(sd)) / a

        # Calculate components of position vector
        s1 = h - sn * cos_x * cos_y
        if self.sweep == 'x':
            s2 = sn * sin_x
            s3 = sn * sin_y * cos_x
        else:
            s2 = sn * sin_x * cos_y
            s3 = sn * sin_y

        # Convert to latitude and longitude
        sxy = np.sqrt(s1 * s1 + s2 * s2)
        lons = np.rad2deg(np.arctan(s2 / s1) + np.deg2rad(self.sub_lon))
        lats = np.rad2deg(np.arctan(latf * s3 / sxy))

        # Adjust longitudes if -180 to 180 range needed
        if maxlon is not None:
            lons[lons >= maxlon] = (lons[lons >= maxlon] - 360.0)

        lons[space_mask] = RMDI
        lats[space_mask] = RMDI
        lons = np.ma.array(lons, mask=space_mask)
        lats = np.ma.array(lats, mask=space_mask)

        return (lons, lats)

    def iris_cube(self, general_info=None, HRV_block=None, **kwargs):
        """
        Return an Iris cube for this Geo projection.

        Add data to the cube by specifying the keyword argument 'data',
        e.g.:

        cube = GeoProjection(satellite).iris_cube(data=array)

        where 'array' is a 2D numpy array

        Data will be filled with MDI otherwise

        Note: this caters for only one block (e.g. North) for
        the SEVIRI HRV grid. A separate cube needs to be created
        for another block (e.g. South).

        Kwargs :

        * general_info (np array with 'general info' dtype)
            GeneralInfo data structure from SPS slotstore. Required
            if this cube is intended for SEVIRI HRV data.

        * HRV_block (string 'Upper' or 'Lower')
            Specify which SEVIRI HRV block to define.

        """
        import iris

        # Define the coordinate reference system for the coords
        CRS = iris.coord_systems.VerticalPerspective(
            latitude_of_projection_origin=0,
            longitude_of_projection_origin=self.sub_lon,
            perspective_point_height=self.satellite_height,
            ellipsoid=iris.coord_systems.GeogCS(semi_major_axis=self.req,
                                                semi_minor_axis=self.rpol))

        # Convert dx/dy into projection coordinates
        dx_coord = np.deg2rad(self.dx) * self.satellite_height
        dy_coord = np.deg2rad(self.dy) * self.satellite_height

        # Now calculate the coords
        if str(self.channel_resolution).upper() == 'HRV':
            if general_info is None:
                print("For a SEVIRI HRV cube supply the GeneralInfo "
                      "data to this routine")
                return

            if str(HRV_block).upper() not in ['U', 'L', 'UPPER', 'LOWER']:
                print("For a SEVIRI HRV cube specify Upper or Lower block")
                return

            if (self.satellite == 'MSG_RSS') and (HRV_block[0] != 'L'):
                print("For RSS only lower block is used, specifying 'L'...")
                HRV_block = 'L'

            if str(HRV_block).upper() in ['U', 'UPPER']:
                nx = (general_info['HRVUpperWestColumn'] -
                      general_info['HRVUpperEastColumn'] + 1)
                ny = (general_info['HRVUpperNorthLine'] -
                      general_info['HRVUpperSouthLine'] + 1)
                col_start = general_info['HRVUpperEastColumn']
                lin_start = general_info['HRVUpperSouthLine']
            else:
                nx = (general_info['HRVLowerWestColumn'] -
                      general_info['HRVLowerEastColumn'] + 1)
                ny = (general_info['HRVLowerNorthLine'] -
                      general_info['HRVLowerSouthLine'] + 1)
                col_start = general_info['HRVLowerEastColumn']
                lin_start = general_info['HRVLowerSouthLine']

            nx, ny = nx[0], ny[0]

        else:
            nx = self.nc
            ny = self.nl
            col_start = 1
            lin_start = 1

        x_start = (col_start - self.coff) * dx_coord
        y_start = (lin_start - self.loff) * dy_coord
        x_pts = x_start + np.arange(nx) * dx_coord
        y_pts = y_start + np.arange(ny) * dy_coord

        # Define the iris dimension coords
        x_coord = iris.coords.DimCoord(x_pts, units='m', coord_system=CRS,
                                       standard_name='projection_x_coordinate')
        y_coord = iris.coords.DimCoord(y_pts, units='m', coord_system=CRS,
                                       standard_name='projection_y_coordinate')

        # Build the cube. Use empty data if none provided.
        try:
            kwargs['data']
            y_inds = lin_start + np.arange(ny) - 1  # 0-indexing
            kwargs['data'] = np.ma.array(
                kwargs['data'][y_inds, :],
                mask=kwargs['data'][y_inds, :] == RMDI)
        except KeyError:
            kwargs['data'] = np.zeros((ny, nx))

        cube = iris.cube.Cube(**kwargs)
        cube.add_dim_coord(x_coord, 1)
        cube.add_dim_coord(y_coord, 0)

        return cube


def find_projections():
    """
    Cartopy projection generator

    Yields:
        generator: cartopy projections
    """
    for obj_name, o in vars(ccrs).copy().items():
        if isinstance(o, type) and issubclass(o, ccrs.Projection) and \
           not obj_name.startswith('_') and obj_name not in ['Projection']:
            yield o


def get_ccrs_list():
    """
    Get a list of all available projections in cartopy

    Returns:
        list: sorted list of available cartopy projections
    """
    return sorted([proj.__name__ for proj in find_projections()])


def print_cartopy_projections():
    """
    Print all available projections in cartopy
    """
    skip_cls = ('_', 'CRS', 'Geocentric', 'Geodetic', 'Projection',
                'RotatedGeodetic')
    for name in dir(ccrs):
        if not name.startswith(skip_cls):
            el = getattr(ccrs, name)
            if isinstance(el, type) and issubclass(el, ccrs.CRS):
                print(name)


def geo_extent(satellite='MSG'):
    """
    Get Geostationary projection extent for specific satellite
    for cartopy coordinate reference system.

    Note: For MSG, there is 1 pixel offset before/after 6 Dec 2017.

    Args:
        satellite (str, optional): Satellite ID, one from the following:
            MSG (default) - Meteosat Second Generation at 0-deg longitude.
            IODC, MSG_IODC - Meteosat Second Generation over (41.5E)
            MSG_RSS - Rapid Scan Sensor (9.5E)
            HIM8 - Himawari-8 (140.7E)
            GOES16 - Geostationary Operational Environmental Satellite (75W)
            GOES15 - Geostationary Operational Environmental Satellite (135W)
            GOES17 - Geostationary Operational Environmental Satellite (137W)


    Returns:
        tuple: (x0, x1, y0, y1) geostationary projection extent in
            X and Y direction.

    """
    sat = satellite.upper()
    return GeoProjection(sat).as_cartopy_crs().extent


def get_dust_rgb(h5handle, stride=(1, 1), satellite='MSG', **kw):
    """
    Scale channel components for dust RGB (range from 0.0 to 1.0)

    Args:
        h5handle (h5Parse object): File handler for hdf5 file
        stride (tuple, optional): number of pixels to stride in (x, y)
            direction
        satellite (str, optional): Satellite ID (default is 'MSG')

        **kw: specify band1, band2 and band3 paths in the HDF file. These
            bands should match closely with SEVIRI 8.7, 10.8 and 12 micron
            channels respectively.

    Returns:
        ndarr: Array containing Dust RGB values

    Raises:
        NotImplementedError: Raise for unknown satellite.
    """
    sat = satellite.upper()

    # Get IR brightness temp data arrays for DustRGB
    if sat == 'HIM8':
        bnd1 = kw.get('band1', '/HIM8/B11/BT')  # 8.6: 8.44--8.76
        bnd2 = kw.get('band2', '/HIM8/IR1/BT')  # 10.4: 10.3--10.6
        bnd3 = kw.get('band3', '/HIM8/IR2/BT')  # 12.3: 12.2--12.5
    elif sat in ('GOES16', 'GOES-E'):
        bnd1 = kw.get('band1', '/GOES16/Ch11/BT')  # 8.5: 8.7--8.7
        bnd2 = kw.get('band2', '/GOES16/Ch13/BT')  # 10.35: 10.1--10.6
        bnd3 = kw.get('band3', '/GOES16/Ch15/BT')  # 12.3: 11.8--12.8
    elif sat in ('GOES17', 'GOES-W'):
        bnd1 = kw.get('band1', '/GOES16/Ch11/BT')  # 8.5: 8.7--8.7
        bnd2 = kw.get('band2', '/GOES16/Ch13/BT')  # 10.35: 10.1--10.6
        bnd3 = kw.get('band3', '/GOES16/Ch15/BT')  # 12.3: 11.8--12.8
    elif sat in ('MSG', 'IODC', 'MSG_IODC'):
        bnd1 = kw.get('band1', '/MSG/IR_087/BT')  # 8.7: 8.3--9.1
        bnd2 = kw.get('band2', '/MSG/IR_108/BT')  # 10.8: 9.8--11.8
        bnd3 = kw.get('band3', '/MSG/IR_120/BT')  # 12.0: 11.0--13.0
    else:
        raise NotImplementedError('%s not implemented yet.', satellite)

    bt_870 = h5handle.get_data(bnd1)[bnd1]
    bt_108 = h5handle.get_data(bnd2)[bnd2]
    bt_120 = h5handle.get_data(bnd3)[bnd3]

    # Reverse and sub-sample arrays
    bt_870 = bt_870[::-1, ::-1][::stride[0], ::stride[1]]
    bt_108 = bt_108[::-1, ::-1][::stride[0], ::stride[1]]
    bt_120 = bt_120[::-1, ::-1][::stride[0], ::stride[1]]

    # Get RGB channels for dust/ash
    red = bt_120 - bt_108
    grn = bt_108 - bt_870
    blu = bt_108

    # Component scaling ranges for each channel:
    smin = np.array([-4.0, 0.0, 261.0])  # Valid minima for [R, G, B]
    smax = np.array([1.0, 15.0, 289.0])  # Valid maxima for [R, G, B]
    gama = np.array([1.0, 2.5, 1.0])  # Gamma values for [R, G, B]
    srng = smax - smin

    # Assign NaN to out-of-disc pixels
    invalid = np.where(bt_108 < -100.)
    red[invalid] = np.NAN
    grn[invalid] = np.NAN
    blu[invalid] = np.NAN

    # Clip RGB arrays to scaled ranges
    red = np.array(red.clip(smin[0], smax[0]))
    grn = np.array(grn.clip(smin[1], smax[1]))
    blu = np.array(blu.clip(smin[2], smax[2]))

    # Initialise RGB image
    rgbImg = np.empty((red.shape[0], red.shape[1], 3))

    # Add scaled RGB components
    with np.errstate(invalid='ignore'):
        rgbImg[:, :, 0] = ((red - smin[0]) / srng[0]) ** (1.0 / gama[0])
        rgbImg[:, :, 1] = ((grn - smin[1]) / srng[1]) ** (1.0 / gama[1])
        rgbImg[:, :, 2] = ((blu - smin[2]) / srng[2]) ** (1.0 / gama[2])

    return rgbImg


def plotx(filename, add_cbar=False,
          dataset_path=None, dust_rgb=False, dust_quality=None,
          satellite='MSG', projection=None, extent=None,
          stride=(1, 1), quick=False, MDI=None, draw_countries=False,
          coastline=True, cres='110m', ccolor='w',
          cb_width=0.02, cb_unit='', cb_nticks=None, cb_bounds=None,
          list_all=False, show_path=False, show_ticks=False,
          xglocs=None, yglocs=None, gl_font=None, glw=0.5,
          tag_fig=None, tag_col='k', title=None, figsize=(6, 6),
          save_fig=None, ax=None, **kw):
    """
    Display 2D arrays from geostationary Slotstore file.

    Args:
        filename (str): Slotstore (hdf5) filename.
        add_cbar (bool, optional): add colour bar to the plot.
        dataset_path (None, optional): Full path to 2D array data in hdf5.
        dust_rgb (bool, optional): Plot Dust RGB, require BT at 3 SEVIRI IR
            window channels in the file. Note: dataset_path and dust_rgb
            keywords are mutually exclusive.
        dust_quality (int, optional): Mask plot based on dust confidence flag
            (value range: 1-7).
        satellite (str, optional): Satellite ID ('MSG')|'GOES16'|'GOES-E'|
            'GOES17'|'GOES-W'|'IODC'|'MSG_IODC'|'HIM8'
        projection (cartopy crs, optional): Cartopy projection to use for
            mapping - one from the following
            Geostationary(), PlateCarree(), Mercator().
        extent (sequence, optional): Longitude and Latitude bounds to plot the
            data (lon0, lon1, lat0, lat1)
        stride (tuple, optional): Skip pixels in x,y dimension, Default (1, 1)
            is to show at full res.
        quick (bool, optional): quick plot. Equivalent to stride=(10, 10).
        MDI (number, optional): Missing Data Indicator.
        draw_countries (bool, optional): Show country boundaries.
        coastline (bool, optional): Show coastlines. Setting this to False
            will plot the 2D array as is (without geo-location information).
        cres (str, optional): Coastline/country boundary resolution.
        ccolor (str, optional): Coastline/boundary colour.
        cb_width (float, optional): Colour bar width.
        cb_unit (str, optional): Data unit to show on top of the colour bar.
        cb_nticks (int, optional): Number of ticks to show in the colour bar.
        cb_bounds (None, optional): Discrete bounds for colour bar.
        list_all (bool, optional): Show list of all available dataset in the
            file (no plot) and exit.
        show_path (bool, optional): Show filename and dataset path on plot.
        show_ticks (bool, optional): plot ticks and tick labels for non-mapped
            display.
        xglocs (array_like, optional): Locations of x grid-lines.
        yglocs (array_like, optional): Locations of y grid-lines.
        gl_font (dict, optional):  Font properties for grid labels default is
            {'family': 'monospace', 'size': 8, 'color': '#333333'}
        glw (float, optional): Width of grid lines.
        tag_fig (None, optional): Add optional string to top left corner of
            figure.
        tag_col (str, optional): Colour of tag_fig text.
        title (str, optional): Figure title.
        figsize (tuple, optional): Figure size. Default is (6, 6)
        save_fig (None, optional): Save plot to a file.
        ax (cartopy GeoAxis, optional): Specify axis to plot.
        **kw: accepted `imshow` keywords.

    Returns:
        matplotlib.pyplot or None: plot object or None if figure saved to disc.

    Raises:
        NotImplementedError: If Projection is not one of the following
            Geostationary, PlateCarree, Mercator
        ValueError: If dataset_path is missing or incorrect or ambiguous
    """
    import matplotlib.pyplot as plt
    import cartopy.feature as cfeature
    import cartopy.mpl.gridliner as cgrid
    from matplotlib import colors, ticker

    # defaults for MSG 0-deg service
    idx = 1  # index of datetime string in hdf filename
    c_lon = 0  # central geostationary longitude
    def_extent = (-80, 80, -80, 80)  # default map extent
    sat = satellite.upper()
    if sat in ('IODC', 'MSG_IODC'):
        # update defaults for IODC service
        idx = 2
        c_lon = 41.5
        def_extent = (-40, 120, -80, 80)
    if sat == 'HIM8':
        # update defaults for Himawari service
        c_lon = 140.7
        def_extent = (60, 220, -80, 80)
    if sat in ('GOES16', 'GOES-E'):
        c_lon = -75
        def_extent = (-155, 5, -80, 80)
    if sat == 'GOES15':
        c_lon = -135
        # def_extent = (-155, 5, -80, 80)
    if sat in ('GOES17', 'GOES-W'):
        c_lon = -137
        # def_extent = (-155, 5, -80, 80)

    h5 = h5Parse(filename)
    if list_all:
        h5.ls()
        return

    if quick:
        stride = (10, 10)

    apply_dust_quality = False
    if dust_quality:
        dconf = '/Product/GM/DustAOD/DustConfidence'
        try:
            dset = h5.get_data(dconf)
            dust_qual = dset[dconf][::-1, ::-1][::stride[0], ::stride[1]]
            apply_dust_quality = True
        except KeyError:
            log.warn('Dust_quality not available - filter wont be applied')

    if gl_font is None:
        gl_font = dict(family='monospace', size=8, color='#333333')

    # Get dust RGB array
    if dust_rgb and dataset_path is None:
        dataset_path = 'DustRGB'
        add_cbar = False
        plot_array = get_dust_rgb(h5, stride=stride, satellite=satellite)

    elif dataset_path and dust_rgb is False:
        dset = h5.get_data(dataset_path)
        if len(list(dset.keys())) == 0:
            return
        plot_array = dset[dataset_path][::-1, ::-1][::stride[0], ::stride[1]]
        # print plot_array.dtype

        if plot_array.dtype in (
                np.int8, np.int16, np.intc, np.uint8, np.uint16, np.intp):
            MDI = [MDI, IMDI][MDI is None]
            plot_array = np.ma.masked_equal(plot_array, MDI)
            if apply_dust_quality:
                plot_array = np.ma.masked_less(dust_qual, dust_quality)
        else:
            MDI = [MDI, RMDI][MDI is None]
            plot_array[plot_array <= (MDI + RMDItol)] = np.nan
            if apply_dust_quality:
                plot_array = np.ma.where(
                    dust_qual >= dust_quality, plot_array, np.nan)

    elif dataset_path and dust_rgb:
        raise ValueError('dust_rgb must be False if dataset_path is given')
    else:
        raise ValueError(
            'Either dust_rgb or dataset_path must be specified.\n\n'
            '  msgview(filename, [dataset_path|dust_rgb=True], ..) OR\n'
            '  msgview(filename, list_all=True) to see all datasets.\n')

    # Get MSG datetime string from hdf5 filename
    msg_datetime = os.path.basename(filename).split('_')[idx]
    title_str = [msg_datetime, os.path.basename(dataset_path)]

    # Reassign 3-channel DustMask values for readable legend/title
    if title_str[1] in 'Dust':
        title_str[1] = 'DustMask'
        plot_array[plot_array < -0.5] = 0
        plot_array[plot_array == 300] = 1
        plot_array[plot_array == 400] = 2
        plot_array[plot_array == 500] = 3

    # Rename DustHeightError to DustPressureError
    if title_str[1] in 'DustHeightError':
        title_str[1] = 'DustPressureError'

    if title is None:
        title = '{}  {}'.format(title_str[0], title_str[1])

    # Start plot
    if ax:
        projection = ax.projection
    else:
        plt.figure(figsize=figsize)
        projection = [projection, ccrs.Geostationary()][projection is None]

    # Normalise colour based on discrete colour bounds (cb_bounds)?
    if cb_bounds:
        bounds = np.array(cb_bounds)
        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds))
    else:
        norm = None

    # combine plot_array and imshow keywords
    plot_array_kw = {
        'extent': geo_extent(satellite=satellite),
        'origin': 'upper',
        'norm': None,
        'interpolation': None}
    plot_array_kw.update(**kw)

    if coastline:
        if projection in (ccrs.PlateCarree(),
                          ccrs.PlateCarree(central_longitude=0),
                          ccrs.PlateCarree(central_longitude=180),
                          ccrs.Mercator(),
                          ccrs.Mercator(central_longitude=0),
                          ccrs.Mercator(central_longitude=180)):

            # Plot using rectangular projection
            if ax is None:
                ax = plt.axes(projection=projection)
            extent = [extent, def_extent][extent is None]
            ax.set_extent(extent, ccrs.PlateCarree())
            gl = ax.gridlines(
                draw_labels=True, alpha=0.5, xlocs=xglocs, ylocs=yglocs,
                linestyle='-', linewidth=glw)
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xformatter = cgrid.LONGITUDE_FORMATTER
            gl.yformatter = cgrid.LATITUDE_FORMATTER
            gl.xlabel_style, gl.ylabel_style = gl_font, gl_font
            im = ax.imshow(plot_array, transform=ccrs.Geostationary(c_lon),
                           **plot_array_kw)

        elif projection == ccrs.Geostationary():
            # Plot in Geostationary projection
            if ax is None:
                ax = plt.axes(projection=ccrs.Geostationary(c_lon))
            ax.set_global()
            im = ax.imshow(plot_array, **plot_array_kw)

            if extent:
                ax.set_extent(extent, ccrs.PlateCarree())
        else:
            raise NotImplementedError(
                '{} projection not implemented.\n'
                'Please use one from the following:\n'
                '\tGeostationary\n\tPlateCarree or\n\tMercator\n'.format(
                    projection.__class__))

        ax.coastlines(resolution=cres, lw=0.5, color=ccolor)
    else:
        if extent:
            log.warn('Ignored extent. For subsets please set coastline=True')
        ax_pos = [0.05, 0.1, 0.8, 0.8]
        if show_ticks:
            ax_pos = [0.1, 0.1, 0.8, 0.8]
        if ax is None:
            ax = plt.axes(ax_pos)

    # Add country borders
    if draw_countries:
        ax.add_feature(cfeature.BORDERS, lw=0.4, edgecolor=ccolor)

    pos = ax.get_position()

    # Tag a figure with additional string at top-left (useful for figure
    # panel numbering)
    if tag_fig:
        # title = '{} {}'.format(tag_fig, title)
        ax.text(0.01, 0.95, tag_fig, color=tag_col, fontsize=14,
                transform=ax.transAxes)

    # Draw title text
    ax.set_title(title, fontsize="13", loc='left')

    # Hide plot ticks and tick labels?
    if show_ticks is False:
        for ticx, ticy in list(
                zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks())):
            ticx.tick1On = ticx.tick2On = False
            ticx.label1On = ticx.label2On = False
            ticy.tick1On = ticy.tick2On = False
            ticy.label1On = ticy.label2On = False

    # Show image without coastline - simple 2D array imshow()
    if coastline is False:
        im = plt.imshow(plot_array, norm=norm, interpolation=None, **kw)
        ax.grid(lw=0.5, alpha=0.5)

    # Attach vertical colour bar
    if add_cbar:
        cax = plt.axes([pos.x1, pos.y0, cb_width, pos.height])
        cb = plt.colorbar(mappable=im, cax=cax, format='%g')
        cb.ax.tick_params(direction='in')
        cb.ax.set_title('{}'.format(cb_unit))
        if cb_nticks:
            tick_locator = ticker.MaxNLocator(nbins=cb_nticks)
            cb.locator = tick_locator
            cb.update_ticks()
        if title_str[1] in 'DustMask':
            cb.ax.set_yticklabels(
                ['Dust', 'Ice cloud', 'Low cloud/\nSurface', 'Uncertain'],
                rotation='90', va='bottom')

    if show_path:
        # Show filename and HDF5 dataset path on plot
        textstr = '{}\n{}'.format(dataset_path, filename)
        plt.gcf().text(0.99, 0.01, textstr, fontsize=8, ha='right')

    if save_fig:
        # save figure as file to disk
        plt.savefig(save_fig)
        plt.close()
        log.info('%s saved.', save_fig)
    else:
        # return pyplot object
        return plt


class MSG(object):

    """
    Warning! The software is for use with MSG data only and will not work in
    the given implementation for Meteosat first generation data.

    Note on CFAC/LFAC and COFF/LOFF:
        The parameters CFAC/LFAC and COFF/LOFF are the scaling coefficients
        provided by the navigation record of the LRIT/HRIT header and used
        by the scaling function given in Ref [1], page 28.

        COFF/LOFF are the offsets for column and line which are basically 1856
        and 1856 for the VIS/IR channels and refer to the middle of the image
        (centre pixel). The values regarding the High Resolution Visible
        Channel (HRVis) will be made available in a later issue of this
        software.

        CFAC/LFAC are responsible for the image "spread" in the NS and EW
        directions. They are calculated as follows:

        CFAC = LFAC = 2**16 / delta, with
        delta = 83.84333 micro Radian (size of one VIS/IR MSG pixel)

        CFAC = LFAC =  781648343.404  1/rad for VIS/IR
        which should be rounded to the nearest integer as stated in Ref [1].

        CFAC = LFAC = 781648343  1/rad for VIS/IR

        The sign of CFAC/LFAC gives the orientation of the image.
        Negative sign give data scanned from south to north as in the
        operational scanning. Positive sign vice versa.

    Reference:
        [1] LRIT/HRIT Global Specification (CGMS, Issue 2.8, 30.10.2013)
        for the parameters used in the program. section 4.3.3.2
        http://www.cgms-info.org/documents/
        cgms-lrit-hrit-global-specification-%28v2-8-of-30-oct-2013%29.pdf

    Attributes:
    Meteosat Second Generation specific values/constants.
        CFAC (float): Scaling coefficients (see note above)
        COFF (float): Column offset (see note above)
        filename (TYPE): MSG slotstore (H5) filename
        LFAC (float): Scaling coefficients (see note above)
        LOFF (float): Line offset (see note above)
        SAT_ALT (float): Distance from earth centre to satellite
        satellite (str): Satellite ID.
        sd_coeff (TYPE): 1737122264.
        SUB_LON (float): Longitude of sub-satellite point
    """

    satellite = 'MSG'

    def __init__(self):
        self.SAT_ALT = 42164.0
        self.SUB_LON = 0.0
        self.COFF = 1856.0
        self.LOFF = 1856.0
        self.CFAC = -781648343.0
        self.LFAC = -781648343.0
        self.sd_coeff = self.SAT_ALT ** 2 - REQ_EARTH ** 2
        self.filename = None

    def geo2pix(self, lon, lat):
        """
        Returns the pixel column and row number of an MSG image for a given
        pair of longitude, latitude values.

        Note: calculation based on the formulae given in Ref [1]

        Args:
            lon (array_like): Longitude value(s)
            lat (arraye_like): Latitude value(s)

        Returns:
            tuple: Tuple of (col, row) where col and row are masked arrays

        Examples:
        >>> x, y = MSG().geo2pix(0, 0)
        >>> print(x, y)
        1856 1856
        """
        lonR = np.radians(lon)
        latR = np.radians(lat)

        # Calculate the geocentric latitude from the geographic one using
        # equations on page 24, Ref [1]
        #     c_lat = np.arctan(0.993243 * np.tan(latR))

        c_lat = np.arctan(rpolsq_reqsq * np.tan(latR))

        # Using c_lat calculate the length form the earth centre to the
        # surface of the earth ellipsoid ;equations on page 24, Ref. [1]
        #     re = rpol / (sqrt(1 - (req^2 - rpol^2) / req^2 * cos^2(c_lat))

        rl = RPOL_EARTH / np.sqrt(1 - rsqdiff_reqsq * np.cos(c_lat) ** 2)

        # Calculate the forward projection using equations on page 24, Ref. [1]
        r1 = self.SAT_ALT - rl * np.cos(c_lat) * np.cos(lonR - self.SUB_LON)
        r2 = -rl * np.cos(c_lat) * np.sin(lonR - self.SUB_LON)
        r3 = rl * np.sin(c_lat)
        rn = np.sqrt(r1 ** 2 + r2 ** 2 + r3 * 2)

        # Check for visibility, whether the point on the earth given by the
        # latitude/longitude pair is visible from the satellite or not. This
        # is given by the dot product between the vectors of:
        #  1) the point to the spacecraft,
        #  2) the point to the centre of the earth.
        # If the dot product is positive the point is visible otherwise it is
        # invisible.

        dotprod = (r1 * (rl * np.cos(c_lat) * np.cos(lonR - self.SUB_LON)) -
                   r2 ** 2 - r3 ** 2 * ((REQ_EARTH / RPOL_EARTH) ** 2))

        # The forward projection is x and y
        xx = np.arctan(-r2 / r1)
        yy = np.arcsin(-r3 / rn)

        # Convert to pixel column and row using the scaling functions on
        # $4.4.4 page 32, Ref. [1]. And finding nearest integer value for them.
        col = self.COFF + np.rint(xx * 2 ** -16 * self.CFAC)
        row = self.LOFF + np.rint(yy * 2 ** -16 * self.LFAC)
        col = np.ma.masked_where(dotprod < 0, col)
        row = np.ma.masked_where(dotprod < 0, row)

        return col.astype(int), row.astype(int)

    def pix2geo(self, col, row, dtype='float64'):
        """
        Calculate the longitude and latitude value of a pixel on MSG disc
        given the column and row index of that pixel.

        Args:
            col (array_like): Column index of MSG/SEVIRI pixel
            row (array_like): Row index of MSG/SEVIRI pixel
            dtype (str, optional): Data type precision for lat and lon arrays
                (default 'float64')

        Returns:
            tuple: retruns a tuple of (lon, lat) where lon and lat are masked
            arrays

        Examples:
            >>> lon, lat = MSG().pix2geo(1856, 1856)
            >>> print(lon, lat)
            0.0 0.0
        """

        # Calculate viewing angle of the satellite by inverting equations
        # on page 32, Ref [1].
        x = (np.array(col) - self.COFF) * 2 ** 16 / self.CFAC
        y = (np.array(row) - self.LOFF) * 2 ** 16 / self.LFAC

        # Calculate the inverse projection - formulas using equations on
        # page 29, Ref. [1]
        # First check for visibility, whether the pixel is located on the
        # earth surface or in space. To do this calculate the argument to
        # sqrt of "sd", which is named "sa". If it is negative then the
        # sqrt will return NaN and the pixel will be located in space,
        # otherwise all is fine and the pixel is located on the earth surface.

        sa = (((self.SAT_ALT * np.cos(x) * np.cos(y)) ** 2) -
              (np.cos(y) ** 2 + reqsq_rpolsq * np.sin(y) ** 2) *
              self.sd_coeff)
        sd = np.sqrt(sa)
        sn = ((self.SAT_ALT * np.cos(x) * np.cos(y) - sd) /
              (np.cos(y) ** 2 + reqsq_rpolsq * np.sin(y) ** 2))
        s1 = self.SAT_ALT - sn * np.cos(x) * np.cos(y)
        s2 = sn * np.sin(x) * np.cos(y)
        s3 = -sn * np.sin(y)
        sxy = np.sqrt(s1 ** 2 + s2 ** 2)

        # Using the previous calculations the inverse projection can be
        # calculated now, which means calculating the lat/long from
        # the pixel row and column by equations on page 29, Ref [1].

        lon = np.arctan(s2 / s1 + self.SUB_LON)
        lat = np.arctan((reqsq_rpolsq * s3) / sxy)

        # Mask invisible (outside SEVIRI field-of-view) pixel locations
        lon = np.ma.masked_where(sa < 0, lon)
        lat = np.ma.masked_where(sa < 0, lat)

        # Convert from radians into degrees
        return np.degrees(lon).astype(dtype), np.degrees(lat).astype(dtype)

    def fulldisc2xyz(self, fd_array, missing_value=RMDI, vname='data',
                     xydtype='float64'):
        """
        Extract valid pixels from a SEVIRI full-disc array and return a
        dictionary of valid data and geographic coordinates for corresponding
        pixel column and rows.


        Args:
            fd_array (2d_array): SEVIRI full disc array of shape (3712, 3712)
            missing_value (int or float, optional): values <= missing_value to
                be remove (default None)
            vname (str, optional): Variable name (default 'data')
            xydtype (str, optional): Data type precision for x and y
                arrays (default 'float64')

        Returns:
            dict: dictionary with valid data points and associated geographic
            coordinates
        """
        w = np.where(fd_array != missing_value)  # returns tuple (rows, cols)
        lon, lat = self.pix2geo(w[1], w[0], dtype=xydtype)
        # print lon.data.dtype, lat.data.dtype
        return {vname: fd_array[w], 'lon': lon.data, 'lat': lat.data}

    def extract(self, filename, dsname, stride=(1, 1)):
        """
        Extract data from a slotstore dataset.

        Args:
            filename (str): Slotstore filename
            dsname (str): H% dataset path
            stride (tuple, optional): stride number of plixels in x and y.

        Returns:
            ndarray: Extracted region
        """
        self.filename = filename
        d = h5Parse(filename)
        data = d.get_data(dsname)[::-1, ::-1][::stride[1], ::stride[0]]
        return data


def wrap_lon(longitude, hemisphere=False):
    """
    Converts Longitude values in an array range from [-180:180] to [0:360]
    (sphere) or vice-versa (hemisphere).

    Args:
        longitude (num or array_like): Longitude array.
        hemisphere (bool, optional): Setting this keyword assumes original
            array ranges from 0 to 360 and therefore returned array will be in
            -180 to 180 range.

    Returns:
        num or array_like: transformed longitude value(s).

    Example:
        >>> print(wrap_lon(-190))
        170
        >>> print(wrap_lon([-170.00, 180, 300]))
        [190. 180. 300.]
        >>> print(wrap_lon([190.00, 300], hemisphere=True))
        [-170.  -60.]
    """
    longitude = np.asarray(longitude)
    if hemisphere:
        wrapped_lon = ((longitude + 180) % 360) - 180  # range 0:360
    else:
        wrapped_lon = (longitude + 360) % 360  # range -180:+180
    return wrapped_lon


if __name__ == "__main__":
    # import doctest
    # doctest.testmod()
    # pass

    filename0 = os.path.expandvars(
        '$SCRATCH/MSG/LiteSlotstore/ICE-D/20150812/'
        'MSG_201508121200_DustAOD.h5')
    plt_kw = dict(
        dust_rgb=True,
        # extent=(-60.5, 55, 0, 70.5),
        # stride=(8, 8),
        extent=(-20, 20, 10, 40),
        draw_countries=True,
        ccolor='k',
        # projection=ccrs.PlateCarree(),
        xglocs=list(range(-60, 55, 20)),
        yglocs=list(range(10, 71, 20)),
        # save_fig=os.path.expandvars('$SCRATCH/dustrgb.png'),
    )
    plotx(filename0, **plt_kw).show()

    # print(wrap_lon([200, 300], hemisphere=True))
    # import doctest
    # doctest.testmod()

    # lon_0 = 180
    # lon = -180 + 5 * np.arange(72)
    # lon_shift = shiftlon(lon, lon_0)
    # print "original lon:", lon
    # print " shifted lon:", lon_shift
    # h5f = '/sample_data/slotstores/MSG_201202151600_lite.h5'
    # msg = MSG()
    # f = msg.geo2pix([0, 2], [0, 10])
    # print f[0], f[1]
    # print(msg.pix2geo([1856, 600, 100], [1857, 600, 100]))
    # x = msg.geo2pix([-20, 0, 10, 20], [0, 0, 0, 0])
    # print(x[0]._data, x[1]._data)
    # print h.extract(h5f, 'MSG/Ch10/Raw', stride=(4, 4)).shape
