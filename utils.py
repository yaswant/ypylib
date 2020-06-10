#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
General purpose utilities.

"""
from __future__ import print_function
import os
import csv
import glob
import errno
import math
import socket
import inspect
import logging as log
import numpy as np
import matplotlib as mpl
import matplotlib.cm as mpl_cm
import collections
from sys import platform
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from scipy.stats import describe, binned_statistic_2d
from matplotlib import pyplot as plt
try:
    from Tkinter import Tk
    from tkFileDialog import askopenfilenames
except ImportError:
    from tkinter import Tk
    from tkinter.filedialog import askopenfilenames
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    Basemap = None
try:
    import cartopy
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
except ImportError:
    cartopy = None
log.basicConfig(
    format='%(asctime)s %(levelname)s: %(message)s',
    datefmt='%F %T', level=log.INFO)
# To log debug messages, call
# >>> log.getLogger().setLevel(log.DEBUG)

__version__ = "1.0"
__author__ = "Yaswant Pradhan"


class Convert(object):
    """
    Simple unit converter

    Attributes
    ----------
    value : number
        Input value to convert.

    """

    def __init__(self, value):
        self.value = np.array(value)

    def __enter__(self):
        return self

    def ft2m(self):
        return self.value * 0.3048

    def ft2km(self):
        """Feet to Kilometre"""
        return self.value * 0.3048 / 1e3

    def kt2mph(self):
        """Knotts to miles-per-hour"""
        return self.value * 1.15077945

    def mph2kt(self):
        """Miles-per-hour to Knotts"""
        return self.value / 1.15077945

    def kg2lb(self):
        """Kilograms to Pounds"""
        return self.value * 2.20462262

    def lb2kg(self):
        """Pounds to Kilograms"""
        return self.value / 2.20462262

    def f2c(self):
        """Fahrenheits to Celcius"""
        return (self.value - 32) / 1.8

    def c2f(self):
        """Celcius to Fahrenheits"""
        return self.value * 1.8 + 32

    def solzen2airmass(self):
        """Solar zenith angle (degrees) to Airmass.

        Air Mass is the path length which light takes through the atmosphere
        normalized to the shortest possible path length (ie, when the sun is
        directly overhead). The Air Mass quantifies the reduction in the power
        of light as it passes through the atmosphere and is absorbed by air
        and dust. For a flat horizontal atmosphere, Airmass is defined as

            1/cos(sol_zen)

        But because of the curvature of the atmosphere, the air mass is not
        quite equal to the atmospheric path length when the sun is close to
        the horizon. At sunrise, the angle of the sun from the vertical
        position is 90 deg and the air mass is infinite, whereas the path
        length clearly is not. An equation which incorporates the curvature of
        the earth is given by F. Kasten and Young, A. T. (1989)
        "Revised optical air mass tables and approximation formula",
        Applied Optics, vol. 28, pp. 4735-4738, 1989.

            1/[cos(sol_zen) + 0.50572(96.07995-sol_zen)^-1.6364]
        """
        corrected_cosine = np.cos(np.radians(self.value)) + \
            0.50572 * (96.07995 - self.value)**-1.6364
        return 1 / corrected_cosine

    def __exit__(self, exc_type, exc_value, traceback):
        self.value = 0


class Integer(object):
    """The usual single-bit operations will work on any Python integer."""

    def __init__(self, int_type):
        super(Integer, self).__init__()
        self.int_type = int_type

    def test_bit(self, offset, mask_type=None):
        """test_bit() returns a non-zero result, 2**offset, if the bit at
        'offset' is one.

        Parameters
        ----------
        offset : int
            Bit position to test.  It is up to the user to make sure that the
            value of offset makes sense in the context of the program
        mask_type : str, optional
            Returns a boolean (True|False) mask, if set to 'bool' or
            a binary (1|0) mask if set to 'bin' or 'int', instead of 2**offset.

        Returns
        -------
        int_type or boolean
            Description

        Examples
        --------
        test 3rd bit for 10

        >>> print(Integer(10).test_bit(3))
        8

        >>> print(Integer(10).test_bit(3, mask_type='bool'))
        True

        >>> print(Integer(10).test_bit(3, mask_type='int'))
        1
        """
        mask = 1 << offset
        if mask_type == 'bool':
            return(self.int_type & mask != 0)
        elif mask_type in ('int', 'bin'):
            return(self.int_type & mask != 0) * 1
        else:
            return(self.int_type & mask)

    def set_bit(self, offset):
        """set_bit() returns an integer with the bit at 'offset' set to 1."""
        mask = 1 << offset
        return(self.int_type | mask)

    def clear_bit(self, offset):
        """clear_bit() returns an integer with the bit at 'offset' cleared."""
        mask = ~(1 << offset)
        return(self.int_type & mask)

    def toggle_bit(self, offset):
        """toggle_bit() returns an integer with the bit at 'offset' inverted,
        0 -> 1 and 1 -> 0.
        """
        mask = 1 << offset
        return(self.int_type ^ mask)


class List(list):
    """Create a custom List object to which user attributes can be added.

    Parameters
    ----------
    list : list, optional
        input list object

    Returns
    -------
    custom List object with user defined attributes

    Example
    -------
        >>> a_list = List()
        >>> a_list.name = 'ListName'
        >>> a_list.append(np.array(10))

        Or,
        >>> a_list = List([np.array(10)], name='ListName')
    """

    def __new__(self, *args, **kwargs):
        return super(List, self).__new__(self, args, kwargs)

    def __init__(self, *args, **kwargs):
        if len(args) == 1 and hasattr(args[0], '__iter__'):
            list.__init__(self, args[0])
        else:
            list.__init__(self, args)
        self.__dict__.update(kwargs)

    def __call__(self, **kwargs):
        self.__dict__.update(kwargs)
        return self


class XYZ(object):
    """Discrete triplet data (x, y, z) analyser."""

    def __init__(self, x, y, z, wrap_lon=False):
        """XYZ Constructor.

        Parameters
        ----------
        x : array_like, shape(N,)
            An array containing the x coordinates of the points to be binned
        y : array_like, shape(N,)
            An array containing the y coordinates of the points to be binned
        z : array_like, shape(N,) f(x,y)
            actual data to be re-sampled (average at each grid cell)
        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)
        if wrap_lon:
            self.x = ((self.x + 180) % 360) - 180
        # bin parameters for griddata
        self.delta = (1, 1)
        self.limit = np.array([[-180, 180], [-90, 90]])
        # plot parameters used in mapdata
        self.figsize = (8, 5)
        self.figheight = self.figsize[1]
        self.gspacing = (30, 30)
        self.xlab = [0, 0, 0, 1]
        self.ylab = [1, 0, 0, 0]
        self.cbpad = '10%'
        self.G, self.xc, self.yc = None, None, None

        # update xlimit based on actual data
        if self.x.max() > 180:
            self.limit[0] = [0, 360]

    def __enter__(self):
        return self

    def _get_extent(self):
        """Get extent for x and y values."""
        return [[self.x.min(), self.x.max()], [self.y.min(), self.y.max()]]

    def _get_latts(self, spacing):
        """Get latitude grid mark locations."""
        return np.arange(-90., 99., spacing)

    def _get_lonts(self, spacing):
        """Get longitude grid mark locations."""
        return np.arange(-180., 181., spacing)

    def bindata(self, **kw):
        """Bin irregular 1D data (triplets) on to 2D plane.

        Deprecated, use griddata with 'mean' method instead.

        Parameters
        ----------
        delta : sequence or [float, float], optional
            Output grid resolution in x and y direction (default (1, 1)).
            The delta specification:
             *  If scalar, the grid resolution for the two dimensions are equal
                (dx=dx=delta)
             *  If [float, float], the grid resolution for the two dimensions
                (dx, dy = delta)
        limit : [[float, float], [float, float]], optional
            Output domain limit [[x0,x1], [y0,y1]] (default calculated from
            actual x and y range)
        globe : bool, optional
            If True, sets the grid x and y limit to [-180,180] and [-90,90],
            respectively. If False, grid x and y limits are taken from input.
            (default False)
        order : bool, optional
            If True, returns a upside-down flip and rotated array (default
            False)

        Returns
        -------
        G : MaskedArray, shape(nxc,nyc)
            The 2-dimensional binned (averaged) array of z
        xc : ndarray, shape(nx,)
            The bin centres along the x dimension
        yc : ndarray, shape(ny,)
            The bin centres along the y dimension

        Example
        -------
        >>> from numpy.random import normal
        >>> from ypylib.utils import XYZ
        >>> x = normal(3, 1, 100)
        >>> y = normal(1, 1, 100)
        >>> z = x * y
        >>> G, xc, yc = XYZ(x, y, z).bindata(delta=[0.1, 0.1])
        """
        delta = kw.get('delta', self.delta)
        globe = kw.get('globe', False)
        order = kw.get('order', False)
        limit = kw.get('limit', self._get_extent())

        if globe is True:
            limit = self.limit
        if len(delta) == 1:
            delta = [delta, delta]
        x, y, z = np.asarray(self.x), np.asarray(self.y), np.asarray(self.z)
        # construct x and y bins
        xs = np.arange(limit[0][0], limit[0][1] + delta[0], delta[0])
        ys = np.arange(limit[1][0], limit[1][1] + delta[1], delta[1])
        # sum and edges of each bin
        Hv, xl, yb = np.histogram2d(x, y, bins=[xs, ys], weights=z)
        # shift bin edged by 0.5*delta to get centre of bins
        xc, yc = xl[:-1] + delta[0] / 2., yb[:-1] + delta[1] / 2.
        # counts for each bin
        Hn, _, _ = np.histogram2d(x, y, bins=[xs, ys])
        # mask sum array where count = 0 i.e. no data in the bin
        Hv = np.ma.masked_where(Hn == 0, Hv)
        if order is True:
            # rotate grid (column major) for display purpose
            return np.flipud(np.rot90(Hv / Hn)), xc, yc
        else:
            return Hv / Hn, xc, yc

    def griddata(self, stat='mean', **kw):
        """Compute a bi-dimensional binned statistic for one or more sets of
        data.

        This is a generalization of a histogram2d function implemented in
        bindata. A histogram divides the space into bins, and returns the
        count of the number of points in each bin. This function allows the
        computation of the sum, mean, median, or other statistic of the values
        (or set of values) within each bin.

        Parameters
        ----------
        delta : sequence or [float, float], optional
            Output grid resolution in x and y direction (default (1, 1)).
            The delta specification:
             *  If scalar, the grid resolution for the two dimensions are equal
                (dx=dx=delta)
             *  If [float, float], the grid resolution for the two dimensions
                (dx, dy = delta)
        globe : bool, optional
            If True, sets the grid x and y limit to [-180,180] and [-90,90],
            respectively. If False, grid x and y limits are taken from input.
            (default False)
        limit : [[float, float], [float, float]], optional
            Output domain limit [[x0,x1], [y0,y1]] (default calculated from
            actual x and y range)
        stat : str, optional
            The statistic to compute (default is 'mean'). Available statistics
            are: 'mean', 'std', median', 'count', 'sum', 'min', 'max'
        norm : bool, optional
            Normalize (rescale values between 0 and 1) gridded output values.
        order : bool, optional
            If True, returns a upside-down flip and rotated array (default
            False)

        Returns
        -------
        G : ndarray, shape(nxc, nyc)
            The 2-dimensional binned (averaged) array of z
        xc : ndarray, shape(nx,)
            The bin centres along the x dimension
        yc : ndarray, shape(ny,)
            The bin centres along the y dimension

        See also
        --------
        bindata (deprecated)
            equivalent to stat='mean'
        """
        try:
            from collections.abc import Sequence
        except ImportError:
            from collections import Sequence

        delta = kw.get('delta', self.delta)
        globe = kw.get('globe', False)
        order = kw.get('order', False)
        limit = kw.get('limit', self._get_extent())
        norm = kw.get('norm', False)
        stat = kw.get('method', stat)
        # limit = np.asarray(limit)

        # parse stat parameter
        stat = stat.lower()
        stat = 'sum' if stat in ('sum', 'total') else stat
        stat = 'mean' if stat in ('avg', 'average', 'mean') else stat
        stat = 'count' if stat in ('count', 'frequency') else stat

        # if stat in ('sum', 'total'):
        #     stat = 'sum'
        # if stat in ('avg', 'average', 'mean'):
        #     stat = 'mean'
        # if stat in ('count', 'freq', 'frequency'):
        #     stat = 'count'

        if globe is True:
            limit = self.limit
        if isinstance(delta, Sequence) is False:
            delta = [delta, delta]
        if len(delta) == 1:
            delta = [delta, delta]

        # construct x and y bins
        xs = np.arange(limit[0][0], limit[0][1] + delta[0], delta[0])
        ys = np.arange(limit[1][0], limit[1][1] + delta[1], delta[1])

        stat4 = binned_statistic_2d(
            self.x, self.y, self.z, statistic=stat, bins=[xs, ys])

        # bin edges to bin centres
        xc = (stat4[1][1:] + stat4[1][:-1]) / 2
        yc = (stat4[2][1:] + stat4[2][:-1]) / 2
        # print(xc.min(), xc.max(), yc.min(), yc.max())

        # mask sum, std and count arrays when no data available in the bin
        if stat in ('count', 'std', 'sum'):
            count4 = binned_statistic_2d(
                self.x, self.y, self.z, statistic='count', bins=[xs, ys])
            G = np.ma.masked_where(count4[0] == 0, stat4[0])
        else:
            G = np.ma.masked_invalid(stat4[0])

        # normalise gridded data if requested
        if norm:
            gmin, gmax = G.min(), G.max()
            G = (G - gmin) / (gmax - gmin)

        self.G, self.xc, self.yc = G, xc, yc

        # transpose grid if requested
        if order:
            self.G = G.transpose()
            return self.G, self.xc, self.yc
        else:
            return self.G, self.xc, self.yc

    def shift_midpoints(self, **kw):
        """Shift x and y vectors by half-grid for pcolormesh

        ``pcolormesh`` patches bottom-left is usually aligned with grid centre
        (x,y) which makes the plot shift by half-grid in top-right direction.
        This method applies the correction so that the patch centres are
        aligned with the grid centre.
        """

        # - expand/and shift x grid centre
        xe = np.zeros(self.xc.size + 2)
        xe[1:-1] = self.xc
        xe[0] = self.xc[0] - np.diff(self.xc)[0]
        xe[-1] = self.xc[-1] + np.diff(self.xc)[-1]
        xc = xe[:-1] + 0.5 * (np.diff(xe))

        # - expand/and shift y grid centre
        ye = np.zeros(self.yc.size + 2)
        ye[1:-1] = self.yc
        ye[0] = self.yc[0] - np.diff(self.yc)[0]
        ye[-1] = self.yc[-1] + np.diff(self.yc)[-1]
        yc = ye[:-1] + 0.5 * (np.diff(ye))
        return xc, yc

    def plot(self, **kw):
        """Line plot (x, y, z) as series."""
        _f, ax = plt.subplots(1, figsize=self.figsize)
        ax.plot(self.x.ravel(), c='r', label='x', **kw)
        ax.plot(self.y.ravel(), c='g', label='y', **kw)
        ax.plot(self.z.ravel(), c='b', label='z', **kw)
        plt.legend(ncol=3, borderaxespad=0)
        plt.tight_layout()
        return plt

    def scatter(self, **kw):
        """Scatter plot (x, y, z) data triplets."""
        _f, ax = plt.subplots(1, figsize=self.figsize)
        im = ax.scatter(
            self.x, self.y, c=self.z, lw=0, rasterized=True, **kw)
        plt.title('xyz.scatter()')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.colorbar(im)
        plt.tight_layout()
        return plt

    def hexbin(self, **kw):
        """Hexagon plot (x, y, z) data triplets.

        hexbin is a bivariate histogram in which the xy-plane is tessellated
        by a regular grid of hexagons.
        """
        _f, ax = plt.subplots(1, figsize=self.figsize)
        im = ax.hexbin(
            self.x.ravel(), self.y.ravel(), C=self.z.ravel(), bins=None, **kw)
        plt.colorbar(im)
        plt.tight_layout()
        return plt

    def pcolormesh(self, delta=(1, 1), **kw):
        """Pcolormesh plot (x, y, z) data triplets."""
        _f, ax = plt.subplots(1, figsize=self.figsize)
        G, xc, yc = self.griddata(delta=delta)
        xx, yy = np.meshgrid(xc, yc)
        im = ax.pcolormesh(xx, yy, G, shading='flat', **kw)
        plt.colorbar(im)
        plt.tight_layout()
        return plt

    def contour(self, delta=(1, 1), **kw):
        """contour plot (x, y, z) data triplets."""
        _f, ax = plt.subplots(1, figsize=self.figsize)
        G, xc, yc = self.griddata(delta=delta)
        im = ax.contour(xc, yc, G, **kw)
        plt.clabel(im, inline=True, fontsize=6, fmt='%g')
        return plt

    def mapdata(self, use_cartopy=True, **kw):
        """Render xyz data on map.

        Parameters
        ----------
        Contour properties:
        c_ecolor : mpl.color, optional
            Set contour line colours (default '0.5')
        c_levels : int or array_like, optional
            Determines the number and positions of the contour lines / regions.
            An integer value indicates number of levels (mpl decides the
            boundaries) where are an array specifies the boundaries
            (default 15)
        c_lines : bool, optional
            Also draw contour lines on top (default: False)
        c_smooth : bool, optional
            Smooth contour plot if plt_type is set to 'contour' (default False)
        c_sm_sigma : The standard deviations of the Gaussian filter applied
            for smoothed contour (increase blurriness of the image.
            Given for each axis as a sequence, or as a single number, in which
            case it is equal for all axes. (default: 0)

        Colorbar properties:
        cb_extend : str, optional
            extend colorbar pointed ends (default: 'neither').
            Accepted values are 'neither'|'both'|'min'|'max'.
        cb_on : bool, optional
            add colorbar on map (default False)
        cb_loc : str, optional
            location of colorbar position (default 'bottom'). available
            options are 'left'|'right'|'top'|'bottom'.
        cb_pad : str, optional
            padding of colorbar from main plot in % (default '10%').
        cb_title : str, optional
            title for colorbar title (default None).
        cbt_ha : str, optional
            colorbar title horizontal alignment (default 'center').
        cbt_size : str, optional
            colorbar size (default '4%').
        cbt_pos : tuple, optional
            colorbar title (x, y) position in normal coordinate
            (default (0.5, 0.75) for `cbt_ha='centre'`
            (1, 0.75) for `cbt_ha='right'`).

        clip : bool, optional
            clip map extent within valid geographic ranges (default False).
        cmap : str, optional
            matplotlib colormap name to use (default "Spectral_r").
        delta : float or (float, float), optional
            resolution specs for binning original data in x and y direction
            (default (1, 1).
        describe_data : bool, optional
            add data statistical description (min, max, mean, ..) for
            un-gridded data
        drawcountries : bool, optional
            draw country boundaries (default False).
        drawstates : bool, optional
            draw USA state administrative boundaries (default False).
        drawstates_ind : bool, optional
            draw India state administrative boundaries (default False).
        figsize : (number, number), optional
            output figure size (default auto adjusted with a base height of 5).
        figheight : number, optional
            control output figure base height (default 5)

        fillcontinents : bool, optional
            Underlay continent with solid colour.  To mask (overlay) continent
            set mask_land=True.
        land_color : str, optional
            Fill land colour (default: grey)
        lake_color : str, optional
            Fill lake colour (default: none)

        Grid(line) properties:
        gcol : str, optional
            colours for parallel and meridian lines (default 'grey').
        gline : (number, number), optional
            grid line style that is dash pattern for meridians and parallels
            (default (None, None) for solid line). For example,
            ``gline = (1, 1)`` will draw 1 pixel on, 1 pixel off.
            Available for Basemap version only.
        gtextcolor : str, optional
            gridline text colour (default '#333333') available for Cartopy only
        gtfamily : str, optional
            gridline text family (default 'monospace').
        gtextsize : number, optional
            gridline text size (default 8).
        globe : bool, optional
            set map extent to global extent (180W:180E, 90S:90N) (default
            False). This keyword overrides limit values.
        gspacing : (number, number), optional
            spacing between grid lines (default (30, 30)).
        gzorder : number, optional
            zorder of gridlines (default 2).  Increase value for visibility
            over other masks/data.

        limit : [[float,float],[float,float]], optional
            output map extent geographic values [[lon0,lon1], [lat0,lat1]]
            (default is calculated from min max values).
        map_res : str, optional
            coastline resolution on map (default 'c' or coarse). available
            options are: 'c'oarse|'l'ow|'h'igh.
        map_buffer : float, optional
            adds extra buffer to map_limit (default None).
        pc_midpoints : bool, optional
            use grid centres to be the midpoints of pcolormesh patches
            (default True).
        plt_type : str, optional
            plot type option (default 'pcolormesh'). available options are
            'hexbin'|'contour'|'tricontour'|pcolormesh'|'scatter'.
        projection : str, optional
            output map projection (default 'cyl'indrical).

        Scatter plot/marker properties:
        s_color : str, optional
            marker colours for scatter points (default z data).
        s_marker : str, optional
            marker style for scatter points (default 's'quare).
        s_ms : int, optional
            marker size for scatter points (default 20).
        s_lw : number, optional
            linewidth of non-filled markers (e.g., '+', 'x')

        show_datapoints : overlay markers at original data points.
        show_gridpoints : overlay markers at grid points.

        Statistic:
        stat : str, optional
            statistic to use for data binning before plotting (default 'mean')
            Other options are
            mean|avg|average, std, median, count, sum|total, min, max.
            Not used when ``plt_type='scatter'``.
        norm : bool, optional
            normalize statistic values between 0 and 1 (default False).
            This also applies for plt_type='scatter'

        title : str, optional
            title string on map (default is no title).
        use_cartopy : bool, optional
            Force maps to use cartopy in place of default basemap.
        verbose : bool, optional
            verbose mode (default False).

        Data scaling:
        vmax : float, optional
            maximum value for plot scaling (default max(data)).
        vmin : float, optional
            minimum value for plot scaling (default min(data))

        Returns
        -------
        matplotlib.pyplt object

        Raises
        ------
        ImportError
            Basemap required for mapping if ``use_cartopy=False``
        NotImplementedError
            Raised when projections other than 'cyl', 'laea', 'lcc', 'merc',
            'mill', 'moll', 'robin', ot 'sinu' used with Basemap

        """
        if cartopy is None and Basemap is None:
            raise ImportError('Neither Basemap nor cartopy available.')
        if Basemap is None:
            use_cartopy = True

        stat = kw.get('stat', 'mean')
        stat = stat.lower()
        describe_data = kw.get('describe_data', False)
        plotype = kw.get('plt_type', 'pcolormesh')
        plotype = plotype.lower()
        # ms = kw.get('markersize', 20)
        cbaron = kw.get('cb_on', True)
        cbloc = kw.get('cb_loc', 'bottom')
        cbpad = kw.get('cb_pad', self.cbpad)
        cbextend = kw.get('cb_extend', 'neither')
        cbsz = kw.get('cb_size', '4%')
        cbtitle = kw.get('cb_title', stat)
        auto_cbtha = ['right', 'center'][cbextend in ('max', 'both')]
        cbtha = kw.get('cbt_ha', auto_cbtha)
        cbtpos = kw.get('cbt_pos', [(0.5, 0.75), (1, 0.75)][cbtha == 'right'])
        cmap = kw.get('cmap', 'Spectral_r')
        delta = kw.get('delta', self.delta)
        drawcountries = kw.get('drawcountries', False)
        drawstates = kw.get('drawstates', False)
        drawstates_ind = kw.get('drawstates_ind', False)
        figsize = kw.get('figsize', self.figsize)
        fight = kw.get('figheight', self.figheight)
        gcol = kw.get('gcol', 'gray')
        gline = kw.get('gline', (None, None))
        gzorder = kw.get('gzorder', 2)
        gtc = kw.get('gtextcolor', '#333333')
        gtf = kw.get('gtfamily', 'monospace')
        gts = kw.get('gtextsize', 8)
        globe = kw.get('globe', False)
        midpoints = kw.get('pc_midpoints', True)
        show_datapoints = kw.get('show_datapoints', False)
        datapoint_size = kw.get('datapoint_size', 2)
        datapoint_color = kw.get('datapoint_size', 'k')
        show_gridpoints = kw.get('show_gridpoints', False)
        fillcontinents = kw.get('fillcontinents', False)
        filloceans = kw.get('filloceans', False)
        land_color = kw.get('land_color', '#DEDEDE')
        lake_color = kw.get('lake_color', '#97B6E1')
        ocean_color = kw.get('ocean_color', '#97B6E1')
        mask_land = kw.get('mask_land', False)
        mask_ocean = kw.get('mask_ocean', False)
        zorder, lm_zorder, om_zorder = 1, 1, 1
        if mask_land:
            fillcontinents = True
            lm_zorder = 2
        if mask_ocean:
            filloceans = True
            om_zorder = 2

        limit = np.array(
            kw.get('limit', np.array(
                [[self.x.min(), self.x.max()], [self.y.min(), self.y.max()]])))

        mbuf = kw.get('map_buffer', None)
        if mbuf:
            limit += np.array([[-mbuf, mbuf], [-mbuf, mbuf]])
            for i in (0, 1):
                limit[0][i] = min(max(limit[0][i], -180), 180)
                limit[1][i] = min(max(limit[1][i], -90), 90)

        # - normalise stat
        norm = kw.get('norm', False)

        if globe is True:
            limit = self.limit
            gspacing = kw.get('gspacing', self.gspacing)
        else:
            auto_gs = (np.floor(limit[0].ptp() / 4),
                       np.floor(limit[1].ptp() / 4))
            gspacing = kw.get('gspacing', auto_gs)

        clip = kw.get('clip', False)
        mres = kw.get('map_res', 'l')
        title = kw.get('title', ' ')
        vmin = kw.get('vmin', self.z.min())
        vmax = kw.get('vmax', self.z.max())
        verb = kw.get('verbose', False)

        if plotype in ('contour', 'hexbin', 'pcolormesh'):
            if self.G is None:
                self.G, self.xc, self.yc = self.griddata(
                    stat=stat, delta=delta, norm=norm, order=True)
            # - 2D grid centres
            xxc, yyc = np.meshgrid(self.xc, self.yc)

            if stat in ('count', 'sum'):
                # vmin, vmax = self.G.min(), self.G.max()
                vmin = kw.get('vmin', self.G.min())
                vmax = kw.get('vmax', self.G.max())

        if norm:
            vmin, vmax = 0, 1
            cbtitle += ' (norm)'

        if clip is True:
            limit[0] = np.clip(limit[0], -180., 180.)
            limit[1] = np.clip(limit[0], -90., 90.)

        clon = np.mean(limit[0])
        clat = np.mean(limit[1])
        dx = np.ptp(np.array(limit[0], 'f'))
        dy = np.ptp(np.array(limit[1], 'f'))
        aspect = dx / dy
        yoff = [1, 0.5][cbloc in ('left', 'right')]
        figsize = kw.get('figsize', (fight * aspect + 0.5, fight + yoff))

        if verb:
            log.info('figsize={}'.format(figsize))

        # - draw map
        ax = kw.get('ax', None)
        if ax is None:
            fig = plt.figure(figsize=figsize)
        else:
            fig = mpl.pyplot.gcf()

        if Basemap is None or use_cartopy is True:  # ------------ use cartopy

            data_proj = cartopy.crs.PlateCarree()
            proj = kw.get('projection', data_proj)
            try:
                proj_name = proj.__class__.__name__
            except AttributeError:
                raise AttributeError(
                    '{} is not a valid ccrs proection class'.format(proj))

            font_kw = {'family': gtf, 'size': gts, 'color': gtc}
            gl_kw = {
                'linestyle': '-', 'color': gcol, 'alpha': 0.5,
                'xlocs': self._get_lonts(gspacing[0]),
                'ylocs': self._get_latts(gspacing[1]),
                'zorder': gzorder,
            }
            if mres in ('h', '10m'):
                mres = '10m'
            elif mres in ('l', '50m'):
                mres = '50m'
            else:
                mres = '110m'

            if 'ax' not in kw:
                ax = fig.add_subplot(111, projection=proj, title=title)
                # ax = plt.axes(projection=proj, title=title)

            # if globe is False:
            #     ax.set_extent(limit.ravel(), crs=data_proj)
            ax.set_extent(limit.ravel(), crs=data_proj)

            if drawcountries:
                countries = cartopy.feature.NaturalEarthFeature(
                    category='cultural', name='admin_0_countries',
                    scale=mres, linewidth=0.5,
                    facecolor='none', edgecolor='k', alpha=0.6,)
                ax.add_feature(countries, zorder=lm_zorder + 3)
                # ax.add_feature(
                #     cartopy.feature.BORDERS, linewidth=0.4, alpha=0.5)

            if drawstates:
                states_provinces = cartopy.feature.NaturalEarthFeature(
                    category='cultural', name='admin_1_states_provinces_lines',
                    scale='10m', linewidth=0.4,
                    facecolor='none', edgecolor='k', alpha=0.5)
                ax.add_feature(states_provinces)

            if drawstates_ind:
                # add Indian states
                import cartopy.io.shapereader as shpreader
                adm_1 = os.path.expandvars(
                    '$HOME/.local/share/cartopy/shapefiles/gadm/cultural/'
                    'gadm36_IND_1.shp'
                )
                adm1_shapes = list(shpreader.Reader(adm_1).geometries())
                ax.add_geometries(
                    adm1_shapes, cartopy.crs.PlateCarree(),
                    edgecolor='k',
                    linewidth=0.5,
                    linestyle='dotted',
                    facecolor='none',
                    # alpha=0.5,
                    zorder=3,
                )

            # add continent, coastline
            if fillcontinents:
                ax.add_feature(cartopy.feature.LAND, zorder=lm_zorder)
                ax.add_feature(cartopy.feature.LAKES, zorder=lm_zorder,
                               edgecolor='k', linewidth=0.5)
            if filloceans:
                ax.add_feature(cartopy.feature.OCEAN, zorder=om_zorder)

            ax.coastlines(lw=0.5, resolution=mres, zorder=lm_zorder + 3)

            if proj_name in ('PlateCarree', 'Mercator'):
                gl = ax.gridlines(draw_labels=True, **gl_kw)
                gl.xlabels_top, gl.ylabels_right = False, False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style, gl.ylabel_style = font_kw, font_kw
            else:
                ax.gridlines(**gl_kw)

        else:  # ------------------------------------------------- use Basemap

            fig.add_axes([0.05, 0.05, 0.9, 0.9], title=title)
            proj = kw.get('projection', 'cyl')
            grid_kw = {
                'family': gtf,
                'fontsize': gts,
                'color': gcol,
                'linewidth': 0.2,
                'dashes': gline,
                'zorder': gzorder,
            }
            valid_proj = (
                'cyl', 'laea', 'lcc', 'merc', 'mill', 'moll', 'robin', 'sinu',
            )

            if proj not in (valid_proj):
                raise NotImplementedError(
                    '\'{}\' not implemented yet.\n'
                    'Available projections are: {}\n'.format(proj, valid_proj))
            try:
                # for pseudo-cylindrical projection set the central coordinate
                # to 0, 0.
                # Note: Geographic subset is not respected for these
                # projections with Basemap (use Cartopy instead).
                if proj in ('moll', 'robin', 'sinu'):
                    clat, clon = 0, 0

                ax = Basemap(
                    projection=proj, lat_0=clat, lon_0=clon, resolution=mres,
                    llcrnrlon=limit[0][0], llcrnrlat=limit[1][0],
                    urcrnrlon=limit[0][1], urcrnrlat=limit[1][1])
            except ValueError as err:
                raise ValueError('{}\nTip: Use mapdata(clip=True)'.format(err))

            ax.drawmapboundary(fill_color='1')
            ax.drawmeridians(
                self._get_lonts(gspacing[0]), labels=self.xlab, **grid_kw)
            ax.drawparallels(
                self._get_latts(gspacing[1]), labels=self.ylab, **grid_kw)
            if fillcontinents:
                ax.fillcontinents(
                    color=land_color, lake_color=lake_color, zorder=lm_zorder)
            if filloceans:
                ax.drawlsmask(
                    land_color='none', ocean_color=ocean_color, lakes=True,
                    zorder=om_zorder)
            if drawcountries:
                ax.drawcountries(
                    linewidth=0.4, color='#333333', zorder=lm_zorder + 3)
            ax.drawcoastlines(linewidth=0.5, zorder=lm_zorder + 3)

        # - adjust subplot position based on color bar location
        if cbloc in ('top', 'bottom'):
            orient = 'horizontal'
            bot = [0.05, 0.15][cbloc == 'bottom']
            top = [0.85, 0.9][cbloc == 'bottom']
            cb_tick_kw = {
                'bottom': True,
                'direction': 'in',
                'labelsize': 9,
                'top': True,
            }
            cb_title_kw = {'fontsize': 10, 'ha': cbtha, 'position': cbtpos}
            fig.subplots_adjust(
                bottom=bot, top=top, left=0.075, right=0.95,
                hspace=0, wspace=0)
        elif cbloc in ('left', 'right'):
            orient = 'vertical'
            left = [0.15, 0.075][cbloc == 'right']
            right = [0.95, 0.85][cbloc == 'right']
            cbpad = kw.get('cb_pad', '5%')
            cb_tick_kw = {
                'direction': 'in',
                'labelsize': 9,
                'left': True,
                'right': True,
            }
            cb_title_kw = {
                'fontsize': 10,
                'position': [-0.25, [0.98, 0.5][cbtha == 'center']],
                'rotation': 'vertical',
                'va': ['top', 'center'][cbtha == 'center'],
            }
            fig.subplots_adjust(
                bottom=0.15, top=0.85, left=left, right=right,
                hspace=0, wspace=0)

        # - add data on map
        # -- hexbin
        if plotype in 'hexbin':
            if use_cartopy:
                xxt, yyt = transform_cartopy_coord(
                    xxc.ravel(), yyc.ravel(), data_proj, proj)
            else:
                xxt, yyt = transform_basemap_coord(
                    xxc.ravel(), yyc.ravel(), ax)

            im1 = ax.hexbin(xxt, yyt, C=self.G.data.ravel(),
                            gridsize=self.G.data.shape[::-1],
                            vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder + 1)

        elif plotype in ('contour', 'tricontour'):
            # - contour properties
            c_levels = kw.get('c_levels', 15)
            c_smooth = kw.get('c_smooth', False)
            c_sm_sigma = kw.get('c_sm_sigma', 0)
            c_lines = kw.get('c_lines', False)
            c_color = kw.get('c_color', '0.5')
            c_lw = kw.get('c_lw', 0.5)
            c_lfs = kw.get('c_lfs', 'x-small')
            norm = mpl.colors.Normalize()  # normalise cmap colours

            # -- tricontour(f): experimental
            if plotype == 'tricontour':
                if use_cartopy:
                    xt, yt = transform_cartopy_coord(
                        self.x, self.y, data_proj, proj)
                else:
                    xt, yt = transform_basemap_coord(
                        self.x.ravel(), self.y.ravel(), ax)

                im1 = plt.tricontourf(xt, yt, self.z, c_levels,
                                      norm=norm, cmap=cmap, zorder=zorder + 1)

                if c_lines:
                    im2 = plt.tricontour(
                        xt, yt, self.z, c_levels, colors=c_color,
                        linewidths=c_lw, zorder=zorder + 2)
                xxt, yyt = xt, yt

            # -- contour(f)
            else:
                if use_cartopy:
                    xxt, yyt = transform_cartopy_coord(
                        xxc, yyc, data_proj, proj)
                else:
                    xxt, yyt = transform_basemap_coord(
                        xxc.ravel(), yyc.ravel(), ax)
                    xxt = np.reshape(xxt, (len(self.yc), len(self.xc)))
                    yyt = np.reshape(yyt, (len(self.yc), len(self.xc)))

                if c_smooth:
                    im1 = ax.contourf(
                        xxt, yyt,
                        gaussian_filter(self.G, sigma=c_sm_sigma, order=0),
                        c_levels, norm=norm, cmap=cmap, zorder=zorder + 1)
                else:
                    im1 = ax.contourf(xxt, yyt, self.G, c_levels,
                                      norm=norm, cmap=cmap, zorder=zorder + 1)

                # overlay contour lines and draw labels
                if c_lines:
                    im2 = plt.contour(
                        im1, levels=im1.levels[::2], colors=c_color,
                        linewidths=c_lw, zorder=zorder + 2)
                    _family = mpl.rcParams['font.family']
                    mpl.rcParams['font.family'] = 'monospace'
                    plt.clabel(im2, inline=True, fontsize=c_lfs, fmt='%g')
                    mpl.rcParams['font.family'] = _family

        # -- scatter
        elif plotype in 'scatter':
            if stat:
                log.warning('stat=%s is not computed for for plt_type=%s',
                            stat, plotype)

            cbtitle = kw.get('cb_title', ' ')
            s_ms = kw.get('s_ms', 20)
            s_marker = kw.get('s_marker', 's')
            s_lw = kw.get('s_lw', 1)
            s_color = kw.get('s_color', self.z)
            # - switch off colour bar for single colour markers
            if len(s_color) < 2:
                cbaron = False

            if use_cartopy:
                xt, yt = transform_cartopy_coord(
                    self.x, self.y, data_proj, proj)
            else:
                xt, yt = transform_basemap_coord(self.x, self.y, ax)
            im1 = ax.scatter(
                xt, yt, c=s_color, marker=s_marker, s=s_ms, lw=s_lw,
                vmin=vmin, vmax=vmax, cmap=cmap, zorder=zorder + 1
            )
            xxt, yyt = xt, yt

        # -- pcolormesh
        elif plotype in 'pcolormesh':
            if midpoints:
                xc, yc = self.shift_midpoints()
                xxc, yyc = np.meshgrid(xc, yc)
            else:
                xc, yc = self.xc, self.yc

            xxx, yyy = np.meshgrid(self.xc, self.yc)

            if use_cartopy:
                xxt, yyt = transform_cartopy_coord(xxc, yyc, data_proj, proj)
                xct, yct = transform_cartopy_coord(xxx, yyy, data_proj, proj)
            else:
                xxt, yyt = transform_basemap_coord(
                    xxc.ravel(), yyc.ravel(), ax)
                xxt = np.reshape(xxt, (len(yc), len(xc)))
                yyt = np.reshape(yyt, (len(yc), len(xc)))
                xct, yct = transform_basemap_coord(
                    xxx.ravel(), yyy.ravel(), ax)

            im1 = ax.pcolormesh(
                xxt, yyt, self.G, vmin=vmin, vmax=vmax, rasterized=True,
                shading='flat', cmap=cmap, zorder=zorder + 1)

        # -- invalid plt_type
        else:
            raise NotImplementedError(plotype)

        # - overlay data points
        if show_datapoints:
            if use_cartopy:
                xt, yt = transform_cartopy_coord(
                    self.x, self.y, data_proj, proj)
            else:
                xt, yt = transform_basemap_coord(self.x, self.y, ax)
            ax.scatter(xt, yt, marker='+', s=datapoint_size, linewidth=0.5,
                       c=datapoint_color, zorder=6, alpha=0.5)

        # - overlay grid centres
        if show_gridpoints:
            sckw = dict(marker='+', linewidth=0.5, c='k', zorder=6, alpha=0.5)
            if plotype == 'pcolormesh':
                ax.scatter(xct, yct, **sckw)
            else:
                ax.scatter(xxt, yyt, **sckw)

        if describe_data:
            s = describe(self.z, nan_policy='omit')
            if plotype == 'scatter':
                stat_str = (
                    'min:{:.4g}|max:{:.4g}|avg:{:.4g}\n'
                    'var:{:.4g}|skw:{:.4g}|nob:{:g}').format(
                        s.minmax[0], s.minmax[1], s.mean,
                        s.variance, s.skewness, s.nobs)
            else:
                sg = describe(self.G.flatten(), nan_policy='omit')
                nobs = self.G.count()

                stat_str = (
                    ' grid [original]\n'
                    'min {:.4g} [{:.4g}]\nmax {:.4g} [{:.4g}]\n'
                    'avg {:.4g} [{:.4g}]\n'
                    'var {:.4g} [{:.4g}]\nskw {:.4g} [{:.4g}]\n'
                    '  N {:g} [{:g}]').format(
                        sg.minmax[0], s.minmax[0], sg.minmax[1], s.minmax[1],
                        sg.mean, s.mean, sg.variance, s.variance,
                        sg.skewness, s.skewness, nobs, s.nobs)

            style = dict(
                bbox=dict(
                    boxstyle='square,pad=0', ec='none', fc=(1, 1, 1, 0.9)
                ),
                color='b',
                family='monospace',
                fontsize='small',
                verticalalignment='top',
            )

            if use_cartopy:
                ax.text(0.01, 0.99, stat_str, transform=ax.transAxes,
                        zorder=11, **style)
            else:
                plt.annotate(stat_str, xy=(0.01, 0.991),
                             xycoords='axes fraction', zorder=11, **style)

        if cbaron:
            if Basemap is None or use_cartopy is True:
                # - Cartopy
                cb_ax = fig.add_axes([0, 0, .1, .1])  # (x0, y0, wide, high)
                def_cb_pad = ['7%', '5%'][cbloc in ('left', 'right')]
                cbpad = kw.get('cb_pad', def_cb_pad)
                fig.canvas.mpl_connect(
                    'resize_event', resize_cbar_function(
                        ax, cb_ax, location=cbloc, pad=cbpad, size=cbsz))

                cb = plt.colorbar(
                    im1, orientation=orient, cax=cb_ax, format='%g',
                    extend=cbextend)
            else:
                # - Basemap
                cb = ax.colorbar(
                    im1, location=cbloc, size=cbsz, pad=cbpad, format='%g',
                    extend=cbextend)

            cb.ax.tick_params(**cb_tick_kw)
            cb.ax.set_title('{}'.format(cbtitle), **cb_title_kw)
            # from matplotlib import ticker
            # tick_locator = ticker.MaxNLocator(nbins=5)
            # cb.locator = tick_locator
            # cb.update_ticks()

            fig.canvas.draw_idle()
        return plt

    def __exit__(self, exc_type, exc_value, traceback):
        return isinstance(exc_value, TypeError)


def aot2pm25(aot, H, S):
    """Convert AOT to PM 2.5.

    AOD and PM2.5 relationship (Hoff and Christopher, 2009):
        `AOT = PM2.5*H*S`

    Parameters
    ----------
    aot : number
        Aerosol optical thickness
    H : number
        Height of the well-mixed planetary boundary layer (PBL) and
    S : Specific extinction efficiency of the aerosol at the ambient
        relative humidity.

    Returns
    -------
    number
        PM2.5 corresponding to AOT loading.
    """

    return aot / (H * S)


def autoswitch_backend(verbose=False):
    """Automatically switch backend to Agg if using a non-interactive or
    non-DISPLAY system.

    """
    if platform.startswith('linux'):
        if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
            plt.switch_backend('Agg')
    if verbose:
        log.info('Using backend: ' + mpl.get_backend())


def csv2nc(csvfile, ncfile=None):
    """Convert csv file to netcdf."""
    import pandas as pd
    from netCDF4 import Dataset

    df = pd.read_csv(csvfile)
    # print(df.dtypes)
    if ncfile is None:
        ncfile = csvfile + '.nc'

    data = df.to_dict('series')
    with Dataset(ncfile, 'w') as f:
        params = f.createGroup('data')
        for k, v in list(data.items()):
            print(k)
            setattr(params, k, v)


def dialog_pickfiles(initial_dir=None, extension=None):
    """Dialog picker for one or more files.

    Parameters
    ----------
    initial_dir : str
        initial directory to select files (default ``os.getcwd()``)
    extension : str or list of str, optional
        list of file extensions (default '*.*')

    Returns
    -------
    str or list of str
        selected filename[s]

    """
    # remove empty application window that pops up behind the file dialogs
    Tk().withdraw()

    # parse initial_dir
    initial_dir = [initial_dir, os.getcwd()][initial_dir is None]
    extension = [extension, '*.*'][extension is None]
    # convert extension to list if required
    if isinstance(extension, list) is False:
        extension = extension.split()

    file_type = list(zip(['file type'] * len(extension), extension))
    files = askopenfilenames(
        initialdir=initial_dir, filetypes=file_type, title='Pick file[s]')

    if files:
        out = [f for f in files] if len(files) > 1 else files[0]
        return out
    else:
        log.warning('No files picked.')


def dict2csv(dictionary, filename, orient='columns', index=False):
    """Write csv file from a dictionary.

    Parameters
    ----------
    dictionary : dict
        Input dictionary
    filename : str
        Output csv filename
    orient : str, optional
        {'columns', 'index'}, default 'columns'
        The "orientation" of the data. If the keys of the passed dict should
        be the columns of the resulting DataFrame, pass 'columns' (default).
        Otherwise if the keys should be rows, pass 'index'.
    index : bool, optional
        Write row names (index)

    """
    import pandas as pd
    # convert dictionary to pandas data-frame
    pd.DataFrame.from_dict(
        dictionary, orient=orient).to_csv(
            filename, index=index)


def dict2nc(dictionary, ncfile, zlib=True, attr=None):
    """Convert a dictionary to simple netCDF (nf4 with compression enabled).

    Warning: it assumes all datasets in the dictionary are 1D floating point
    arrays and of equal length. Do not use this function for complex
    multi-dimensional arrays.

    TODO: add dimensions, variables, local attributes optional parameters

    Parameters
    ----------
    dictionary : dict
        Input dictionary of 1D data arrays. all variables are of equal length
    ncfile : str
        output netCDF filename
    zlib : bool, optional
        turn on or off data compression (default True)
    attr : dict, optional
        global attributes. Standard CF-1.6 global attributes are
        'Conventions', 'title', 'institution', 'source', 'references',
        'history', 'comment' (default is None which adds the file creation
        history only)

    Returns
    -------
    None

    """
    with Dataset(ncfile, 'w', format='NETCDF4') as f:
        # create record number dimensions
        f.createDimension('N', len(list(dictionary.values())[0]))

        # create variables from keys and fill data
        for key, value in list(dictionary.items()):
            var = f.createVariable(key, 'f4', ('N',), zlib=zlib)
            var[:] = value

        # add global attributes to file
        f.history = "[{}] Created netCDF4 zlib={} dataset.".format(
            datetime.today().strftime('%F %T'), zlib)
        if attr is not None:
            for key, value in list(attr.items()):
                f.temp = value
                f.renameAttribute('temp', key)


def doy(year=None, month=None, day=None):
    """Calculate serial day of year (1-366).

    Parameters
    ----------
    year : int, optional
        Year
    month : int, optional
        Month
    day : Int, optional
        Day in calendar month

    Returns
    -------
    int
        Serial day of year.
    """
    if all([year, month, day]):
        dt = date(year, month, day)
        return dt.timetuple().tm_yday
    else:
        return date.today().timetuple().tm_yday


def doy2ymd(year, dayofyear):
    """Get date object from year and seq day of year."""
    return datetime(year, 1, 1) + timedelta(dayofyear - 1)


def dt2et(date_time, start_date, date_format='%Y-%m-%d %H:%M:%S'):
    """Get list of datetime objects for elapsed seconds from a given start
    date.

    Note: this function ignores microseconds.

    Parameters
    ----------
    date_time : list of str
        input datetime
    start_date : str
        input reference datetime
    date_format : str, optional
        format specification for date_time and start_date (default
        '%Y-%m-%d %H:%M:%S')

    Examples
    --------
    >>> dt2et("2010-01-03 00:00:00", "2010-01-01 00:00:00")
    [172800.0]
    >>> dt2et(["2010-01-03 00:00:00", "2010-01-03 10:00:00"],
    ...       "2010-01-01 00:00:00")
    [172800.0, 208800.0]

    Returns
    -------
    array_like
        Elapsed time from a reference time
    """
    dt = [date_time] if np.isscalar(date_time) else date_time
    return [(datetime.strptime(d, date_format) -
             datetime.strptime(start_date, date_format)).total_seconds()
            for d in dt]


def dt2yjd(in_date, in_fmt='%Y%m%d', sep=''):
    """Convert Date to Year and Dayofyear.

    Parameters
    ----------
    in_date : str or list of str
    in_fmt : str, optional
        format of input date (default '%Y%m%d')
    sep : str, optional
        separator between output Year and Day of year (default '')

    Returns
    -------
    str or list of str
        Year and Doy combined

    Examples
    --------
    >>> dt2yjd(['20170301', '20160919'])
    ['2017060', '2016263']
    >>> dt2yjd(['20170301', '20160919'], sep='-')
    ['2017-060', '2016-263']
    >>> dt2yjd(['2017-03-01', '2016-09-19'], in_fmt='%Y-%m-%d', sep='-')
    ['2017-060', '2016-263']

    """
    return format_date(in_date, in_fmt, '%Y' + sep + '%j')


def et2dt(elapsed_seconds, start_date, date_format='%Y-%m-%d'):
    """Get list of datetime objects for elapsed seconds from a given start
    date.

    Parameters
    ----------
    elapsed_seconds : float64 or sequence of floats64
        input elapsed seconds
    start_date : str
        reference date from which elapsed time is measured
    date_format : str, optional
        format of start_date (default '%Y-%m-%d')
    """
    es = [elapsed_seconds] if np.isscalar(elapsed_seconds) else elapsed_seconds
    return [datetime.strptime(start_date, date_format) +
            timedelta(seconds=td) for td in es]


def ext2vis(ext):
    """Visibility V, as a function of the extinction coefficient ext (using
    Koschmieder's equation)

        Vk = 3.912/ext

    ext: average extinction coefficient of the aerosol over the visual range

        Vis = (1/ext) * ln(abs(C0/eps))
        C0 = inherent contracst of the object ~0.02
        eps = liminal/inherent theoretical contrast of 'black object' -1

    This formula is widely used, especially for its simplicity, and it's
    main assumptions for its derivation are:

    a) The extinction coefficient is constant along the  path of sight (this
    also means a horizontal path of sight and neglecting the curvature of the
    earth).

    b) The amount of light scattered by a volume element is proportional to
    its volume and the extinction coefficient, and also is constant along the
    path of sight.

    c) The target is absolutely black, with the horizon as background.

    d) The eye has a contrast threshold of 0.02.
    """
    return 3.912 / ext


def format_date(in_date=None, in_format='%Y%m%d', out_format='%Y/%j'):
    """Convert date string to desired format.

    Parameters
    ----------
    in_date : str or list of str, optional
        date [and time] in_format (default None)
    in_format : str, optional
        input date format (default %Y%m%d)
    out_format : str, optional
        desired output format (default %Y/%j)

    Returns
    -------
    str
        Date in desired format.

    """
    if in_date is None:
        in_date = datetime.utcnow().strftime(in_format)
    in_date = [in_date] if np.isscalar(in_date) else in_date

    return [datetime.strptime(d, in_format).strftime(out_format) for
            d in in_date]


def frange(start=0, stop=None, step=1.0):
    """List of floating point values.

    Parameters
    ----------
    start : number
        Start of interval. The interval includes this value. The default
        start value is 0.
    stop : number, optional
        End of interval. The interval does not include this value, except in
        some cases where step is not an integer and floating point round-off
        affects the length of out.
    step : number, optional
        Spacing between values. For any output out, this is the distance
        between two adjacent values, out[i+1] - out[i]. The default step size
        is 1. If step is specified as a position argument, start must also be
        given.

    Yields
    ------
    generator
        Array of evenly spaced values.

        Return a generator containing an arithmetic progression of numbers.
        frange(i, j) returns [i, i+1, i+2, ..., j-1]; start (!) defaults to 0.
        When step is given, it specifies the increment (or decrement).
        For example, list(frange(4)) returns [0., 1.0, 2.0, 3.0].
        The end point is omitted! These are exactly the valid indices for a
        list of 4 elements.

    """
    if stop is None:
        stop = start
        start = 0
    i = start
    if start < stop and step > 0:
        while i < stop:
            yield i
            i += step
    elif start > stop and step < 0:
        while i > stop:
            yield i
            i += step


def freq2wl(freq):
    """Convert frequency (cm-1) to wavelength (microns)."""
    return 1e4 / freq


def getnn(data1, data2, r, k=5, p=2.0, eps=0.0, n_proc=1):
    """Search nearest neighbours between two coordinate catalogues.

    See https://docs.scipy.org/doc/scipy/reference/generated/
    scipy.spatial.cKDTree.query.html

    Parameters
    ----------
    data1, data2 : array_like, shape (n, m)
        Arrays of data points.  These are the n data points of dimension m
        to be indexed and queried (in (n, 2)).

    r : nonnegative float
        Return only neighbours within this distance. This is used to prune
        tree searches, so if you are doing a series of nearest-neighbour
        queries, it may help to supply the distance to the nearest neighbour
        of the most recent point.

    k : int of list of ints
        The list of k-th nearest neighbours to return. If k is an integer
        it is treated as a list of [1, ... k] (range(1, k+1)). Note that
        the counting starts from 1.

    p : float, 1<=p<=infinity
        Which Minkowski p-norm to use. 1 is the sum-of-absolute-values
        "Manhattan" distance, 2 is the usual Euclidean distance, infinity is
        the maximum-coordinate-difference distance.

    eps : non-negative float
        Return approximate nearest neighbours; the k-th returned value is
        guaranteed to be no further than (1+eps) times the distance to the
        real k-th nearest neighbour.

    n_proc : int, optional
        Number of jobs to schedule for parallel processing. If -1 is given
        all processors are used. Default: 1.  Requires scipy version >= 0.16.0

    Returns
    -------
    indices : ndarray of ints
        The locations of the neighbours in self.data. If x has shape
        tuple+(self.m,), then i has shape tuple+(k,). When k == 1, the last
        dimension of the output is squeezed. Missing neighbours are indicated
        with self.n.

    distances : array of floats
        The distances to the nearest neighbours. If x has shape
        tuple+(self.m,), then d has shape tuple+(k,). When k == 1, the last
        dimension of the output is squeezed. Missing neighbours are indicated
        with infinite distances.

    Example:
        >>> import numpy
        >>> from ypylib.utils import getnn
        >>> Lon1 = numpy.random.random(2000)
        >>> Lat1 = numpy.random.random(2000)
        >>> Lon2 = numpy.random.random(20)
        >>> Lat2 = numpy.random.random(20)
        >>> data1 = numpy.array(list(zip(Lon1, Lat1)))
        >>> data2 = numpy.array(list(zip(Lon2, Lat2)))
        >>> i, d = getnn(data1, data2, 0.1, k=3)
    """
    tree = cKDTree(data1)
    try:
        distances, indices = tree.query(
            data2, k=k, eps=eps, p=p, distance_upper_bound=r, n_jobs=n_proc)
    except TypeError as err:
        # for older version of scipy n_jobs is not implemented
        msg = ("which was introduced in scipy-0.16.0; hence getnn() "
               "will ignore this option.")
        log.warn('%s - %s', err, msg)
        distances, indices = tree.query(
            data2, k=k, eps=eps, p=p, distance_upper_bound=r)

    return indices, distances


def get_public_methods(class_name):
    """Return public method names in a class.

    Parameters
    ----------
    class_name : class
        input class

    Returns
    -------
    tuple
        available public method names
    """
    return list(zip(*inspect.getmembers(
        class_name, predicate=inspect.ismethod)))[0][1:]


def ll_parse(ll_string):
    """Parse a coordinate string to equivalent number.

    Parameters
    ----------
    ll_string : str
        Coordinate string annotated with E|W|N|S, e.g., 20N = 20.0, 10W = -10.0

    Returns
    -------
    real
        Equivalent real
    """
    mult = [-1, 1][ll_string[-1].upper() in ['N', 'E']]
    return mult * float(ll_string[:-1])


def load_nc(filename, variables=None, verb=False, gattr=False, order=False):
    """Load variables from a netCDF file to a dictionary.

    Parameters
    ----------
    filename : str
        input netCDF filename
    variables : sequence of str, optional
        Sequence of variable names to read from the file. default is to read
        all variables.
    verb : bool, optional
        Verbose mode - print attributes of loaded variables.
    gattr : bool, optional
        Print global attributes of the netCDF file.
    order : bool, optional
        Retain variable sequence as requested in the output dictionary
        (ordered).

    Returns
    -------
    dict or OrderedDict
        (Ordered) Dictionary of requested or all variables from the file.

    """
    out = collections.OrderedDict() if order else {}

    basestring = str  # python 2/3 compatible

    with Dataset(filename, 'r') as nc:

        # Dimensions
        if verb:
            dims = nc.dimensions
            print('dimensions:')
            for j in dims:
                print('   {} = {}'.format(dims[j].name, dims[j].size))
            print('loaded variables:')

        # Parse variable names
        if variables is None:
            variables = list(nc.variables.keys())
        elif isinstance(variables, basestring):
            variables = (variables,)

        # Iterate through requested variables and update out dictionary
        for item in variables:
            try:
                out.update({item: nc.variables[item][:]})

                # be verbose?
                if verb:
                    # Variable name, dims, type, data type
                    print('   {} {} {} {}'.format(
                        item,
                        nc.variables[item].dimensions,
                        type(nc.variables[item][:]),
                        nc.variables[item].dtype))

                    # Variable Attributes
                    for key in sorted(nc.variables[item].ncattrs()):
                        print('      {}: {}'.format(
                            key, getattr(nc.variables[item], key)))

            except KeyError:
                log.error('%s: No such variable in %s', item, filename)

        # Global attributes
        if gattr:
            print('\n// global attributes:')
            if len(nc.ncattrs()) == 0:
                print(None)
            else:
                for attr in sorted(nc.ncattrs()):
                    print('      :{} = {}'.format(attr, getattr(nc, attr)))
    # Return data dict
    return out


def locate(array, value, epsilon=1.0):
    """Locate the index of a value in an array."""
    inarray = np.asanyarray(array, dtype=array.dtype)
    loc = (np.abs(value - inarray)).argmin()
    if np.abs(value - inarray[loc]) < epsilon:
        return loc
    else:
        return -1


def merge_dicts(dict_seq):
    """Merge dictionaries while keeping the original keys and values of each
    dictionary.

    Parameters
    ----------
    dict_seq : sequence of dicts
        Input dictionary sequence (list or tuple) to merge.  Each dict in
        dict_seq are expected to be simple dict class. Mixing `dicts` and
        `defaultdicts` in dict_seq may produce inconsistent result.

    Returns
    -------
    collections.defaultdict
        Merged dictionary will all key values. All values are lists.
    """
    # from collections import defaultdict

    out = collections.defaultdict(list)
    for i in range(len(dict_seq)):
        if dict_seq[i].__class__ is not dict:
            raise TypeError('dict_seq {} is not a dict'.format(i + 1))
        for k, v in dict_seq[i].items():
            out[k].append(v)
    return out


def mkdirp(path):
    """Make parent directories as needed, no error if existing.

    Parameters
    ----------
    path : str
        name of directory (tree) to create

    Raises
    ------
    err
        Description
    """
    try:
        os.makedirs(path)
    except OSError as err:
        if err.errno != errno.EEXIST:
            # ignore if path already exists, else raise error.
            raise err


def modify_cmap(index_list=None, rgba_list=None, N=None, in_cmap=None,
                name='Custom Cmap'):
    """Modify rgba values of an existing colormap.

    Parameters
    ----------
    index_list : list of ints
        list of indices in colour table for which needs customisation
    rgba_list : list of tuples
        list of input (r,g,b,a) tuples corresponding to index. Note that
        starting with matplotlib 2, provide rgba values in (0-1) range.
    N : int, optional
        the number of entries in the cmap. The default is None, in which case
        there is one colormap entry for each element in the list of colours
    in_cmap : mpl colormap, optional
        Use this colormap as input instead of default cmap.
    name : str, optional
        string to identify the customised colormap (default 'Custom Cmap')

    Returns
    -------
    matplotlib.colors.ListedColormap object

    Examples
    --------
        >>> my_cmap1 = modify_cmap([0], [(255, 255, 255, 0.5)])

        >>> mpl.rc('image', cmap='gray')  # change default cmap to gray
        >>> my_cmap2 = modify_cmap([0], [(0, 0, 0, 0.5)])
    """
    if in_cmap is None:
        ccm = mpl.cm.get_cmap()  # get default colormap
    else:
        ccm = mpl.cm.get_cmap(in_cmap)

    ccl = [ccm(i) for i in range(ccm.N)]  # build custom color list
    if index_list and rgba_list:
        for i, _ in enumerate(index_list):
            ccl[index_list[i]] = rgba_list[i]  # modify colors

    return mpl.colors.ListedColormap(ccl, name, N)


def nearest(array, value, location=False):
    """Returns the nearest value in an array.

    Parameters
    ----------
    array : array_like
        search array from which the nearest value to be found
    value : scalar
        to search
    location : bool, optional
        return locations of nearest values (default False)

    Returns
    -------
    index, value : tuple
        index or position of the value in the array and the actual value in the
        array if location set to True, otherwise returns the closest value.
    """
    inarray = np.asarray(array)
    loc = (np.abs(inarray - value)).argmin()
    if location is True:
        return (loc, inarray[loc])
    else:
        return inarray[loc]


def pct2num(string):
    """Convert percent string to number."""
    return float(string.strip('%')) / 100


def purge(dsrc, ext):
    """Delete all files with matching extension in a directory.

    Parameters
    ----------
    dsrc : str
        source directory
    ext : str
        matching file extension
    """
    exts = '*' + str(ext).strip()
    for f in glob.glob(os.path.join(dsrc, exts)):
        os.remove(f)


def read_nc3(ncfile, var=[], version=1):
    """Read variables from a NetCDF-3 (classic) file.

    Parameters
    ----------
    ncfile : str
        netCDF-3 (classic) filename
    var : list of str, optional
        List of variable names to read
    version : {1, 2}, optional
        version of netcdf to read, where 1 means Classic format and
        2 means 64-bit offset format. Default is 1.

    Returns
    -------
    dict
        Dictionary with requested vars (as keys)

    """
    from scipy.io import netcdf
    with netcdf.netcdf_file(ncfile, 'r', mmap=False, version=version) as nc:
        if var == []:
            var = nc.variables.keys()

        return dict((v, nc.variables[v].data) for v in var)


def resize_cbar_function(ax, cax, location='bottom', pad='10%', size='4%'):
    """Returns a function to automatically resize the colorbar for cartopy
    plots.

    Parameters
    ----------
    ax : axis
    cax : colorbar axis
    location : str, optional
        location of colorbar ('top' | 'bottom' | 'left' |'right')
    pad : str, optional
        colorbar padding from main plot (default '10%')
    size : str, optional
        relative size of colorbar (default '4%')

    Examples
    --------
        >>> import cartopy.crs as ccrs
        >>> import matplotlib.pyplot as plt

        >>> fig, ax = plt.subplots(
        ...     figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})
        >>> cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])

        [... your code generating a scalar mappable ...]

        >>> resize_colorbar = resize_cbar_function(ax, cbar_ax)
        >>> _ = fig.canvas.mpl_connect('resize_event', resize_colorbar)

    Credits
    -------
    Solution by pelson at http://stackoverflow.com/a/30077745/512111
    """
    def resize_colorbar(event):
        # Tell matplotlib to re-draw everything, so that we can get
        # the correct location from get_position.
        plt.draw()
        pos = ax.get_position()
        if location == 'bottom':
            cax_pos = [
                pos.xmin, pos.ymin - (pct2num(pad) + pct2num(size) / 1.5),
                pos.width, pct2num(size) / 1.5]
        if location == 'top':
            cax_pos = [
                pos.xmin, pos.ymax + (pct2num(pad)),
                pos.width, pct2num(size) / 1.5]
        elif location == 'left':
            cax_pos = [
                pos.xmin - pct2num(pad), pos.ymin,
                pct2num(size) / 1.5, pos.height]
        elif location == 'right':
            cax_pos = [
                pos.xmin + pos.width + pct2num(pad), pos.ymin,
                pct2num(size) / 1.5, pos.height]

        cax.set_position(cax_pos)

    resize_colorbar(None)
    plt.sca(ax)
    return resize_colorbar


def reverse_cmap(cmap, N=256):
    """Reverse a matplotlib colormap.

    Parameters
    ----------
    cmap : LinearSegmentedColormap
        original matplotlib colormap to reverse
    N : int, optional
        number of rgb quantization levels (default 256)

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap

    """
    rev_cm_name = '_'.join([cmap.name, 'r'])
    try:
        rev_cm = mpl.colors.LinearSegmentedColormap(
            rev_cm_name, mpl_cm.revcmap(cmap._segmentdata), N=N)
    except AttributeError:
        seg_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            'temp', cmap.colors)
        rev_cm = mpl.colors.LinearSegmentedColormap(
            rev_cm_name, mpl_cm.revcmap(seg_cmap._segmentdata), N=N)

    return rev_cm


def save_notebook_figures(url):
    """
    Extraction image files from a Jupyter notebook and save them to disk.

    See https://nbconvert.readthedocs.io/en/latest/nbconvert_library.html

    Parameters
    ----------
    url : str
        Notebook URL. Example
        'https://path/to/jupyter.notebook'
        'file:///path/to/jupyter.notebook'
    """

    from urllib.request import urlopen
    import nbformat
    from nbconvert import RSTExporter

    response = urlopen(url).read().decode()
    notebook = nbformat.reads(response, as_version=4)

    # Extract the base64 encoded figures as files
    rst_exporter = RSTExporter()
    (body, resources) = rst_exporter.from_notebook_node(notebook)

    # write image hex-dump to indivisual file
    for imgfile in sorted(resources['outputs'].keys()):
        with open(imgfile, "wb") as f:
            print("writing {}...".format(imgfile))
            f.write(resources['outputs'][imgfile])


def strmatch(patterns, filename, mutually_inclusive=False):
    """Saerch lines with smatching strings in a text file.

    Parameters
    ----------
    patterns : str or list of str
        List of strings to search in a line. Note that patterns being
        searched are not mutually inclusive, i.e. find lines where any of
        the pattenrs are matched.
    filename : str
        Filename f the text file to search from
    mutually_inclusive : bool, optional
        When patterns is a list of strings, setting mutually_inclusive to True
        means user wish to find lines in file where all string patterns are
        matched.  This keyword have no effect of patterns contains only one
        string.

    Yields
    ------
    generator
        Count, line-number in file and matched line

    Examples
    --------
    >>> pattern = ['ftp.metoffice.gov.uk', 'mo_goes']
    >>> filename = os.path.expanduser('~/.netrc')
    >>> for cnt, row, line in strmatch(pattern, filename):
    ...     # commening out sensintive info from testing
    ...     # print(cnt, row, line, end='')
    ...     pass

    """
    cnt = 0
    row = 0
    set1 = set(patterns)

    with open(os.path.expanduser(filename), 'r') as f:
        for line in f:
            row += 1
            set2 = set(word for word in line.split())

            if mutually_inclusive:
                if len(set2.intersection(set1)) == len(set1):
                    cnt += 1
                    yield cnt, row, line
            else:
                if len(set2.intersection(set1)) > 0:
                    cnt += 1
                    yield cnt, row, line


def symlinkf(fsrc, fdst):
    """Create a symbolic link pointing to fsrc named fdst.

    This is a os.symlink wrapper which removes existing destination files

    Availability: Unix.

    Paremeters
    ----------
    fsrc : str
        source file
    fdst : str
        destination symbolic link
    """
    try:
        os.symlink(fsrc, fdst)
    except OSError as err:
        if err.errno == errno.EEXIST:
            os.remove(fdst)
            os.symlink(fsrc, fdst)


def seqdate(start_date, end_date, in_fmt='%Y-%m-%d', out_fmt='%Y%m%d'):
    """Generate sequence of dates between two dates.

    Parameters
    ----------
    start_date : str
        Start date
    end_date : str
        End date
    in_fmt : str, optional
        Format specification for start and end dates. default %Y-%m-%d
    out_fmt : str, optional
        Format specification for output sate sequences. default %Y%m%d

    Returns
    -------
    list
        List of date strings
    """
    d1 = datetime.strptime(start_date, in_fmt)
    d2 = datetime.strptime(end_date, in_fmt)
    delta = d2 - d1
    return [(d1 + timedelta(days=i)).strftime(out_fmt) for
            i in range(delta.days + 1)]


def transform_basemap_coord(lon, lat, to_proj):
    """Transform latitude, longitude values to specific basemap projection.

    Parameters
    ----------
    lon : array_like
        Longitude array in degrees.
    lat : array_like
        Latitude array in degrees
    to_proj : class
        mpl_toolkits.basemap.Basemap (used to define a map axis), to which
        the input coordinates will be converted.

    Returns
    -------
    tuple
        Transformed longitude, latitude arrays

    Example
    -------
    >>> ax = Basemap(projection='robin', lon_0=0)
    >>> longitudes = [10, 20, 30]
    >>> latitudes = [13, 0, -10]
    >>> tax = transform_basemap_coord(longitudes, latitudes, ax)
    >>> tax[0]
    array([17923407.08812449, 18874217.96108178, 19804905.68141247])
    >>> tax[1]
    array([10004318.19545184,  8615499.70605575,  7547177.81737116])

    """
    tran = np.array([to_proj.projtran(x, y) for x, y in zip(lon, lat)])
    return tran[..., 0], tran[..., 1]


def transform_cartopy_coord(x, y, from_proj, to_proj, truncate_xy=None):
    """Convert between projections using Cartopy tranform_points.

    Parameters
    ----------
    x : TYPE
        Description
    y : TYPE
        Description
    from_proj : TYPE
        Description
    to_proj : TYPE
        Description

    Returns
    -------
    TYPE
        Description
    """
    if truncate_xy is not None:
        x[abs(x) > truncate_xy[0]] = truncate_xy[0]
        y[abs(y) > truncate_xy[1]] = truncate_xy[1]

    tran = to_proj.transform_points(from_proj, x, y)
    return tran[..., 0], tran[..., 1]


def v_locate(array, value):
    """v_locate function finds the intervals within a given monotonic vector
    that brackets a given set of one value. Much faster than locate when the
    array is sorted (or monotonically increasing).
    """
    idx = np.searchsorted(array, value, side="left")
    if math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx]):
        return idx - 1  # array[idx - 1]
    else:
        return idx  # array[idx]


def wl2freq(wl):
    """Convert wavelength (microns) to frequency (cm-1)."""
    return 1e4 / wl


def write_csv(filename, data=None, header=None, append=False):
    mode = 'ab' if append is True else 'wb'
    with open(filename, mode) as fp:
        writer = csv.writer(fp)
        # write header row
        if append is False:
            try:
                writer.writerow(header)
            except UserWarning:
                log.warning('No header.')
        # write row row(s)
        try:
            writer.writerow(data)
        except UserWarning:
            log.warning('No data.')


def yjd2dt(yjd, fmt='%Y%j'):
    """Converts date string in YYYYjjj to datetime object.

    Parameters
    ----------
    yjd : str
        Year and day-of-year in yyyyddd form
    in_fmt: str, optional
        exact format of yjd (default '%Y%j)

    Returns
    -------
    list of datetime.date objects

    """
    yjd = [yjd] if np.isscalar(yjd) else yjd
    return [datetime.strptime(i, fmt).date() for i in yjd]


def testBit(int_type, offset):
    """testBit() returns a non-zero result, 2**offset, if the bit at 'offset'
    is one.
    """
    mask = 1 << offset
    return(int_type & mask)


def setBit(int_type, offset):
    """setBit() returns an integer with the bit at 'offset' set to 1."""
    mask = 1 << offset
    return(int_type | mask)


def clearBit(int_type, offset):
    """clearBit() returns an integer with the bit at 'offset' cleared."""
    mask = ~(1 << offset)
    return(int_type & mask)


def toggleBit(int_type, offset):
    """toggleBit() returns an integer with the bit at 'offset' inverted,
    0 -> 1 and 1 -> 0.
    """
    mask = 1 << offset
    return(int_type ^ mask)


if __name__ == "__main__":
    pass
    # import doctest
    # doctest.testmod()
    # print("done.")
