#!/usr/bin/env python2.7
"""
:author: yaswant.pradhan
:copyright: Crown copyright. Met Office.
"""
import os
import csv
import glob
import errno
import math
import socket
import numpy as np
import matplotlib as mpl
from sys import platform
from netCDF4 import Dataset
from collections import Sequence
from datetime import datetime, timedelta
from datetime import date
from scipy.spatial import cKDTree
from scipy.stats import binned_statistic_2d
from matplotlib import pyplot as plt
from Tkinter import Tk
from tkFileDialog import askopenfilenames
try:
    from mpl_toolkits.basemap import Basemap
except ImportError:
    Basemap = None
try:
    import cartopy
    from cartopy.mpl.gridliner import \
        LONGITUDE_FORMATTER, \
        LATITUDE_FORMATTER
except ImportError:
    cartopy = None
if platform in ('linux', 'linux2'):
    if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
        plt.switch_backend('Agg')


def transform_cartopy_coord(x, y, from_proj, to_proj):
    """
    Convert between projections using cartopy tranform_points
    """
    tran = to_proj.transform_points(from_proj, x, y)
    return tran[..., 0], tran[..., 1]


def transform_basemap_coord(lon, lat, to_proj):
    """
    Conver Lat Lon values to meters using basemap projection transformation
    """
    tran = np.array([to_proj.projtran(x, y) for x, y in zip(lon, lat)])
    return tran[..., 0], tran[..., 1]


def resize_cbar_function(ax, cax, location='bottom', pad='10%', size='4%'):
    """
    Returns a function to automatically resize the colorbar for cartopy plots

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
                figsize=(10,5), subplot_kw={'projection': ccrs.PlateCarree()})
        >>> cbar_ax = fig.add_axes([0, 0, 0.1, 0.1])

        [... your code generating a scalar mappable ...]

        >>> resize_colorbar = get_resize_event_function(ax, cbar_ax)
        >>> fig.canvas.mpl_connect('resize_event', resize_colorbar)

    Credits
    -------
    Solution by pelson at http://stackoverflow.com/a/30077745/512111
    """
    def resize_colorbar(event):
        # Tell matplotlib to re-draw everything, so that we can get
        # the correct location from get_position.
        plt.draw()
        pos = ax.get_position()
        if location is 'bottom':
            cax_pos = [
                pos.xmin, pos.ymin - (pct2num(pad) + pct2num(size) / 1.5),
                pos.width, pct2num(size) / 1.5]
        if location is 'top':
            cax_pos = [
                pos.xmin, pos.ymax + (pct2num(pad)),
                pos.width, pct2num(size) / 1.5]
        elif location is 'left':
            cax_pos = [
                pos.xmin - pct2num(pad), pos.ymin,
                pct2num(size) / 1.5, pos.height]
        elif location is 'right':
            cax_pos = [
                pos.xmin + pos.width + pct2num(pad), pos.ymin,
                pct2num(size) / 1.5, pos.height]

        cax.set_position(cax_pos)

    resize_colorbar(None)
    return resize_colorbar


def pct2num(string):
    """
    Convert percent string to number
    """
    return float(string.strip('%')) / 100


def dialog_pickfiles(initial_dir=None, extension=None):
    """
    Dialog picker for one or more files.

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

    file_type = zip(['file type'] * len(extension), extension)
    files = askopenfilenames(
        initialdir=initial_dir, filetypes=file_type, title='Pick file[s]')

    if files:
        out = [f for f in files] if len(files) > 1 else files[0]
        return out
    else:
        print 'Warning! Nothing picked, None returned.'


def dict2nc(dictionary, ncfile, zlib=True, attr=None):
    """
    Convert a dictionary to simple netCDF (nf4 with compression enabled).

    Warning: it assumes all datasets in the dictionary are 1D and of equal
    length. Do not use this function for complex multi-dimensional arrays.

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
        f.createDimension('N', len(dictionary.values()[0]))

        # create variables from keys and fill data
        for key, value in dictionary.items():
            var = f.createVariable(key, 'f4', ('N',), zlib=zlib)
            var[:] = value

        # add global attributes to file
        f.history = "[{}] Created netCDF4 zlib={} dataset.".format(
            datetime.today().strftime('%F %T'), zlib)
        if attr is not None:
            for key, value in attr.items():
                f.temp = value
                f.renameAttribute('temp', key)


def seqdate(start_date, end_date, in_fmt='%Y-%m-%d', out_fmt='%Y%m%d'):
    d1 = datetime.strptime(start_date, in_fmt)
    d2 = datetime.strptime(end_date, in_fmt)
    delta = d2 - d1
    return [(d1 + timedelta(days=i)).strftime(out_fmt) for
            i in range(delta.days + 1)]


def modify_cmap(index_list, rgba_list, N=None, name='Custom Cmap'):
    """
    Modify rgba values of an existing colormap.

    Parameters
    ----------
    index_list : list of ints
        list of position in colour table where the values need customisation
    rgba_list : list of tuples
        list of input (r,g,b,a) tuples corresponding to index
    N : int, optional
        the number of entries in the map. The default is None, in which case
        there is one colormap entry for each element in the list of colours
    name : str, optional
        a string to identify the colormap (default 'Custom Cmap')

    Returns
    -------
    matplotlib.colors.ListedColormap object

    Examples
    --------
        >>> my_cmap = modify_cmap([0], [(255, 255, 255, 0.5)])

    """
    cm = mpl.cm.get_cmap()  # get default colormap
    ccl = [cm(i) for i in range(cm.N)]  # build custom color list
    for i, _ in enumerate(index_list):
        ccl[index_list[i]] = rgba_list[i]  # modify colors

    return mpl.colors.ListedColormap(ccl, name, N)


class List(list):
    """
    Create a custom List object to which user attributes can be added.

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
    """Discrete triplet data (x, y, z) analyser"""
    def __init__(self, x, y, z):
        """
        XYZ Constructor

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
        # bin parameters for griddata
        self.delta = (1, 1)
        self.limit = [[-180., 180.], [-90., 90.]]
        # plot parameters used in mapdata
        self.figsize = (8, 5)
        self.figheight = self.figsize[1]
        self.gspacing = (30, 30)
        self.xlab = [0, 0, 0, 1]
        self.ylab = [1, 0, 0, 0]
        self.cbpad = '10%'
        self.G, self.xc, self.yc = None, None, None

    def __enter__(self):
        return self

    def _get_extent(self):
        """get extent for x and y values"""
        return [[self.x.min(), self.x.max()], [self.y.min(), self.y.max()]]

    def _get_latts(self, spacing):
        """get latitude grid mark locations"""
        return np.arange(-90., 99., spacing)

    def _get_lonts(self, spacing):
        """get longitude grid mark locations"""
        return np.arange(-180., 181., spacing)

    def bindata(self, **kw):
        """
        Bin irregular 1D data (triplets) on to 2D plane.

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

    def griddata(self, method='mean', **kw):
        """
        Compute a bidimensional binned statistic for one or more sets of data.

        This is a generalization of a histogram2d function implemented in
        bindata. A histogram divides the space into bins, and returns the
        count of the number of points in each bin. This function allows the
        computation of the sum, mean, median, or other statistic of the values
        (or set of values) within each bin.

        Parameters
        ----------
        method : str, optional
            The statistic to compute (default is 'mean'). Available statistics
            are: 'mean', 'median', 'count', 'sum', 'min', 'max'
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
            equivalent to method='mean'
        """
        delta = kw.get('delta', self.delta)
        globe = kw.get('globe', False)
        order = kw.get('order', False)
        limit = kw.get('limit', self._get_extent())
        # parse method parameters
        if method.lower() == 'total':
            method = 'sum'
        if method.lower() == 'avg' or method.lower() == 'average':
            method = 'mean'

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
            self.x, self.y, self.z, statistic=method.lower(), bins=[xs, ys])
        # bin edges to bin centres
        xc = (stat4[1][1:] + stat4[1][:-1]) / 2
        yc = (stat4[2][1:] + stat4[2][:-1]) / 2
        # return masked (masks for NaNs) array
        if method == 'sum':
            count4 = binned_statistic_2d(
                self.x, self.y, self.z, statistic='count', bins=[xs, ys])
            G = np.ma.masked_where(count4[0] == 0, stat4[0])
        else:
            G = np.ma.masked_invalid(stat4[0])
        #
        self.G, self.xc, self.yc = G, xc, yc
        if order is True:
            # return np.flipud(np.rot90(G)), xc, yc
            self.G = G.transpose()
            return self.G, self.xc, self.yc
            # return G.transpose(), xc, yc
        return self.G, self.xc, self.yc
        # return G, xc, yc

    def plot(self, **kw):
        """Plot (x, y, z) series"""
        _f, ax = plt.subplots(1, figsize=self.figsize)
        ax.plot(self.x.ravel(), c='r', label='x', **kw)
        ax.plot(self.y.ravel(), c='g', label='y', **kw)
        ax.plot(self.z.ravel(), c='b', label='z', **kw)
        plt.legend(ncol=3, borderaxespad=0)
        plt.tight_layout()
        return plt

    def scatter(self, **kw):
        """
        Scatter plot (x, y, z) data triplets

        Parameters
        ----------
            see matplotlib scatter.

        """
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
        """
        Hexagon plot (x, y, z) data triplets.
        hexbin is a bivariate histogram in which the xy-plane is tessellated
        by a regular grid of hexagons.

        """
        _f, ax = plt.subplots(1, figsize=self.figsize)
        im = ax.hexbin(
            self.x.ravel(), self.y.ravel(), C=self.z.ravel(), bins=None, **kw)
        plt.colorbar(im)
        plt.tight_layout()
        return plt

    def mapdata(self, use_cartopy=False, **kw):
        """
        Render xyz data on map.

        Parameters
        ----------
        use_cartopy : bool, optional
            Force maps to use cartopy in place of default Basemap.
        cb_extend : str, optional
            extend colorbar pointed ends (default: 'both').
            Accepted values are 'neither'|'both'|'min'|'max'.
        cb_on : bool, optional
            add colorbar on map (default False)
        cb_loc : str, optional
            location of colorbar position (default 'bottom'). available
            options are 'left'|'right'|'top'|'bottom'.
        cb_pad : str, optional
            paddding of colorbar from main plot in % (default '10%').
        cb_title : str, optional
            title for colorbar title (default None).
        cbt_ha : str, optional
            colorbar title horizontal alignment (default 'center').
        cbt_size : str, optional
            colorbar size (default '4%').
        cbt_pos : tuple, optional
            colorbar title (x, y) position in normal coordinate
            (default (0.9, 0.75) for all methods (1, 0.75) for counts).
        clip : bool, optional
            clip map extent within valid geographic ranges (default False).
        cmap : str, optional
            matplotlib colormap name to use (default "Spectral_r").
        delta : float or (float, float), optional
            resolution specs for binning original data in x and y direction
            (default (1, 1).
        drawcountries : bool, optional
            draw country boundaries (default False).
        figsize : (number, number), optional
            output figure size (default auto adjusted with a base height of 5).
        figheight : number, optional
            control output figure base height (default 5)
        gcol : str, optional
            colors for parallel and meridian lines (default 'gray').
        gline : (number, number), optional
            grid line style that is dash pattern for meridians and parallels
            (default (None, None) for solid line). For example,
            ``gline = (1, 1)`` will draw 1 pixel on, 1 pixel off.
            Available for Basemap version only.
        gtextcolor : str, optional
            gridline text color (default '#333333') available for cartopy only.
        gtfamily : str, optional
            gridline text family (default monospace).
        gtextsize : number, optional
            gridline text size (default 7).
        globe : bool, optional
            set map extent to global extent (180W:180E, 90S:90N) (default
            False). This keyword overrides limit values.
        gspacing : (number, number), optional
            spacing between grid lines (default (30, 30)).
        limit : [[float,float],[float,float]], optional
            output map extent geographic values [[lon0,lon1], [lat0,lat1]]
            (default is calculated from min max values).
        map_res : str, optional
            coastline resolution on map (default 'c' or coarse). available
            options are: 'c'oarse|'l'ow|'h'igh.
        map_buffer : float, optional
            adds extra buffer to map_limit (default None).
        markersize : int, optional
            marker size for scatter points (default 20).
        method : str, optional
            statistic to use for data binning before plotting (default 'mean')
            see griddata. Not used when ``plt_type='scatter'``.
        plt_type : str, optional
            plot type option (default 'pcolormesh'). available options are
            'hexbin'|'pcolormesh'|'scatter'.
        projection : str, optional
            output map projection (default 'cyl'indrical).
        title : str, optional
            title string on map (default is no title).
        vmax : float, optional
            maximum value for plot scaling (default max(data)).
        vmin : float, optional
            minimum value for plot scaling (default min(data))
        verbose : bool, optional
            verbose mode (default False).

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

        stat = kw.get('method', 'mean')
        plotype = kw.get('plt_type', 'pcolormesh')
        ms = kw.get('markersize', 20)
        cbaron = kw.get('cb_on', False)
        cbloc = kw.get('cb_loc', 'bottom')
        cbpad = kw.get('cb_pad', self.cbpad)
        cbextend = kw.get('cb_extend', 'neither')
        cbsz = kw.get('cb_size', '4%')
        cbtitle = kw.get('cb_title', stat)
        auto_cbtha = ['right', 'center'][cbextend in ('max', 'both')]
        cbtha = kw.get('cbt_ha', auto_cbtha)
        cbtpos = kw.get('cbt_pos', [(0.5, 0.75), (1, 0.75)][cbtha is 'right'])
        cmap = kw.get('cmap', 'Spectral_r')
        delta = kw.get('delta', self.delta)
        drawcountries = kw.get('drawcountries', False)
        figsize = kw.get('figsize', self.figsize)
        fight = kw.get('figheight', self.figheight)
        gcol = kw.get('gcol', 'gray')
        gline = kw.get('gline', (None, None))
        gspacing = kw.get('gspacing', self.gspacing)
        gtc = kw.get('gtextcolor', '#333333')
        gtf = kw.get('gtfamily', 'monospace')
        gts = kw.get('gtextsize', 7)
        globe = kw.get('globe', False)
        limit = kw.get('limit', [[self.x.min(), self.x.max()],
                                 [self.y.min(), self.y.max()]])
        mbuf = kw.get('map_buffer', None)
        if mbuf:
            limit += np.array([[-mbuf, mbuf], [-mbuf, mbuf]])

        if globe is True:
            limit = self.limit

        clip = kw.get('clip', False)
        mres = kw.get('map_res', 'c')
        title = kw.get('title', ' ')
        vmin = kw.get('vmin', self.z.min())
        vmax = kw.get('vmax', self.z.max())
        verb = kw.get('verbose', False)

        if plotype.lower() in ('hexbin', 'pcolormesh'):
            if self.G is None:
                self.G, self.xc, self.yc = self.griddata(
                    method=stat, delta=delta, order=True)

        if stat.lower() is 'count':
            vmin, vmax = self.G.min(), self.G.max()

        if clip is True:
            limit[0] = np.clip(limit[0], -180., 180.)
            limit[1] = np.clip(limit[0], -90., 90.)

        clon = np.mean(limit[0])
        clat = np.mean(limit[1])
        dx = np.ptp(np.array(limit[0], 'f'))
        dy = np.ptp(np.array(limit[1], 'f'))
        aspect = dx / dy
        yoff = [1, 0.5][cbloc in ('left', 'right')]
        figsize = kw.get('figsize', (fight * aspect, fight + yoff))
        if verb:
            print 'figsize={}'.format(figsize)

        # draw map
        fig = plt.figure(figsize=figsize)

        if Basemap is None or use_cartopy is True:  # ------------- use cartopy
            data_proj = cartopy.crs.PlateCarree()
            proj = kw.get('projection', data_proj)
            try:
                proj_name = proj.project_geometry.im_class.__name__
            except AttributeError:
                raise AttributeError(
                    '{} is not a valid ccrs proection class'.format(proj))

            font_kw = {'family': gtf, 'size': gts, 'color': gtc}
            gl_kw = {
                'linestyle': '-', 'color': gcol, 'alpha': 0.5,
                'xlocs': self._get_lonts(gspacing[0]),
                'ylocs': self._get_latts(gspacing[1])}
            if mres is 'h':
                mres = '10m'
            elif mres is 'l':
                mres = '50'
            else:
                mres = '110m'

            ax = fig.add_subplot(111, projection=proj, title=title)
            # ax = plt.axes(projection=proj, title=title)
            cb_ax = fig.add_axes([0, 0, 0.1, .1])

            if globe is False:
                ax.set_extent(limit[0] + limit[1], crs=data_proj)
            ax.coastlines(linewidth=0.5, resolution=mres)

            if drawcountries:
                ax.add_feature(
                    cartopy.feature.BORDERS, linewidth=0.4, alpha=0.5)

            if proj_name in ('PlateCarree', 'Mercator'):
                gl = ax.gridlines(draw_labels=True, **gl_kw)
                gl.xlabels_top, gl.ylabels_right = False, False
                gl.xformatter = LONGITUDE_FORMATTER
                gl.yformatter = LATITUDE_FORMATTER
                gl.xlabel_style, gl.ylabel_style = font_kw, font_kw
            else:
                ax.gridlines(**gl_kw)

        else:  # -------------------------------------------------- use Basemap
            fig.add_axes([0.05, 0.05, 0.9, 0.9], title=title)
            proj = kw.get('projection', 'cyl')
            grid_kw = {'family': gtf, 'fontsize': gts, 'color': gcol,
                       'linewidth': 0.2, 'dashes': gline}
            valid_proj = (
                'cyl', 'laea', 'lcc', 'merc', 'mill', 'moll', 'robin', 'sinu')

            if proj not in (valid_proj):
                raise NotImplementedError(
                    '\'{}\' not implemented yet.\n'
                    'Available projections are: {}\n'.format(proj, valid_proj))
            try:
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
            if drawcountries:
                ax.drawcountries(linewidth=0.4, color='#333333')
            ax.drawcoastlines(linewidth=0.5)

        # adjust subplot position based on color bar location
        if cbloc in ('top', 'bottom'):
            orient = 'horizontal'
            bot = [0.05, 0.15][cbloc is 'bottom']
            top = [0.85, 0.9][cbloc is 'bottom']
            cb_tick_kw = {
                'labelsize': 9, 'direction': 'in', 'bottom': True, 'top': True}
            cb_title_kw = {'fontsize': 10, 'ha': cbtha, 'position': cbtpos}
            fig.subplots_adjust(
                bottom=bot, top=top, left=0.075, right=0.95,
                hspace=0, wspace=0)
        elif cbloc in ('left', 'right'):
            orient = 'vertical'
            left = [0.15, 0.075][cbloc is 'right']
            right = [0.95, 0.85][cbloc is 'right']
            cbpad = kw.get('cb_pad', '5%')
            cb_tick_kw = {
                'labelsize': 9, 'direction': 'in', 'left': True, 'right': True}
            cb_title_kw = {
                'fontsize': 10, 'rotation': 'vertical',
                'position': [-0.25, 0.98], 'va': 'top'}
            fig.subplots_adjust(
                bottom=0.15, top=0.85, left=left, right=right,
                hspace=0, wspace=0)

        # add data on map
        if plotype is 'hexbin':
            xx, yy = np.meshgrid(self.xc, self.yc)
            if use_cartopy:
                xxt, yyt = transform_cartopy_coord(
                    xx.ravel(), yy.ravel(), data_proj, proj)
            else:
                xxt, yyt = transform_basemap_coord(xx.ravel(), yy.ravel(), ax)

            im1 = ax.hexbin(
                xxt, yyt, C=self.G.data.ravel(), vmin=vmin, vmax=vmax,
                cmap=cmap)

        elif plotype is 'scatter':
            if use_cartopy:
                xt, yt = transform_cartopy_coord(
                    self.x, self.y, data_proj, proj)
            else:
                xt, yt = transform_basemap_coord(self.x, self.y, ax)
            im1 = ax.scatter(
                xt, yt, c=self.z, lw=0, zorder=5, marker='o', s=ms,
                rasterized=True, cmap=cmap)

        else:  # ## pcolormesh ##
            xx, yy = np.meshgrid(self.xc, self.yc)
            if use_cartopy:
                xxt, yyt = transform_cartopy_coord(xx, yy, data_proj, proj)
            else:
                xxt, yyt = transform_basemap_coord(xx.ravel(), yy.ravel(), ax)
                xxt = np.reshape(xxt, (len(self.yc), len(self.xc)))
                yyt = np.reshape(yyt, (len(self.yc), len(self.xc)))

            im1 = ax.pcolormesh(
                xxt, yyt, self.G, vmin=vmin, vmax=vmax, rasterized=True,
                shading='flat', cmap=cmap)

        if cbaron:
            if Basemap is None or use_cartopy is True:
                # cb = plt.colorbar(
                #     im1, orientation='horizontal', aspect=40,
                #     pad=pct2num(cbpad), shrink=0.9, format='%g',
                #     extend=cbextend)
                def_cb_pad = ['7%', '5%'][cbloc in ('left', 'right')]
                cbpad = kw.get('cb_pad', def_cb_pad)
                fig.canvas.mpl_connect(
                    'resize_event', resize_cbar_function(
                        ax, cb_ax, location=cbloc, pad=cbpad, size=cbsz))

                cb = plt.colorbar(
                    im1, orientation=orient, cax=cb_ax, format='%g',
                    extend=cbextend)
            else:
                cb = ax.colorbar(
                    im1, location=cbloc, size=cbsz, pad=cbpad, format='%g',
                    extend=cbextend)

            cb.ax.tick_params(**cb_tick_kw)
            cb.ax.set_title('{}'.format(cbtitle), **cb_title_kw)

            fig.canvas.draw_idle()
        return plt

    def __exit__(self, exc_type, exc_value, traceback):
        return isinstance(exc_value, TypeError)


class Convert(object):
    """Simple unit converter."""
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
        """
        Solar zenith angle (degrees) to Airmass.

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


def mkdirp(directory):
    """
    Make parent directories as needed, no error if existing.

    Parameters
    ----------
    directory : str
        name of directory (tree) to create
    """
    try:
        os.makedirs(directory)
    except OSError as err:
        if err.errno != errno.EEXIST:
            raise err


def symlinkf(fsrc, fdst):
    """
    Create a symbolic link pointing to fsrc named fdst.
    This is a os.symlink wrapper which removes existing destination files

    Availability: Unix.

    Paremeters
    ----------
    fsrc : str
        source file
    fdst : str
        symbolic link
    """
    try:
        os.symlink(fsrc, fdst)
    except OSError, err:
        if err.errno == errno.EEXIST:
            os.remove(fdst)
            os.symlink(fsrc, fdst)


def purge(dsrc, ext):
    """
    Delete all files with matching extension in a directory.

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


def et2dt(elapsed_seconds, start_date, date_format='%Y-%m-%d'):
    """
    Get list of datetime objects for elapsed seconds from a given start date.

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


def dt2et(date_time, start_date, date_format='%Y-%m-%d %H:%M:%S'):
    """
    Get list of datetime objects for elapsed seconds from a given start date.

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
                   "2010-01-01 00:00:00")
            [172800.0, 208800.0]
    """
    dt = [date_time] if np.isscalar(date_time) else date_time
    return [(datetime.strptime(d, date_format) -
             datetime.strptime(start_date, date_format)).total_seconds()
            for d in dt]


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

    See also
    --------
    datetime format codes
    https://docs.python.org/2/library/datetime.html#strftime-and-strptime-
    behavior

    """
    if in_date is None:
        in_date = datetime.utcnow().strftime(in_format)
    in_date = [in_date] if np.isscalar(in_date) else in_date

    return [datetime.strptime(d, in_format).strftime(out_format) for
            d in in_date]


def dt2yjd(in_date, in_fmt='%Y%m%d', sep=''):
    """
    Convert Date to Year and Dayofyear.

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


def doy(year=None, month=None, day=None):
    """Calculate serial day of year"""
    if day is not None and month is not None and year is not None:
        dt = date(year, month, day)
        return dt.timetuple().tm_yday
    else:
        return date.today().timetuple().tm_yday


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


def locate(array, value, epsilon=1.0):
    """Locate the index of a value in an array."""
    inarray = np.asanyarray(array, dtype=array.dtype)
    loc = (np.abs(value - inarray)).argmin()
    if np.abs(value - inarray[loc]) < epsilon:
        return loc
    else:
        return -1


def v_locate(array, value):
    """
    v_locate function finds the intervals within a given monotonic vector
    that brackets a given set of one value. Much faster than locate when the
    array is sorted (or monotonically increasing).
    """
    idx = np.searchsorted(array, value, side="left")
    if math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx]):
        return idx - 1  # array[idx - 1]
    else:
        return idx  # array[idx]


def getnn(d1, d2, r, k=5):
    """
    Search nearest neighbors between two coordinate catalogues. See
    https://docs.scipy.org/doc/scipy/reference/generated/
    scipy.spatial.cKDTree.query.html

    Parameters
    ----------
    d1, d2 : array_like, last dimension self.m
        Arrays of points to query (in (n, 2)).

    r : nonnegative float
        Return only neighbors within this distance. This is used to prune
        tree searches, so if you are doing a series of nearest-neighbor
        queries, it may help to supply the distance to the nearest neighbor of
        the most recent point.

    k : list of integer or integer
        The list of k-th nearest neighbors to return. If k is an integer
        it is treated as a list of [1, ... k] (range(1, k+1)). Note that
        the counting starts from 1.

    Returns
    -------
    idx : ndarray of ints
        The locations of the neighbors in self.data. If x has shape
        tuple+(self.m,), then i has shape tuple+(k,). When k == 1, the last
        dimension of the output is squeezed. Missing neighbors are indicated
        with self.n.
    d : array of floats
        The distances to the nearest neighbors. If x has shape tuple+(self.m,),
        then d has shape tuple+(k,). When k == 1, the last dimension of the
        output is squeezed. Missing neighbors are indicated with infinite
        distances.

    Example:
        >>> Lon1 = numpy.random.random(2000)
        >>> Lat1 = numpy.random.random(2000)
        >>> Lon2 = numpy.random.random(20)
        >>> Lat2 = numpy.random.random(20)
        >>> d1 = numpy.asarray(zip(Lon1, Lat1))
        >>> d2 = numpy.asarray(zip(Lon2, Lat2))
        >>> i, d = getnn(d1, d2, 0.1, k=3)
    """
    t = cKDTree(d1)
    d, idx = t.query(d2, k=k, eps=0, p=2, distance_upper_bound=r)
    return idx, d


def write_csv(filename, data=None, header=None, append=False):
    mode = 'ab' if append is True else 'wb'
    with open(filename, mode) as fp:
        writer = csv.writer(fp)
        # write header row
        if append is False:
            try:
                writer.writerow(header)
            except UserWarning:
                print 'No header.'
        # write row row(s)
        try:
            writer.writerow(data)
        except UserWarning:
            print "No data."


def ext2vis(ext):
    """
    Visibility V, as a function of the extinction coefficient ext (using
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


def aot2pm25(aot, H, S):
    # AOD and PM2.5 relationship (Hoff and Christopher, 2009)
    # AOT = PM2.5*H*S
    # H: the height of the well-mixed planetary boundary layer (PBL) and
    # S: is the specific extinction efficiency of the aerosol at the ambient
    # relative humidity.
    return aot / (H * S)


def doy2ymd(year, dayofyear):
    """Get date object from year and seq day of year."""
    return datetime(year, 1, 1) + timedelta(dayofyear - 1)


def freq2wl(freq):
    """Convert frequency (cm-1) to wavelength (microns)."""
    return 1e4 / freq


def wl2freq(wl):
    """Convert wavelength (microns) to frequency (cm-1)."""
    return 1e4 / wl
