#!/usr/bin/env python2.7
"""
Created on 29 Jun 2015
:author: yaswant.pradhan
"""
from __future__ import print_function
import os
import numpy as np
import cartopy.crs as ccrs
from ypylib.hdf import h5Parse
from ypylib import log


# Common constants
IMDI = -32768
RMDI = -1073741824.0
RMDItol = 1e-3
RPOL_EARTH = 6356.7523  # radius in km from earth centre to pole
REQ_EARTH = 6378.1370  # radius in km from earth centre to equator
rpolsq_reqsq = RPOL_EARTH ** 2 / REQ_EARTH ** 2
reqsq_rpolsq = REQ_EARTH ** 2 / RPOL_EARTH ** 2
rsqdiff_reqsq = (REQ_EARTH ** 2 - RPOL_EARTH ** 2) / REQ_EARTH ** 2


def find_projections():
    """Cartopy projection generator"""
    for obj_name, o in vars(ccrs).copy().items():
        if isinstance(o, type) and issubclass(o, ccrs.Projection) and \
           not obj_name.startswith('_') and obj_name not in ['Projection']:

            yield o


def get_ccrs_list():
    """Get a list of all available projections in cartopy"""
    return sorted([proj.__name__ for proj in find_projections()])


def print_cartopy_projections():
    """Print all available projections in cartopy"""
    skip_cls = ('_', 'CRS', 'Geocentric', 'Geodetic', 'Projection',
                'RotatedGeodetic')
    for name in dir(ccrs):
        if not name.startswith(skip_cls):
            el = getattr(ccrs, name)
            if isinstance(el, type) and issubclass(el, ccrs.CRS):
                print(name)


class MSG:
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

    """
    satellite = 'MSG'

    def __init__(self):
        """Meteosat Second Generation specific values/constants."""
        self.SAT_ALT = 42164.0  # distance from earth centre to satellite
        self.SUB_LON = 0.0  # longitude of sub-satellite point
        self.COFF = 1856.0  # Column offset (see note above)
        self.LOFF = 1856.0  # Line offset (see note above)
        self.CFAC = -781648343.0  # scaling coefficients (see note above)
        self.LFAC = -781648343.0  # scaling coefficients (see note above)
        self.sd_coeff = self.SAT_ALT ** 2 - REQ_EARTH ** 2  # 1737122264
        self.filename = None

    def geo2pix(self, lon, lat):
        """
        Returns the pixel column and row number of an MSG image for a
        given pair of longitude, latitude values.
        Note: calculation based on the formulae given in Ref [1]

        Parameters
        ----------
        lon : array_like
            Longitude value(s)
        lat : arraye_like
            Latitude value(s)

        Returns
        -------
        tuple
            retruns a tuple of (col, row) where col and row are masked arrays

        Example
        -------
        >>> from ypylib.geo import MSG
        >>> x, y = MSG().geo2pix(0, 0)
        >>> print x, y
        1856 1856
        """
        lonR = np.radians(lon)
        latR = np.radians(lat)

        # Calculate the geocentric latitude from the geographic one using
        # equations on page 24, Ref [1]
        # c_lat = np.arctan(0.993243 * np.tan(latR))

        c_lat = np.arctan(rpolsq_reqsq * np.tan(latR))

        # Using c_lat calculate the length form the earth centre to the
        # surface of the earth ellipsoid ;equations on page 24, Ref. [1]
        # re = rpol / (sqrt(1 - (req^2-rpol^2)/req^2 * cos^2(c_lat) )

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
        Return the longitude and latitude value of a pixel on MSG disc
        given the column and row index of that pixel.

        Parameters
        ----------
        col : array_like
            column index of MSG/SEVIRI pixel
        row : arraye_like
            row index of MSG/SEVIRI pixel
        dtype : dtype
            data type precision for lat and lon arrays (default 'float64')

        Returns
        -------
        tuple
            retruns a tuple of (lon, lat) where lon and lat are masked arrays

        Example
        -------
        >>> from ypylib.geo import MSG
        >>> lon, lat = MSG().pix2geo(1856, 1856)
        >>> print lon, lat
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

        Parameters
        ----------
        fd_array : 2D array
            SEVIRI full disc array of shape (3712, 3712)
        missing_value : int or float, optional
            values <= missing_value to remove (default None)
        vname : str
            variable name (default 'data')
        dtype : dtype
            data type precision for lat and lon arrays (default 'float64')

        Returns
        -------
        dict
            returns a dictionary with valid data points and associated
            geographic coordinates

        Examples
        --------
            >>> dict = fulldisc2xyz(fdarray, missing=0.0)
            >>> print dict.keys()

        """
        w = np.where(fd_array != missing_value)  # returns tuple (rows, cols)
        lon, lat = self.pix2geo(w[1], w[0], dtype=xydtype)
        # print lon.data.dtype, lat.data.dtype
        return {vname: fd_array[w], 'lon': lon.data, 'lat': lat.data}

    def extract(self, filename, dsname, stride=(1, 1)):
        """
        Extract data from a h5 slotstore dataset.
        """
        self.filename = filename
        d = h5Parse(filename)
        data = d.get_data(dsname)[::-1, ::-1][::stride[1], ::stride[0]]
        return data

    # def extract_geo(self, msgarray, **kw):
    #     # TODO: incomplete
    #     tlLonLat = kw.get('tlLonLat', (0.0, 0.0))  # top-left Lon, Lat
    #     brLonLat = kw.get('brLonLat', (0.0, 0.0))  # bottom-right Lon, Lat
    #     tl = self.geo2pix(tlLonLat[0], tlLonLat[1])
    #     br = self.geo2pix(brLonLat[0], brLonLat[1])
    #     # trLonLat = kw.get('trLonLat', (0.0, 0.0))  # top-right Lon, Lat
    #     # blLonLat = kw.get('blLonLat', (0.0, 0.0))  # bottom-left Lon, Lat
    #     # tr = self.geo2pix(trLonLat[0], trLonLat[1])
    #     # bl = self.geo2pix(blLonLat[0], blLonLat[1])
    #     print(tl, br)
# msg = MSG


class Plot(object):

    """docstring for ClassName"""
    def __init__(self, h5file, satellite='MSG', **kw):

        super(Plot, self).__init__()
        self.h5file = h5file
        self.satellite = satellite.upper()
        self.stride = kw.get('stride', (1, 1))
        self.dust_quality = kw.get('dust_quality', None)
        self.MDI = None,
        self.projection = kw.get('projection', None)
        self.extent = kw.get('extent', None)
        self.coastline = kw.get('coastline', True)
        self.cres = kw.get('cres', '50m')
        self.ccolor = kw.get('ccolor', 'w')
        self.drawcountries = kw.get('drawcountries', False)
        self.cbar = kw.get('cbar', False)
        self.cb_width = kw.get('cb_width', 0.02)
        self.cb_unit = kw.get('cb_unit', '')
        self.cb_nticks = kw.get('cb_nticks', None)
        self.cb_bounds = kw.get('cb_bounds', None)
        self.quick = kw.get('quick', False)
        self.list_all = kw.get('list_all', False)
        self.show_path = kw.get('show_path', False)
        self.show_ticks = kw.get('show_ticks', False)
        self.xglocs = kw.get('xglocs', None)
        self.yglocs = kw.get('yglocs', None)
        self.gl_font = kw.get(
            'gl_font', dict(family='monospace', size=8, color='#333333'))
        self.glw = kw.get('glw', 0.5)
        self.tag_fig = kw.get('tag_fig', None)
        self.tag_col = kw.get('tag_col', 'k')
        self.figsize = kw.get('figsize', (6, 6))

    def dust_rgb(self, h5p=None):
        """
        """
        h5p = h5p or h5Parse(self.h5file)

        if self.satellite == 'HIM8':
            bnd1 = '/HIM8/B11/BT'  # 8.6: 8.44--8.76
            bnd2 = '/HIM8/IR1/BT'  # 10.4: 10.3--10.6
            bnd3 = '/HIM8/IR2/BT'  # 12.3: 12.2--12.5
        elif self.satellite in ('GOES16', 'GOES-E'):
            bnd1 = '/GOES16/Ch11/BT'  # 8.5: 8.7--8.7
            bnd2 = '/GOES16/Ch13/BT'  # 10.35: 10.1--10.6
            bnd3 = '/GOES16/Ch15/BT'  # 12.3: 11.8--12.8
        elif self.satellite in ('GOES17', 'GOES-W'):
            bnd1 = '/GOES16/Ch11/BT'  # 8.5: 8.7--8.7
            bnd2 = '/GOES16/Ch13/BT'  # 10.35: 10.1--10.6
            bnd3 = '/GOES16/Ch15/BT'  # 12.3: 11.8--12.8
        elif self.satellite in ('MSG', 'IODC', 'MSG_IODC'):
            bnd1 = '/MSG/IR_087/BT'  # 8.7: 8.3--9.1
            bnd2 = '/MSG/IR_108/BT'  # 10.8: 9.8--11.8
            bnd3 = '/MSG/IR_120/BT'  # 12.0: 11.0--13.0
        else:
            raise NotImplementedError('%s not implemented yet.',
                                      self.satellite)

        bt_870 = h5p.get_data(bnd1)[bnd1]
        bt_108 = h5p.get_data(bnd2)[bnd2]
        bt_120 = h5p.get_data(bnd3)[bnd3]

        # Reverse and sub-sample arrays
        bt_870 = bt_870[::-1, ::-1][::self.stride[0], ::self.stride[1]]
        bt_108 = bt_108[::-1, ::-1][::self.stride[0], ::self.stride[1]]
        bt_120 = bt_120[::-1, ::-1][::self.stride[0], ::self.stride[1]]

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

    def geo_extent(self, msg_datetime, satellite='MSG'):
        """
        Get MSG Geostationary projection extent.
        Note: there is 1 pixel offset before/after 6 Dec 2017

        Parameters
        ----------
        msg_datetime : str
            MSG datetime string in CCYYmmddHHMM format

        Returns
        -------
        tuple
            MSG extent in X and Y direction (x0, x1, y0, y1)

        """
        from SpsMod_Coordinates import GeoProjection
        return GeoProjection(self.satellite).as_cartopy_crs().extent

    def show(self, dataset=None, dust_rgb=False, save_fig=None, **ims_kw):
        import matplotlib.pyplot as plt
        from matplotlib import colors, ticker
        from cartopy.mpl import gridliner
        import cartopy.feature as cfeature

        stride = self.stride
        idx = 1  # index of datatime string in hdf filename
        c_lon = 0  # central geostationary longitude
        def_extent = (-80, 80, -80, 80)  # default map extent
        sat = self.satellite
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

        h5 = h5Parse(self.h5file)
        if self.list_all:
            h5.ls()

        if self.quick:
            stride = (10, 10)

        apply_dust_quality = False
        if self.dust_quality:
            dconf = '/Product/GM/DustAOD/DustConfidence'
            try:
                dset = h5.get_data(dconf)
                dust_qual = dset[dconf][::-1, ::-1][::stride[0], ::stride[1]]
                apply_dust_quality = True
            except KeyError:
                log.warn('Dust_quality not available - filter wont be applied')

        # Get dust RGB array
        if dust_rgb and dataset is None:
            dataset = 'DustRGB'
            cbar = False
            plot_array = self.dust_rgb(h5)

        elif dataset and dust_rgb is False:
            dset = h5.get_data(dataset)
            plot_array = dset[dataset][::-1, ::-1][::stride[0], ::stride[1]]
            # print plot_array.dtype

            if plot_array.dtype in (
                    np.int8, np.int16, np.intc, np.uint8, np.uint16, np.intp):
                self.MDI = [self.MDI, IMDI][self.MDI is None]
                plot_array = np.ma.masked_equal(plot_array, self.MDI)
                if apply_dust_quality:
                    plot_array = np.ma.masked_less(
                        dust_qual, self.dust_quality)
            else:
                self.MDI = [self.MDI, RMDI][self.MDI is None]
                plot_array[plot_array <= (self.MDI + RMDItol)] = np.nan
                if apply_dust_quality:
                    plot_array = np.ma.where(
                        dust_qual >= self.dust_quality, plot_array, np.nan)

        elif dataset and dust_rgb:
            raise Exception('Cant plot dust_rgb and other dataset together.')
        else:
            raise Exception(
                'Dont know what to plot, please specify dataset path.\n\n' +
                '  msgview(filename, [dataset_path|dust_rgb=True], .. OR\n' +
                '  msgview(filename, list_all=True) to see all datasets.\n')
            return

        # Get MSG datetime string from hdf5 filename
        msg_datetime = os.path.basename(self.h5file).split('_')[idx]
        title_str = [msg_datetime, os.path.basename(dataset)]

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

        title = '{}  {}'.format(title_str[0], title_str[1])

        # Start plot
        plt.figure(figsize=self.figsize)

        # normalise colour based on discrete color bounds (cb_bounds)?
        if self.cb_bounds:
            bounds = np.array(self.cb_bounds)
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds))
        else:
            norm = None

        if self.coastline:
            if self.projection in (ccrs.PlateCarree(), ccrs.Mercator()):
                # Plot using rectangular projection
                ax = plt.axes(projection=self.projection)
                extent = [self.extent, def_extent][self.extent is None]
                ax.set_extent(extent, ccrs.PlateCarree())
                gl = ax.gridlines(draw_labels=True, alpha=0.5,
                                  xlocs=self.xglocs, ylocs=self.yglocs,
                                  linestyle='-', linewidth=self.glw)
                gl.xlabels_top = False
                gl.ylabels_right = False
                gl.xformatter = gridliner.LONGITUDE_FORMATTER
                gl.yformatter = gridliner.LATITUDE_FORMATTER
                gl.xlabel_style, gl.ylabel_style = self.gl_font, self.gl_font
                im = ax.imshow(
                    plot_array, transform=ccrs.Geostationary(c_lon),
                    extent=self.geo_extent(
                        msg_datetime, satellite=self.satellite),
                    origin='upper', norm=None, interpolation=None, **ims_kw)
            else:
                # Plot in Geostationary projection
                ax = plt.axes(projection=ccrs.Geostationary(c_lon))
                ax.set_global()
                im = ax.imshow(
                    plot_array,
                    extent=self.geo_extent(
                        msg_datetime, satellite=self.satellite),
                    origin='upper', norm=None, interpolation=None, **ims_kw)

                if self.extent:
                    ax.set_extent(extent, ccrs.PlateCarree())

            ax.coastlines(resolution=self.cres, lw=0.5, color=self.ccolor)
        else:
            if self.extent:
                log.warn('Ignored extent. For subsets set coastline=True')
            ax_pos = [0.05, 0.1, 0.8, 0.8]
            if self.show_ticks:
                ax_pos = [0.1, 0.1, 0.8, 0.8]
            ax = plt.axes(ax_pos)

        # Add country borders
        if self.drawcountries:
            ax.add_feature(cfeature.BORDERS, lw=0.4, edgecolor=self.ccolor)

        pos = ax.get_position()

        # Tag a figure with additional string at top-left (useful for figure
        # panel numbering)
        if self.tag_fig:
            # title = '{} {}'.format(tag_fig, title)
            ax.text(0.01, 0.95, self.tag_fig, color=self.tag_col, fontsize=14,
                    transform=ax.transAxes)

        # Draw title text
        # ax.set_title(title, family='FreeSans', fontsize="14", loc='left')
        ax.set_title(title, fontsize="13", loc='left')

        # Hide plot ticks and tick labels?
        if self.show_ticks is False:
            for ticx, ticy in list(
                zip(ax.xaxis.get_major_ticks(),
                    ax.yaxis.get_major_ticks())):
                ticx.tick1On = ticx.tick2On = False
                ticx.label1On = ticx.label2On = False
                ticy.tick1On = ticy.tick2On = False
                ticy.label1On = ticy.label2On = False

        # Show image without coastline - simple 2D array imshow()
        if self.coastline is False:
            im = plt.imshow(plot_array, norm=norm, interpolation=None,
                            **ims_kw)
            ax.grid(lw=0.5, alpha=0.5)

        # Attach vertical colour bar
        if cbar:
            cax = plt.axes([pos.x1, pos.y0, self.cb_width, pos.height])
            cb = plt.colorbar(mappable=im, cax=cax, format='%g')
            cb.ax.tick_params(direction='in')
            # cb.ax.get_xaxis().labelpad = -115
            cb.ax.set_title('{}'.format(self.cb_unit))
            if self.cb_nticks:
                tick_locator = ticker.MaxNLocator(nbins=self.cb_nticks)
                cb.locator = tick_locator
                cb.update_ticks()
            if title_str[1] in 'DustMask':
                cb.ax.set_yticklabels(
                    ['Dust', 'Ice cloud', 'Low cloud/\nSurface', 'Uncertain'],
                    rotation='90', va='bottom')

        # without colorbar, data value will be shown at cursor position
        # my_imshow(plot_array, ax=ax, interpolation=None, **imshow_kw)
        # plt.imshow(plot_array, interpolation=None, **imshow_kw)

        # Show slotstore file and dataset path on plot?
        if self.show_path:
            textstr = '{}\n{}'.format(dataset, self.h5file)
            plt.gcf().text(0.99, 0.01, textstr, fontsize=8, ha='right')

        # save figure as file to disk?
        if save_fig:
            plt.savefig(save_fig)
            plt.close()
        else:
            plt.show()
            pass


# def dust_rgb(h5p, stride=(1, 1), satellite='MSG'):
#     """
#     Scale channel components for dust RGB (range from 0.0 to 1.0)

#     Parameters
#     ----------
#     h5p : h5Parse object
#         File handler for hdf5 file
#     stride : tuple, optional
#         number of pixels to stride in (x, y) direction

#     Returns
#     -------
#     numpy array
#         Array containing RGB values

#     """
#     sat = satellite.upper()

#     # Get IR brightness temp data arrays for DustRGB
#     if sat == 'HIM8':
#         bnd1 = '/HIM8/B11/BT'  # 8.6: 8.44--8.76
#         bnd2 = '/HIM8/IR1/BT'  # 10.4: 10.3--10.6
#         bnd3 = '/HIM8/IR2/BT'  # 12.3: 12.2--12.5
#     elif sat in ('GOES16', 'GOES-E'):
#         bnd1 = '/GOES16/Ch11/BT'  # 8.5: 8.7--8.7
#         bnd2 = '/GOES16/Ch13/BT'  # 10.35: 10.1--10.6
#         bnd3 = '/GOES16/Ch15/BT'  # 12.3: 11.8--12.8
#     elif sat in ('GOES17', 'GOES-W'):
#         bnd1 = '/GOES16/Ch11/BT'  # 8.5: 8.7--8.7
#         bnd2 = '/GOES16/Ch13/BT'  # 10.35: 10.1--10.6
#         bnd3 = '/GOES16/Ch15/BT'  # 12.3: 11.8--12.8
#     elif sat in ('MSG', 'IODC', 'MSG_IODC'):
#         bnd1 = '/MSG/IR_087/BT'  # 8.7: 8.3--9.1
#         bnd2 = '/MSG/IR_108/BT'  # 10.8: 9.8--11.8
#         bnd3 = '/MSG/IR_120/BT'  # 12.0: 11.0--13.0
#     else:
#         raise NotImplementedError('%s not implemented yet.', satellite)

#     bt_870 = h5p.get_data(bnd1)[bnd1]
#     bt_108 = h5p.get_data(bnd2)[bnd2]
#     bt_120 = h5p.get_data(bnd3)[bnd3]

#     # Reverse and sub-sample arrays
#     bt_870 = bt_870[::-1, ::-1][::stride[0], ::stride[1]]
#     bt_108 = bt_108[::-1, ::-1][::stride[0], ::stride[1]]
#     bt_120 = bt_120[::-1, ::-1][::stride[0], ::stride[1]]

#     # Get RGB channels for dust/ash
#     red = bt_120 - bt_108
#     grn = bt_108 - bt_870
#     blu = bt_108

#     # Component scaling ranges for each channel:
#     smin = np.array([-4.0, 0.0, 261.0])  # Valid minima for [R, G, B]
#     smax = np.array([1.0, 15.0, 289.0])  # Valid maxima for [R, G, B]
#     gama = np.array([1.0, 2.5, 1.0])  # Gamma values for [R, G, B]
#     srng = smax - smin

#     # Assign NaN to out-of-disc pixels
#     invalid = np.where(bt_108 < -100.)
#     red[invalid] = np.NAN
#     grn[invalid] = np.NAN
#     blu[invalid] = np.NAN

#     # Clip RGB arrays to scaled ranges
#     red = np.array(red.clip(smin[0], smax[0]))
#     grn = np.array(grn.clip(smin[1], smax[1]))
#     blu = np.array(blu.clip(smin[2], smax[2]))

#     # Initialise RGB image
#     rgbImg = np.empty((red.shape[0], red.shape[1], 3))

#     # Add scaled RGB components
#     with np.errstate(invalid='ignore'):
#         rgbImg[:, :, 0] = ((red - smin[0]) / srng[0]) ** (1.0 / gama[0])
#         rgbImg[:, :, 1] = ((grn - smin[1]) / srng[1]) ** (1.0 / gama[1])
#         rgbImg[:, :, 2] = ((blu - smin[2]) / srng[2]) ** (1.0 / gama[2])

#     return rgbImg


# def plot(hdf5file,
#          dataset_path=None,
#          dust_rgb=False,
#          satellite='MSG',
#          stride=(1, 1),
#          dust_quality=None,
#          MDI=None,
#          projection=None,
#          extent=None,
#          coastline=True,
#          cres='50m',
#          ccolor='w',
#          drawcountries=False,
#          cbar=False,
#          cb_width=0.02,
#          cb_unit='',
#          cb_nticks=None,
#          cb_bounds=None,
#          quick=False,
#          list_all=False,
#          show_path=False,
#          show_ticks=False,
#          xglocs=None,
#          yglocs=None,
#          gl_font=None,
#          glw=0.5,
#          tag_fig=None,
#          tag_col='k',
#          figsize=(6, 6),
#          save_fig=None,
#          **imshow_kw):
#     """
#     Display 2D array from MSG slotstore file.

#     Parameters
#     ----------
#     hdf5file : str
#         MSG slotstore filename (hdf5).

#     dataset_path : str
#         Full path to 2D array data in hdf5.

#     dust_rgb : bool
#         Plot Dust RGB, require BT at 3 SEVIRI IR window channels in the file.

#     dust_quality : int
#         Mask data based on dust confidence flag.

#     MDI : number
#         Missing Data Indicator.

#     stride : sequence
#         Skip pixels in x,y dimension, Default (1, 1) is to show at full res.

#     projection : cartopy crs
#         Cartopy projection to use for mapping - one of the following
#         PlateCarree, Mercator, or Geostationary (default).

#     extent : sequence
#         Longitude and Latitude bounds to plot the data (lon0, lon1, lat0, lat1)

#     coastline : bool
#         Show coastlines. Default is True. Setting this to False will plot the
#         2D array as is (without geolocation information).

#     cres : str
#         Coastline resolution. Default is '50m'.

#     ccolor : str
#         Coastline colour. Default is white.

#     drawcountries : bool, optional
#         Description

#     cbar : bool
#         Show colour bar. Default is False.

#     cb_width : number
#         Colour bar width. Default is 0.02.

#     cb_unit : str
#         Data unit to show on top of the colour bar.

#     cb_nticks : number
#         Number of ticks to show in the colour bar.

#     cb_bounds : array_like
#         Discrete bounds for colour bar.

#     list_all : bool
#         Print list of all available dataset in the file (no plot).

#     show_path : bool
#         Show slotstore filename and dataset path on plot. Default is False.

#     show_ticks : bool
#         plot ticks and tick labels for non-mapped display.

#     xglocs, yglocs : array_like
#         Locations of x grid-lines and y grid-lines.

#     gl_font : dict, optional
#         Font properties for grid labels default is
#         {'family': 'monospace', 'size': 8, 'color': '#333333'}

#     glw : number
#         Width of grid lines.

#     tag_fig : str
#         Add optional string to top left corner of figure.

#     tag_col : str
#         Colour of tag_fig text.

#     figsize : array_like
#         Figure size. Default is (6, 6)

#     save_fig : str
#         Save plot to a file instead of show. Default is to show plot.
#     **imshow_kw
#         optional keywords for imshow()

#     See also
#     --------
#         See pyplot.imshow() keywords.

#     Returns
#     -------
#     TYPE
#         Description

#     """
#     # defaults for MSG 0-deg service
#     idx = 1  # index of datatime string in hdf filename
#     c_lon = 0  # central geostationary longitude
#     def_extent = (-80, 80, -80, 80)  # default map extent
#     sat = satellite.upper()
#     if sat in ('IODC', 'MSG_IODC'):
#         # update defaults for IODC service
#         idx = 2
#         c_lon = 41.5
#         def_extent = (-40, 120, -80, 80)
#     if sat == 'HIM8':
#         # update defaults for Himawari service
#         c_lon = 140.7
#         def_extent = (60, 220, -80, 80)
#     if sat in ('GOES16', 'GOES-E'):
#         c_lon = -75
#         def_extent = (-155, 5, -80, 80)

#     h5 = h5Parse(hdf5file)
#     if list_all:
#         h5.ls()

#     if quick:
#         stride = (10, 10)

#     apply_dust_quality = False
#     if dust_quality:
#         dconf = '/Product/GM/DustAOD/DustConfidence'
#         try:
#             dset = h5.get_data(dconf)
#             dust_qual = dset[dconf][::-1, ::-1][::stride[0], ::stride[1]]
#             apply_dust_quality = True
#         except KeyError:
#             log.warn('Dust_quality not available - filter wwont be applied')

#     if gl_font is None:
#         gl_font = dict(family='monospace', size=8, color='#333333')

#     # Get dust RGB array
#     if dust_rgb and dataset_path is None:
#         dataset_path = 'DustRGB'
#         cbar = False
#         # plot_array = msg_dustrgb(h5, stride=stride, him8=him8, goes16=goes16)
#         plot_array = msg_dustrgb(h5, stride=stride, satellite=satellite)

#     elif dataset_path and dust_rgb is False:
#         dset = h5.get_data(dataset_path)
#         plot_array = dset[dataset_path][::-1, ::-1][::stride[0], ::stride[1]]
#         # print plot_array.dtype

#         if plot_array.dtype in (
#                 np.int8, np.int16, np.intc, np.uint8, np.uint16, np.intp):
#             MDI = [MDI, IMDI][MDI is None]
#             plot_array = np.ma.masked_equal(plot_array, MDI)
#             if apply_dust_quality:
#                 plot_array = np.ma.masked_less(dust_qual, dust_quality)
#         else:
#             MDI = [MDI, RMDI][MDI is None]
#             plot_array[plot_array <= (MDI + RMDItol)] = np.nan
#             if apply_dust_quality:
#                 plot_array = np.ma.where(
#                     dust_qual >= dust_quality, plot_array, np.nan)

#     elif dataset_path and dust_rgb:
#         raise Exception('Cant plot dust_rgb and other dataset together.')
#     else:
#         raise Exception(
#             'Dont know what to plot, please specify dataset path.\n\n' +
#             '  msgview(filename, [dataset_path|dust_rgb=True], .. OR\n' +
#             '  msgview(filename, list_all=True) to see all datasets.\n')
#         return

#     # Get MSG datetime string from hdf5 filename
#     msg_datetime = os.path.basename(hdf5file).split('_')[idx]
#     title_str = [msg_datetime, os.path.basename(dataset_path)]

#     # Reassign 3-channel DustMask values for readable legend/title
#     if title_str[1] in 'Dust':
#         title_str[1] = 'DustMask'
#         plot_array[plot_array < -0.5] = 0
#         plot_array[plot_array == 300] = 1
#         plot_array[plot_array == 400] = 2
#         plot_array[plot_array == 500] = 3

#     # Rename DustHeightError to DustPressureError
#     if title_str[1] in 'DustHeightError':
#         title_str[1] = 'DustPressureError'

#     title = '{}  {}'.format(title_str[0], title_str[1])

#     # Start plot
#     plt.figure(figsize=figsize)

#     # normalise colour based on discrete color bounds (cb_bounds)?
#     if cb_bounds:
#         bounds = np.array(cb_bounds)
#         norm = colors.BoundaryNorm(boundaries=bounds, ncolors=len(bounds))
#     else:
#         norm = None

#     if coastline is True:
#         if projection in (ccrs.PlateCarree(), ccrs.Mercator()):
#             # Plot using rectangular projection
#             ax = plt.axes(projection=projection)
#             extent = [extent, def_extent][extent is None]
#             ax.set_extent(extent, ccrs.PlateCarree())
#             gl = ax.gridlines(draw_labels=True, alpha=0.5,
#                               xlocs=xglocs, ylocs=yglocs,
#                               linestyle='-', linewidth=glw)
#             gl.xlabels_top = False
#             gl.ylabels_right = False
#             gl.xformatter = gridliner.LONGITUDE_FORMATTER
#             gl.yformatter = gridliner.LATITUDE_FORMATTER
#             gl.xlabel_style, gl.ylabel_style = gl_font, gl_font
#             im = ax.imshow(plot_array,
#                            transform=ccrs.Geostationary(c_lon),
#                            # extent=msg_extent(msg_datetime,
#                            #                   him8=him8,
#                            #                   goes16=goes16),
#                            extent=msg_extent(msg_datetime,
#                                              satellite=satellite),
#                            origin='upper',
#                            norm=None,
#                            interpolation=None, **imshow_kw)
#         else:
#             # Plot in Geostationary projection
#             ax = plt.axes(projection=ccrs.Geostationary(c_lon))
#             ax.set_global()
#             im = ax.imshow(plot_array,
#                            # extent=msg_extent(msg_datetime,
#                            #                   him8=him8,
#                            #                   goes16=goes16),
#                            extent=msg_extent(msg_datetime,
#                                              satellite=satellite),
#                            origin='upper',
#                            norm=None,
#                            interpolation=None, **imshow_kw)

#             if extent:
#                 ax.set_extent(extent, ccrs.PlateCarree())

#         ax.coastlines(resolution=cres, lw=0.5, color=ccolor)
#     else:
#         if extent:
#             log.warn('Ignored extent. For subsets please set coastline=True')
#         ax_pos = [0.05, 0.1, 0.8, 0.8]
#         if show_ticks:
#             ax_pos = [0.1, 0.1, 0.8, 0.8]
#         ax = plt.axes(ax_pos)

#     # Add country borders
#     if drawcountries:
#         ax.add_feature(cfeature.BORDERS, lw=0.4, edgecolor=ccolor)

#     pos = ax.get_position()

#     # Tag a figure with additional string at top-left (useful for figure
#     # panel numbering)
#     if tag_fig:
#         # title = '{} {}'.format(tag_fig, title)
#         ax.text(0.01, 0.95, tag_fig, color=tag_col, fontsize=14,
#                 transform=ax.transAxes)

#     # Draw title text
#     # ax.set_title(title, family='FreeSans', fontsize="14", loc='left')
#     ax.set_title(title, fontsize="13", loc='left')

#     # Hide plot ticks and tick labels?
#     if show_ticks is False:
#         for ticx, ticy in list(
#                 zip(ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks())):
#             ticx.tick1On = ticx.tick2On = False
#             ticx.label1On = ticx.label2On = False
#             ticy.tick1On = ticy.tick2On = False
#             ticy.label1On = ticy.label2On = False

#     # Show image without coastline - simple 2D array imshow()
#     if coastline is False:
#         im = plt.imshow(plot_array, norm=norm, interpolation=None, **imshow_kw)
#         ax.grid(lw=0.5, alpha=0.5)

#     # Attach vertical colour bar
#     if cbar:
#         cax = plt.axes([pos.x1, pos.y0, cb_width, pos.height])
#         cb = plt.colorbar(mappable=im, cax=cax, format='%g')
#         cb.ax.tick_params(direction='in')
#         # cb.ax.get_xaxis().labelpad = -115
#         cb.ax.set_title('{}'.format(cb_unit))
#         if cb_nticks:
#             tick_locator = ticker.MaxNLocator(nbins=cb_nticks)
#             cb.locator = tick_locator
#             cb.update_ticks()
#         if title_str[1] in 'DustMask':
#             cb.ax.set_yticklabels(
#                 ['Dust', 'Ice cloud', 'Low cloud/\nSurface', 'Uncertain'],
#                 rotation='90', va='bottom')

#     # without colorbar, data value will be shown at cursor position
#     # my_imshow(plot_array, ax=ax, interpolation=None, **imshow_kw)
#     # plt.imshow(plot_array, interpolation=None, **imshow_kw)

#     # Show slotstore file and dataset path on plot?
#     if show_path:
#         textstr = '{}\n{}'.format(dataset_path, hdf5file)
#         plt.gcf().text(0.99, 0.01, textstr, fontsize=8, ha='right')

#     # save figure as file to disk?
#     if save_fig:
#         plt.savefig(save_fig)
#         plt.close()
#     else:
#         return plt
#     pass


def wrap_lon(longitude, hemisphere=False):
    """
    Converts Longitude values in an array range from [-180:180] to [0:360]
    (shpere) or vice-versa (hemisphere).

    :param longitude: (sequence) original longitude array
    :param hemisphere: (bool) setting this keyword assumes original array
     ranges from 0 to 360 and therefore returned array will be in
     -180 to 180 range.

    """
    longitude = np.asarray(longitude)
    if hemisphere is True:
        wrapped_lon = ((longitude + 180) % 360) - 180  # range 0:360
    else:
        wrapped_lon = (longitude + 360) % 360  # range -180:+180
    return wrapped_lon


def ll_vec2arr(xv, yv):
    """
    Get Longitude and Latitude 2D arrays from monotonically increasing
    Longitude and Latitude vectors (1D arrays).

    :param xv: (sequence) monotonically increasing x (longitude) array
    :param yv: (sequence) monotonically increasing y (latitude) array

    """
    return np.meshgrid(xv, yv)


if __name__ == "__main__":
    # lon_0 = 180
    # lon = -180 + 5 * np.arange(72)
    # lon_shift = shiftlon(lon, lon_0)
    # print "original lon:", lon
    # print " shifted lon:", lon_shift
    h5f = '$SCRATCH/sample_data/slotstores/MSG_201202151600_lite.h5'
    msg = MSG()
#     f = msg.geo2pix([0, 2], [0, 10])
#     print f[0], f[1]
    print(msg.pix2geo([1856, 600, 100], [1857, 600, 100]))
    x = msg.geo2pix([-20, 0, 10, 20], [0, 0, 0, 0])
    print(x[0]._data, x[1]._data)

    # print h.extract(h5f, 'MSG/Ch10/Raw', stride=(4, 4)).shape
