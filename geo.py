#!/usr/bin/env python2.7
"""
Created on 29 Jun 2015

:author: fra6 (yaswant.pradhan)
:copyright: Crown copyright. Met Office
"""
import numpy as np
import cartopy.crs as ccrs


# Common constants
IMDI = -32768
RMDI = -1073741824.0
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
                print name


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
        from ypylib.hdf import h5Parse
        self.filename = filename
        d = h5Parse(filename)
        data = d.get_data(dsname)[::-1, ::-1][::stride[1], ::stride[0]]
        return data

    def extract_geo(self, msgarray, **kw):
        # TODO: incomplete
        tlLonLat = kw.get('tlLonLat', (0.0, 0.0))  # top-left Lon, Lat
        brLonLat = kw.get('brLonLat', (0.0, 0.0))  # bottom-right Lon, Lat
#         trLonLat = kw.get('trLonLat', (0.0, 0.0))
#         blLonLat = kw.get('blLonLat', (0.0, 0.0))
        tl = self.geo2pix(tlLonLat[0], tlLonLat[1])
        br = self.geo2pix(brLonLat[0], brLonLat[1])
#         tr = self.geo2pix(trLonLat[0], trLonLat[1])
#         bl = self.geo2pix(blLonLat[0], blLonLat[1])
        print tl, br
# msg = MSG


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
    h5f = '/data/local/fra6/sample_data/slotstores/MSG_201202151600_lite.h5'
    msg = MSG()
#     f = msg.geo2pix([0, 2], [0, 10])
#     print f[0], f[1]
    print msg.pix2geo([1856, 600, 100], [1857, 600, 100])
    x = msg.geo2pix([-20, 0, 10, 20], [0, 0, 0, 0])
    print x[0]._data, x[1]._data

    # print h.extract(h5f, 'MSG/Ch10/Raw', stride=(4, 4)).shape
