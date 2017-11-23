"""
:author: yaswant.pradhan
:copyright: Crown copyright. Met Office.
"""

import numpy as np
from scipy import stats
from math import log
from random import randint
# from .utils import v_locate
# from matplotlib import cm
# import matplotlib.pyplot as plt


def stat2(x, y, lognormal=False):
    """
    Statistical comparison between two series.

    For verification of prediction against truth use x as truth and y as
    predicted variable.

    Note:
        - Multi-dimensional input arrays are flattened
        - Invalid numbers are removed

    Parameters
    ----------
    x : array of numbers, reference
        independent or truth variable containing series of points.
    y : array of numbers, predicted
        dependent or predicted variable containing series of points.
        Multi-dimensional arrays will be flattened
    lognormal: bool, optional
        Not implemented

    Returns
    -------
    dict
        dictionary containing statistical metrics (min, max, median, variance,
        skewness, kurtosis) for each series and bias, rmse, urmse,
        fge, slope, stderr for slope, intercept, skill score, correlation,
        p-value, rank correlation and p-value between x and y, and counts for
        number of points used in cross-comparison.

    """
    # TODO
    x = np.array(x).flatten()
    y = np.array(y).flatten()

    # remove invalid points from both series
    valid = np.where(np.isfinite(x) & np.isfinite(y))
    x, y = x[valid], y[valid]

    # describe each series
    # 0: nobs, 1: (min, max), 2: mean, 3: var, 4: skew, 5: kurtosis (Fisher)
    xd = stats.describe(x, axis=None)
    yd = stats.describe(y, axis=None)

    # get spearman's rank-correlation coefficient (and p value) between
    # the two series
    rs_p = stats.spearmanr(x, y)

    # get linear regression (Ordinary Least Square) fit between two series
    # (x being independent). Returns 5 values:
    # 0: slope (gradient), 1: intercept, 2: correlation, 3: p-value,
    # 4: stderr (of the estimated gradient)
    fit = stats.linregress(x, y)

    rmse = np.sqrt(np.mean((y - x)**2))
    bias = np.mean(y - x)

    return {
        'x_min': xd[1][0], 'x_max': xd[1][1], 'x_med': np.median(x),
        'x_avg': xd[2], 'x_var': xd[3], 'x_skew': xd[4], 'x_kurt': xd[5],
        'y_min': yd[1][0], 'y_max': yd[1][0], 'y_med': np.median(y),
        'y_avg': yd[2], 'y_var': yd[3], 'y_skew': yd[4], 'y_kurt': yd[5],
        'bias': bias,                           # mean bias
        'rmse': rmse,                           # RMSE
        'urmse': np.sqrt(rmse**2 - bias**2),    # unbiased RMSE
        'fge': 2.0 * np.mean(np.abs((y - x) / (y + x))),
        'r_spearman': rs_p[0],          # Rank correlation (Spearman)
        'p_spearman': rs_p[1],          # p-value for Rank correlation
        'slope': fit[0],                # Slope of linear fit
        'intercept': fit[1],            # intercept og linear fit
        'r_value': fit[2],              # Linear correlation (Pearson)
        'p_value': fit[3],              # p-value for Linear correlation
        'std_err': fit[4],              # Error of the estimated slope
        'ss': 1 - (rmse**2 / xd[3]),    # Murphy Skill Score
        'count': xd[0]                  # Number of pairs
    }


def bin_xyz(x, y, z, delta=(1., 1.), limit=None, globe=False, order=False):
    """
    Bin (average) irregular 1D data (triplets) on to 2D plane

    *** DEPRECATED *** Use ypylub.utils.XYZ(...).griddata()

    Args:
     * x array_like, shape(N,) An array containing the x coordinates of the
       points to be binned.
     * y array_like, shape(N,) An array containing the y coordinates of the
       points to be binned.
     * z: array_like, shape(N,) f(x,y) actual data to be re-sampled (average
       at each grid cell).
     * delta: float or [float, float], optional Output grid resolution in x
       and y direction.

    The delta specification:
     * If scalar, the grid resolution for the two dimensions (dx=dx=delta).
     * If [float, float], the grid resolution for the two dimensions (dx,
       dy = delta).
     * limit [[float,float],[float,float]], optional Output domain limit
       [lon0,lon1], [lat0,lat1]
     * globe: bool, optional If True, sets the grid x and y limit to
       [-180,180] and [-90,90], respectively. If False, grid x and y
       limits are taken from input.
     * order: bool, optional If True, returns a upside-down flip and rotated
       array

    Returns:
     * G: MaskedArray, shape(nxc,nyc) The bi-dimensional binned (averaged)
       array of z.
     * xc: ndarray, shape(nx,) The bin centres along the x dimension
     * yc: ndarray, shape(ny,) The bin centres along the y dimension.

    Example:
    ::
        >>> from numpy.random import normal
        >>> from ypylib.stats import bin_xyz
        >>> x = normal(3, 1, 100)
        >>> y = normal(1, 1, 100)
        >>> z = x * y
        >>> dd, xc, yc = bin_xyz(x, y, z, delta=[0.1,0.1])

    """
    if np.size(delta) == 1:
        delta = [delta, delta]
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)

    limit = [[[x.min(), x.max()], [y.min(), y.max()]],
             [[-180., 180.], [-90., 90.]]][globe]

    xs = np.arange(limit[0][0], limit[0][1] + delta[0], delta[0])
    ys = np.arange(limit[1][0], limit[1][1] + delta[1], delta[1])

    Hv, xl, yb = np.histogram2d(x, y, weights=z, bins=[xs, ys])  # sum and edge
    xc, yc = xl[:-1] + delta[0] / 2., yb[:-1] + delta[1] / 2.  # centre of bins
    Hn, _, _ = np.histogram2d(x, y, bins=[xs, ys])  # counts
    Hvm = np.ma.masked_where(Hn == 0, Hv)  # mask sum array where count = 0

    if order is True:
        # order the data to column major format for display
        return np.flipud(np.rot90(Hvm / Hn)), xc, yc
    else:
        return Hvm / Hn, xc, yc


def creategrid(x1, x2, y1, y2, dx, dy, mesh=True):
    '''Output grid within geo-bounds and specific cell size.

    Args:
     * x1, x2 (real) values of lower, upper limits of X (longitude)
     * y1, y2 (real) values of lower, upper limits of Y (latitude)
     * dx, dy (real) X (longitude), Y (latitude) grid size
    '''
    x1, x2 = np.floor(x1), np.ceil(x2)
    y1, y2 = np.floor(y1), np.ceil(y2)

    nx = (np.ceil((x2 - x1) / dx)).astype(int)
    ny = (np.ceil((y2 - y1) / dy)).astype(int)

    x_grid = np.zeros(nx)  # fill with lon_min
    y_grid = np.zeros(ny)  # fill with lon_max
    x_grid = x_grid + (np.asarray(range(nx)) * dx)
    y_grid = y_grid + (np.asarray(range(ny)) * dy)

    x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    if mesh is not True:
        x_grid = np.ravel(x_grid)
        y_grid = np.ravel(y_grid)

    return x_grid, y_grid


# def grid_xyz(x, y, z, xo, yo):
#     '''Bin (average) 1d scattered data on a plane

#     Args:
#      * x: array_like, shape(N,) An array containing the x coordinates of the
#        points to be gridded.
#      * y: array_like, shape(N,) An array containing the y coordinates of the
#        points to be gridded.
#      * z: array_like, shape(N,) f(x,y) actual data to be re-sampled (average
#        at each grid cell).
#      * xo: array_like, shape(X,) An array (monotonically increasing)
#        containing the x coordinates of the output grid.
#      * yo: array_like, shape(Y,) An array (monotonically increasing)
#        containing the y coordinates of the output grid.

#     Returns:
#      * D: ndarray, shape(nxo, nyo) The 2D gridded array using sample z

#     See also:
#      * bin_xyz
#     '''
#     x = np.asarray(x, x.dtype)
#     y = np.asarray(y, y.dtype)
#     z = np.asarray(z, z.dtype)

#     # TODO: do some sanity check for 1d arrays with equal number of elements
#     nrow = yo.shape[0]
#     ncol = xo.shape[0]

#     # step size (rectangular cells)
#     dx, dy = xo[1] - xo[0], yo[1] - yo[0]
#     ex, ey = dx / 2.0, dy / 2.0
#     # print xo.min(), xo.max()

#     w = np.where((x >= xo.min()) & (x <= xo.max()) &
#                  (y >= yo.min()) & (y <= yo.max()))

#     # Shrink array
#     x = x[w]
#     y = y[w]
#     z = z[w]

#     # store sum of z values falling in each bin and count n_sample in each
#     #  bin
#     tot = np.empty((nrow, ncol), z.dtype)
#     tot[:, :] = np.NAN
#     cnt = np.ones((nrow, ncol), 'int')
#     # ones in order to avoid division error

#     print('Binning {} samples on {}x{} grid'.format(
#           z.size, xo.size, yo.size))

#     for i in np.arange(z.size):
#         px = v_locate(xo, x[i])
#         py = v_locate(yo, y[i])
#         if np.abs((xo[px] - x[i]) <= ex and np.abs(yo[py] - y[i]) <= ey):
#             tot[py, px] = np.nansum([tot[py, px], z[i]])
#             cnt[py, px] += 1
#     return tot / cnt  # mean


def josephus(n, k):
    '''Josephus circular elimination (eliminate every kth item) from a sample
    of n.

    f(n,k) = (f(n-1,k)+k) mod n), with f(1,k)=0
    f(n,k) = ((f(n-1,k)+ k-1) mod n) + 1, with f(1,k)=1

    Args:
     * n number of samples
     * k number of items
    '''

    if n == 1:
        return 1
    elif n >= 1000:
        r = 0
        i = 1
        while i <= n:
            r = (r + k) % i
            i += 1
        return r + 1
    else:
        # use recursive function for small n
        return ((josephus(n - 1, k) + k - 1) % n) + 1


def josephus_2(n):
    '''Josephus circular elimination (eliminate every 2nd item) from a sample
    of n.

    f(n) = 2(n - 2^log2(n)) + 1

    Args:
     * n number of sample
    '''
    return 2 * (int(n) - 2 ** (int(log(n, 2)))) + 1


def bias(targets, predictions):
    '''mean bias between two series
    '''
    return np.nanmean(predictions - targets)


def rmse(targets, predictions):
    '''root-mean-squared error between two series
    '''
    return np.sqrt(((predictions - targets) ** 2).mean())


def normalise(data):
    '''Normalise original data between 0 and 1 prange

    Args:
     * data array or list of original data
    '''
    normalised = np.asarray(data)
    return (normalised - normalised.min()) / \
        (normalised.max() - normalised.min())


def nrand(n=5, low=1, high=49):
    '''Creates n random integers between low and high

    Args:
     * low, integer, lower limit
     * high, integer, upper limit
     * n, integer, number of random integers to generate

    Returns:
     * n random numbers in [low, high] range
   '''
    return [randint(low, high) for _ in range(0, n)]


# def bindata(x, y, z, xi, yi, ppbin=False, method='median'):
#     '''Bin irregularly spaced data on a regular grid (centre of the bins).
#     Computes the median (default) or mean value within bins defined by
#     regularly spaced xi and yi coordinates (the grid defining the bins).

#     Args:
#     * x, y: ndarray (1D) The independent variables x- and y-axis of the grid.
#     * z: ndarray (1D) The dependent variable in the form z = f(x,y).
#     * xi, yi: ndarray (1D) The coordinates defining the x- and y-axis of the
#       grid.

#     Kwargs:
#      * ppbin: boolean, optional The function returns `bins` variable
#        (see below for description): [False | True].
#      * method: string, optional The statistical operator used to compute the
#        value of each bin: ['median' | 'mean'].

#     Returns:
#      * grid: ndarray (2D) The evenly binned data. The value of each cell is
#        the median (or mean) value of the contents of the bin.
#      * bins: ndarray (2D) A grid the same shape as `grid`, except the value of
#        each cell is the number of points per bin. Returns only if `ppbin` is
#        set to True.

#     Revisions:
#     Implemented from Fernando Paolo's initial version (2010-11-06).
#     '''

#     if x.ndim != y.ndim != z.ndim != 1 \
#             or x.shape[0] != y.shape[0] != z.shape[0]:
#         raise TypeError('input x,y,z must be all 1D arrays of the same length')

#     if method == 'median':
#         median = True
#     else:
#         median = False

#     # make the grid
#     print('Binning {} samples on {}x{} grid'.format(z.size, xi.size, yi.size))
#     nrow = yi.shape[0]
#     ncol = xi.shape[0]
#     grid = np.empty((nrow, ncol), dtype=xi.dtype)
#     if ppbin:
#         bins = np.copy(grid)

#     # step size (rectangular cells)
#     dx = xi[1] - xi[0]
#     dy = yi[1] - yi[0]
#     hx = dx / 2.
#     hy = dy / 2.

#     # bin the data
#     for row in xrange(nrow):
#         for col in xrange(ncol):
#             xc = xi[col]  # (xc,yc) = center of the bin
#             yc = yi[row]
#             ind, = np.where((xc - hx <= x) & (x < xc + hx) &
#                             (yc - hy <= y) & (y < yc + hy))
#             npts = len(ind)
#             if npts > 0:
#                 if median:
#                     grid[row, col] = np.median(z[ind])
#                 else:
#                     grid[row, col] = np.mean(z[ind])
#                 if ppbin:
#                     bins[row, col] = npts
#             else:
#                 grid[row, col] = np.nan
#                 if ppbin:
#                     bins[row, col] = 0

#     # return the grid
#     if ppbin:
#         return grid, bins
#     else:
#         return grid


# def griddata(x, y, z, binsize=1, retbin=True, retloc=True):
#     # taken from:
#     # http://wiki.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
#     """"Place unevenly spaced 2D data on a grid by 2D binning (nearest neighbour
#     interpolation).

#     Args:
#      * x: ndarray (1D) The independent data x-axis of the grid.
#      * y: ndarray (1D) The independent data y-axis of the grid.
#      * z: ndarray (1D) The dependent data in the form z = f(x,y).

#     Kwargs:
#      * binsize: scalar, optional The full width and height of each bin on the
#        grid. If each bin is a cube, then this is the x and y dimension. This
#        is the step in both directions, x and y.
#      * retbin: boolean, optional Function returns `bins` variable (see below
#        for description) if set to True.
#      * retloc: boolean, optional Function returns `wherebins` variable (see
#        below for description) if set to True.

#     Returns:
#      * grid: ndarray (2D) The evenly gridded data.  The value of each cell is
#        the median value of the contents of the bin.
#      * bins: ndarray (2D) A grid the same shape as `grid`, except the value of
#        each cell is the number of points in that bin.  Returns only if `retbin`
#        is set to True.
#      * wherebin: list (2D) A 2D list the same shape as `grid` and `bins` where
#      each cell contains the indicies of `z` which contain the values stored in
#      the particular bin.

#     Revisions:
#     2010-07-11  ccampo  Initial version

#     """
#     # get extrema values.
#     xmin, xmax = x.min(), x.max()
#     ymin, ymax = y.min(), y.max()

#     # make coordinate arrays.
#     xi = np.arange(xmin, xmax + binsize, binsize)
#     yi = np.arange(ymin, ymax + binsize, binsize)
#     xi, yi = np.meshgrid(xi, yi)

#     # make the grid.
#     grid = np.zeros(xi.shape, dtype=x.dtype)
#     nrow, ncol = grid.shape
#     if retbin:
#         bins = np.copy(grid)

#     # create list in same shape as grid to store indices
#     if retloc:
#         wherebin = np.copy(grid)
#         wherebin = wherebin.tolist()

#     # fill in the grid.
#     for row in range(nrow):
#         for col in range(ncol):
#             xc = xi[row, col]  # x coordinate.
#             yc = yi[row, col]  # y coordinate.

#             # find the position that xc and yc correspond to.
#             posx = np.abs(x - xc)
#             posy = np.abs(y - yc)
#             ibin = np.logical_and(posx < binsize / 2., posy < binsize / 2.)
#             ind = np.where(ibin is True)[0]

#             # fill the bin.
#             my_bin = z[ibin]
#             if retloc:
#                 wherebin[row][col] = ind
#             if retbin:
#                 bins[row, col] = my_bin.size
#             if my_bin.size != 0:
#                 binval = np.median(my_bin)
#                 grid[row, col] = binval
#             else:
#                 grid[row, col] = np.nan  # fill empty bins with nans.

#     # return the grid
#     if retbin:
#         if retloc:
#             return grid, bins, wherebin
#         else:
#             return grid, bins
#     else:
#         if retloc:
#             return grid, wherebin
#         else:
#             return grid


# def plotbins(xi, yi, grid, cmap='Spectral_r'):
#     '''Plots data binned with bin_xyz'''

#     if xi.shape[0] < 2 or yi.shape[0] < 2:
#         raise TypeError('x- or y-axis too small: N data < 2')
#     dx = xi[1] - xi[0]
#     dy = yi[1] - yi[0]
#     left, right, bottom, top = xi.min(), xi.max(), yi.min(), yi.max()
    # extent = (left - dx / 2., right + dx / 2.,
    #           bottom - dy / 2., top + dy / 2.)

#     plt.imshow(grid, extent=extent, aspect='auto', origin='lower',
#                cmap=cmap, interpolation='none')
#     print "plotting.."
#     plt.colorbar()  # draw colour bar
#     plt.title('Binned data')
#     plt.show()


def main():
    pass
    # from ypylib.utils import doy
    # from ypylib.geo import ll_vec2arr
    # print doy(day=1, month=3, year=2000)
    # print doy(2013, 10, 10)
    # xv = np.arange(10)
    # yv = np.arange(5)
    # xx, yy = ll_vec2arr(xv, yv)
    # print xx.shape, yy.shape
    # print xx, yy


if __name__ == '__main__':
    main()
