"""
:author: yaswant.pradhan
:copyright: Crown copyright. Met Office.
"""
import numpy as np
from scipy import stats
from math import log
from random import randint
from scipy.ndimage.filters import uniform_filter1d
# from .utils import v_locate
# from matplotlib import cm
# import matplotlib.pyplot as plt


def ecdf(data):
    """Empirical cumulative distribution function.

    Map every data point in the dataset to a quantile, which is a number
    between 0 and 1 that indicates the cumulative fraction of data points
    smaller than that data point itself.

    Interpretation:
    median: drawn from 0.5 on y-axis


    Parameters
    ----------
    data : array_like
        point dataset

    Returns
    -------
    tuple
        Sorted data, and cumulative fraction of data points
    """
    x, y = np.sort(data), np.arange(1, len(data) + 1) / len(data)
    return x, y


def scale_range(input, min, max):
    """Scale an input array-like to a minimum and maximum number
    the input array must be of a floating point array
    if you have a non-floating point array, convert to floating using
    `astype('float')`
    this works with n-dimensional arrays
    it will mutate in place
    min and max can be integers
    """

    input += -(np.min(input))
    input /= np.max(input) / (max - min)
    input += min
    return input


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    *** Output array is longer by window_len -1 ***

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the beginning and end part of the output signal.

    Parameters
    ----------
    x: the input signal
    window_len: the dimension of the smoothing window;
                should be an odd integer
    window: the type of window from 'flat', 'hanning', 'hamming',
            'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    Output
    ------
        the smoothed signal

    Examples
    --------
    >>> t = linspace(-2,2,0.1)
    >>> x = sin(t)+randn(len(t))*0.1
    >>> y = smooth(x)

    See also
    --------

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve, scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array
          instead of a string
    NOTE: length(output) != length(input), to correct this:
    return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if window not in ('flat', 'hanning', 'hamming', 'bartlett', 'blackman'):
        raise ValueError(
            "Window not 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y


def running_mean(x, N):
    """Calculate running average from a series (x) using a window size (N).

    Parameters
    ----------
    x : array of numbers
    N : scalar integer
    """
    # cumsum = np.cumsum(np.insert(x, 0, 0))
    # return (cumsum[N:] - cumsum[:-N]) / float(N)
    return uniform_filter1d(x, size=N)


def stat2(x, y):
    """Statistical comparison between two series.

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

    Returns
    -------
    dict
        Series stats
        ------------
        x_min, y_min :
            Series minimum value of reference (x), predicted (y)
        x_max, y_max :
            Series maximum value of reference (x), predicted (y)
        x_med, y_med :
            Series median of reference (x), predicted (y)
        x_avg, y_avg :
            Average of reference (x), predicted (y)
        x_var, y_var :
            Sample (n-1) variance of reference (x), predicted (y)
        x_skew, y_skew :
            Skewness of reference (x), predicted (y)
        x_kurt, y_kurt :
            Kurtosis of reference (x), predicted (y)

        Bias vs reference
        -----------------
        bias :
            Mean bias
            mean(y - x)
        med_bias :
            Median bias
            median(y - x)
        rel_bias :
            Relative bias normalised by reference mean
            bias / mean(x)
        nmb :
            Normalised mean bias
            sum(y - x) / sum(x)  or bias / mean(x)
        nmbf :
            Normalised mean bias factor
            bias / mean(x|y) for bias >=0|<0

        Error vs reference
        ------------------
        rmse :
            Root mean square difference
            sqrt(mean( (y - x)**2) )
        rel_rmse :
            Relative rmse normalised by reference mean
            rmse / mean(x)
        nme :
            Normalised mean absolute error
            |sum(y - x)| / sum(x)
        nmef :
            Normalised mean absolute error factor (doi:10.1002/asl.125)
            sum(|y - x|) / sum(x|y) for bias >=0|<0
        urmse :
            Unbiased root mean square difference
            sqrt(rmse**2 - bias**2)
        fge :
            Fractional gross error
            2 * mean(|(y - x) / (y + x)|)

        Goodness of fit
        ---------------
        slope :
            Slope of linear regression fit
        intercept :
            Intercept og linear fit
        r_value :
            Linear correlation (Pearson) coefficient
        p_value :
            p-value for Linear correlation
        std_err :
            Error of the estimated slope
        r_spearman:
            Rank correlation (Spearman) non-parametric
        p_spearman:
            p-value for Rank correlation

        Skill metrics
        -------------
        ss :
            Murphy skill score
            1 - (rmse**2 / var(x))
        d :
            Willmott accuracy index
            1 - ( sum(y - x) / sum( |y - mean(x)| + |x - mean(x)|)**2)

        count :
            Number of pairs

    """

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
    med_bias = np.median(y - x)
    # normalised mean bias factor (nmbf) and normalised mean absolute error
    # factor (nmaef). see doi:10.1002/asl.125
    # result of sum of indiv. factor bias with obs (or model) conc. as a
    # weighting function. this metric avoids undue influence of small numbers
    # in denominator
    # nmbf = [1 - pred.sum/obs.sum, pred.sum/obs.sum - 1][bias >= 0]

    return {
        'x_min': xd[1][0], 'y_min': yd[1][0],
        'x_max': xd[1][1], 'y_max': yd[1][1],
        'x_med': np.median(x), 'y_med': np.median(y),
        'x_avg': xd[2], 'y_avg': yd[2],
        'x_var': xd[3], 'y_var': yd[3],
        'x_skew': xd[4], 'y_skew': yd[4],
        'x_kurt': xd[5], 'y_kurt': yd[5],
        'bias': bias,
        'med_bias': med_bias,
        'rel_bias': bias / xd[2],
        'nmb': np.sum(y - x) / np.sum(x),
        'nmbf': bias / [yd[2], xd[2]][bias >= 0],
        'rmse': rmse,
        'rel_rmse': rmse / xd[2],
        'urmse': np.sqrt(rmse**2 - bias**2),
        'nme': np.abs(np.sum(y - x)) / np.sum(x),
        'nmef': np.sum(np.abs(y - x)) / [np.sum(y), np.sum(x)][bias >= 0],
        'fge': 2.0 * np.mean(np.abs((y - x) / (y + x))),
        'r_spearman': rs_p[0],
        'p_spearman': rs_p[1],
        'slope': fit[0],
        'intercept': fit[1],
        'r_value': fit[2],
        'p_value': fit[3],
        'std_err': fit[4],
        'ss': 1 - (rmse**2 / xd[3]),
        'd': 1 - (np.sum(y - x) /
                  np.sum((np.abs(y - xd[2]) + np.abs(x - xd[2]))**2)),
        'count': xd[0]
    }


def compare(x, y):
    import collections

    x = np.array(x).flatten()
    y = np.array(y).flatten()

    result = collections.namedtuple(
        'Stat2',
        ['bias', 'bias_median', 'bias_norm_mean', 'bias_norm_mean_factor',
         'bias_relative', 'count', 'fractional_gross_error', 'intercept',
         'normalised_mean_error', 'normalised_mean_error_factor',
         'murphy_skill_score', 'p_linear', 'p_rank', 'r_linear', 'r_rank',
         'rmse', 'rmse_relative', 'rmse_unbiased', 'slope', 'standard_error',
         'willmott_index',
         'xmin', 'xmedian', 'xmean', 'xmax', 'xvar', 'xskew', 'xkurt',
         'ymin', 'ymedian', 'ymean', 'ymax', 'yvar', 'yskew', 'ykurt'])

    # remove invalid points from both series
    valid = np.where(np.isfinite(x) & np.isfinite(y))
    x, y = x[valid], y[valid]
    xd = stats.describe(x, axis=None)
    yd = stats.describe(y, axis=None)
    rs_p = stats.spearmanr(x, y)
    fit = stats.linregress(x, y)
    rmse = np.sqrt(np.mean((y - x)**2))
    bias = np.mean(y - x)
    med_bias = np.median(y - x)

    result.bias = bias
    result.bias_median = med_bias
    result.bias_norm_mean = np.sum(y - x) / np.sum(x)
    result.bias_norm_mean_factor = bias / [yd[2], xd[2]][bias >= 0]
    result.bias_relative = bias / xd[2]
    result.count = xd[0]
    result.fractional_gross_error = 2.0 * np.mean(np.abs((y - x) / (y + x)))
    result.intercept = fit[1]
    result.normalised_mean_error = np.abs(np.sum(y - x)) / np.sum(x)
    result.normalised_mean_error_factor = np.sum(np.abs(y - x)) / [np.sum(y), np.sum(x)][bias >= 0]  # noqa
    result.murphy_skill_score = 1 - (rmse**2 / xd[3])
    result.p_linear = fit[3]
    result.p_rank = rs_p[0]
    result.r_linear = fit[2]
    result.r_rank = rs_p[1]
    result.rmse = rmse
    result.rmse_relative = rmse / xd[2]
    result.rmse_unbiased = np.sqrt(rmse**2 - bias**2)
    result.slope = fit[0]
    result.standard_error = fit[4]
    result.willmott_index = 1 - (np.sum(y - x) / np.sum((np.abs(y - xd[2]) + np.abs(x - xd[2]))**2))  # noqa
    result.xmin = xd[1][0]
    result.xmedian = np.median(x)
    result.xmean = xd[2]
    result.xmax = xd[1][1]
    result.xvar = xd[3]
    result.xskew = xd[4]
    result.xkurt = xd[5]
    result.ymin = yd[1][0]
    result.ymedian = np.median(y)
    result.ymean = yd[2]
    result.ymax = yd[1][1]
    result.yvar = yd[3]
    result.yskew = yd[4]
    result.ykurt = yd[5]

    return result


def bin_xyz(x, y, z, delta=(1., 1.), limit=None, globe=False, order=False):
    """Bin (average) irregular 1D data (triplets) on to 2D plane.

    *** DEPRECATED *** Use ypylib.utils.XYZ(...).griddata()

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

    Examples:
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
    """Output grid within geo-bounds and specific cell size.

    Args:
     * x1, x2 (real) values of lower, upper limits of X (longitude)
     * y1, y2 (real) values of lower, upper limits of Y (latitude)
     * dx, dy (real) X (longitude), Y (latitude) grid size
    """
    x1, x2 = np.floor(x1), np.ceil(x2)
    y1, y2 = np.floor(y1), np.ceil(y2)

    nx = (np.ceil((x2 - x1) / dx)).astype(int)
    ny = (np.ceil((y2 - y1) / dy)).astype(int)

    x_grid = np.zeros(nx)  # fill with lon_min
    y_grid = np.zeros(ny)  # fill with lon_max
    x_grid = x_grid + (np.asarray(list(range(nx))) * dx)
    y_grid = y_grid + (np.asarray(list(range(ny))) * dy)

    x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    if mesh is not True:
        x_grid = np.ravel(x_grid)
        y_grid = np.ravel(y_grid)

    return x_grid, y_grid


def josephus(n, k):
    """Josephus circular elimination (eliminate every kth item) from a sample
    of n.

    f(n,k) = (f(n-1,k)+k) mod n), with f(1,k)=0
    f(n,k) = ((f(n-1,k)+ k-1) mod n) + 1, with f(1,k)=1

    Args:
     * n number of samples
     * k number of items
    """

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
    """Josephus circular elimination (eliminate every 2nd item) from a sample
    of n.

    f(n) = 2(n - 2^log2(n)) + 1

    Args:
     * n number of sample
    """
    return 2 * (int(n) - 2 ** (int(log(n, 2)))) + 1


def bias(targets, predictions):
    """Mean bias between two series."""
    return np.nanmean(predictions - targets)


def rmse(targets, predictions):
    """Root-mean-squared error between two series."""
    return np.sqrt(((predictions - targets) ** 2).mean())


def normalise(data):
    """Normalise original data between 0 and 1 prange.

    Args:
     * data array or list of original data
    """
    normalised = np.asarray(data)
    return (normalised - normalised.min()) / \
        (normalised.max() - normalised.min())


def nrand(n=5, low=1, high=49):
    """Create N random integers between low and high.

    Args:
     * low, integer, lower limit
     * high, integer, upper limit
     * n, integer, number of random integers to generate

    Returns:
     * n random numbers in [low, high] range
    """
    return [randint(low, high) for _ in range(0, n)]


def main():
    from ypylib.utils import doy
    import numpy as np
    assert doy(day=1, month=3, year=2000) == 61, "doy assertion failed."
    assert doy(2013, 10, 10) == 283, "doy assertion failed."


if __name__ == '__main__':
    main()
