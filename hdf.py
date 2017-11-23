"""
Created on 29 Jun 2015

:author: yaswant.pradhan
:copyright: Crown copyright. Met Office
"""
import os
import numpy as np
import h5py
from ypylib.geo import MSG
from ypylib import dialog
try:
    from pyhdf.HDF import HDF
    from pyhdf.SD import SD, SDC
    from pyhdf.VS import VS
    from pyhdf.error import HDF4Error
    h4err = None
except ImportError as h4err:
    pass
import matplotlib.pyplot as plt
from netCDF4 import Dataset
# from datetime import datetime
# import pdb
__version__ = '1.0.0'
# Real Missing Data Indicator
RMDI = -1073741824.0


class h4Parse(object):
    """
    A pyhdf interface to parse hdf4 file.

    Examples
    --------
    >>> d = h4_parse('file.hdf')
    >>> print d.items  # print available datasets in hdf file

    Author: yaswant.pradhan
    """

    def __init__(self, filename=None):
        # if hdf4import is False:
        if h4err:
            raise ImportError(
                "{}, which is required to read '{}'".format(
                    h4err, os.path.basename(filename)))
        self.sds = ''
        self.items = []
        self.attr = []
        self.filename = filename
        if filename:
            self._populate_SD()

    def set_filename(self, filename):
        """Set or update hdf filename"""
        self.filename = filename
        self._populate_SD()

    def _populate_SD(self):
        """Populate SDs and their shape attributes"""

        try:
            h4 = SD(self.filename, mode=SDC.READ)
            self.sds = sorted(h4.datasets().keys())
            self.attr.append(h4.attributes())
            for k, v in sorted(h4.datasets().viewitems()):
                self.items.append((k, v[1]))
            h4.end()
        except HDF4Error as e:
            raise HDF4Error('{}: {}'.format(e, self.filename))

    def get_sds(self, fieldnames=[]):
        """
        Returns specific or all SDS in the hdf file as dictionary.

        SDS arrays can be accessed using the 'data' key. Note that no scaling
        is applied to the data in get() method (use get_scaled() to achieve
        that). However, the scaling and missing data information can be
        accessed using the following keys:
            'scale_factor'
            'add_offset'
            '_FillValue'
        """
        # Convert scalar fieldnames to list
        if not isinstance(fieldnames, list):
            fieldnames = [fieldnames]
        # Open file to read SDs
        try:
            h4 = SD(self.filename, mode=SDC.READ)
            sclinfo = None
            if 'Slope_and_Offset_Usage' in h4.attributes():
                sclinfo = 'Slope_and_Offset_Usage'
            # Get all available SDS from file if fieldnames in not given
            if len(fieldnames) == 0:
                fieldnames = []
                for key in sorted(h4.datasets()):
                    fieldnames.append(key)
            # Create and empty dataset dictionary with all available
            # fields fill in data from SDS
            sds = dict.fromkeys(fieldnames, {})
            for key in sds:
                attrs = h4.select(key).attributes()
                if sclinfo:
                    attrs[sclinfo] = h4.attributes()[sclinfo]

                sds[key] = attrs
                sds[key]['data'] = h4.select(key).get()
            # Close hdf interface
            h4.end()
        except HDF4Error as e:
            raise HDF4Error(e)

        # Return raw (possibly un-calibrated) SDS/attributes dictionary
        return sds

    def get_vdata(self, VDataName):
        """
        Return VData (binary table) from hdf4.

        Parameters
        ----------
        VDataName : str
            Name of the VData (stored as binary table in hdf) field

        Returns
        -------
        dict
            returns VData dictionary
        """
        try:
            h4 = HDF(self.filename)
            vs_handle = h4.vstart()
            # in the following vs_handle.vdatainfo() should give information
            # about all vdata, but this does not function correctly with MO
            # installation.
            # print vs_handle.vdatainfo()
            vd = vs_handle.attach(VDataName)
            vdi = vd.fieldinfo()
            vd.detach()
            vdata = {}
            for i in vdi:
                vd = vs_handle.attach(VDataName)
                vd.setfields(i[0])
                vdata[i[0]] = vd.read()
                vd.detach()
            vs_handle.end()
            h4.close()
        except HDF4Error as e:
            raise HDF4Error(e)
        return vdata

    def get_scaled(self, fieldnames=[]):
        """
        Return scaled data asuming that scale_factor and add_offset are
        available in dataset attributes.

        Not a general purpose method, so should be used with caution.
        """
        temp = self.get_sds(fieldnames)
        print fieldnames
        # print temp[fieldnames].keys()
        # print dir(temp)
        # print temp.keys()
        scaled = dict.fromkeys(temp.keys(), None)
        fillvalue = {}
        for k in scaled:
            # see h4.attributes()['Slope_and_Offset_Usage']
            fillvalue[k] = temp[k]['_FillValue']
            scaled[k] = temp[k]['data'] * (
                temp[k]['scale_factor'] - temp[k]['add_offset'])

            w = np.where(temp[k]['data'] == fillvalue[k])
            scaled[k][w] = fillvalue[k]

        # Add FillValues information
        scaled['_FillValues'] = fillvalue

        # Return scaled datasets dictionary
        return scaled


class h5Parse(object):
    """
    Represents structure of a single (simple) HDF5 file.

    Last update: Jun 2017 yaswant.pradhan

    TODO:
        - test with complex HDF5 files.
        - netCDF4 files:
            reading attributes - partially/not functional yet.
    """
    def __init__(self, filename=None, **kw):
        """
        Initialise the HDF5 parser object.

        Parameters
        ----------
        filename : str, optional
            input hdf5 filename (default None). This can be updated using
            .filename attribute.
        verbose : bool, optional
            if true prints all items when instantiated

        Examples
        --------
            >>> h5 = h5Parse("filename.h5")
            # # OR
            >>> h5 = h5Parse()
            >>> h5.filename = "filename.h5"
            >>> h5.items    # quick list items
            >>> h5.ls()     # list groups, datasets (shape attributes)
            >>> data = h5.get_data('/Path/to/Dataset')

            # Image show 2D array (show data values interactively)
            >>> h5.imshow('/Path/to/2D-Dataset', fillvalue=hdf.RMDI,
                          flipx=False, flipy=False, stride=(1, 1))

            # Line plot 1D array
            >>> h5.plot('/Path/to/1D-Dataset')

        See also
        --------
        methods ls(), get_data(), get_dslist(), imshow()
        """
        self.verbose = kw.get('verbose', None)
        self.items = []
        self.groups = []
        self.datasets = []
        self.filename = filename
        self.h5f = None
        if filename:
            self._filetest()
            self._populate_items()

    def __enter__(self):
        return self

    def set_filename(self, filename):
        """Set or update hdf5 filename"""
        self.filename = filename
        self._filetest()
        self._populate_items()

    def _filetest(self):
        """Check filename is a valid (hdf5) file"""
        try:
            open(self.filename)
        except IOError as e:
            raise IOError(e)
        finally:
            if not h5py.is_hdf5(self.filename):
                err = "HDF5Error: not a valid HDF5: '{}'".format(self.filename)
                raise IOError(err)

    def _populate_items(self):
        """Self contained function to populate items in hdf5 file."""

        def list_objects(name, obj):
            if isinstance(obj, h5py.Group):
                self.items.append(name)
                self.groups.append(name)
            if isinstance(obj, h5py.SoftLink) or \
                isinstance(obj, h5py.ExternalLink) or \
                    isinstance(obj, h5py.HardLink):
                self.items.append(name)
            elif isinstance(obj, h5py.Dataset):
                self.items.append(name)
                self.datasets.append(name)
            else:
                pass
        with h5py.File(self.filename, mode='r') as h5f:
            h5f.visititems(list_objects)
        # sort items, groups and datasets
        self.items = sorted(self.items)
        self.groups = sorted(self.groups)
        self.datasets = sorted(self.datasets)
        if self.verbose:
            for item in self.items:
                print item

    def _print_items(self, name, obj):
        """
        Callable function to visititems()

        Note: This is the official approach to walk through h5 datatree and
        get attributes. However, failing on nc4 (h5 model) with current
        version of h5py (1.8.11), so _print_h5_dsets() below was written to
        list h5/nc4 data structure.

        """
        if isinstance(obj, h5py.Group):
            print name
            self.items.append(obj.name)
        elif isinstance(obj, h5py.Dataset):
            print name, obj.shape, obj.dtype
            self.datasets.append(name)
            self.items.append(obj.name)
        # and attributes
        for key, val in obj.attrs.iteritems():
            print "    %s: %s" % (key, val)

    def _print_h5_dsets(self, obj, offset=''):
        """Print data structure of a h5/nc4 file."""

        if isinstance(obj, h5py.File):
            if self.verbose is True:
                print obj.file, '(File)', obj.name
        elif isinstance(obj, h5py.Group) or \
            isinstance(obj, h5py.SoftLink) or \
            isinstance(obj, h5py.ExternalLink) or \
                isinstance(obj, h5py.HardLink):
            pass
        elif isinstance(obj, h5py.Dataset):
            # self.items.append(obj.name)
            if self.verbose is True:
                print obj.name, "\t", obj.shape, obj.dtype
        else:
            print 'WARNING: Unknown item in HDF5 file', obj.name
            # sys.exit("Execution terminated.")
            raise Exception

        if isinstance(obj, h5py.File) or isinstance(obj, h5py.Group):
            for _key, val in sorted(dict(obj).iteritems()):
                # if self.verbose is True: print offset,
                # self._print_h5_dsets(val, offset + '')
                try:
                    self._print_h5_dsets(val, offset + '')
                except UserWarning:
                    print "** Skipping: {} **".format(_key)
                    pass

    def ls(self):
        """
        Recursively list all items, silently ignoring links to external
        files. This should work with h5 files with or without external links.

        Examples
        --------
            >>> h5_parse(h5file).ls() # OR preferably
            >>> h5_parse(h5file).items # for all items

        """
        self.verbose = True
        with h5py.File(self.filename, mode='r') as h5f:
            h5f.visititems(self._print_items)

    def lsd(self):
        """
        Recursively list only datasets in the file. This should work with
        both nc4 and h5 files.

        Examples
        --------
            >>> h5_parse(h5file).lsd()  # OR preferably
            >>> h5_parse(h5file).datasets  # for all datasets
        """
        self.verbose = True
        with h5py.File(self.filename, mode='r') as h5f:
            self._print_h5_dsets(h5f)

    def lsattr(self):
        """
        Equivalent to ls(), but checks for h5py version and cleanly exit
        for h5py<2.3.

        Note: Requires pyhdf >= 2.3

        Examples
        --------
            >>> h5_parse(h5file).lsattr()
        """
        h5py_version = h5py.version.version_tuple
        if h5py_version[0] >= 2 and h5py_version[1] >= 3:
            with h5py.File(self.filename, mode='r') as h5f:
                h5f.visititems(self._print_items)
        else:
            print "VersionError: lsattr() requires h5py version >=2.3 " + \
                "but installed version is", h5py.version.version

    def get_dslist(self):
        """
        Return all valid datasets in hdf5 file as a list.
        Note: This is a redundant method now (h5_parse().dataset stores list
        of all valid datasets) but kept for backward compatibility.

        Examples
        --------
            >>> h5 = h5_parse(h5file)
            >>> dataset_list = h5.get_dslist()  # (old) OR
            >>> dataset_list = h5.datasets  # (prefered)

            # Or in one line
            >>> dataset_list = h5_parse(h5file).get_dslist()
            >>> dataset_list = h5_parse(h5file).datasets
        """
        self.verbose = False
        with h5py.File(self.filename, mode='r') as h5f:
            self._print_h5_dsets(h5f)
        return self.datasets

    def get_data(self, dsname=None, verbose=False):
        """
        Get specific datasets from hdf5 file.

        TODO: A better model would be to return [list of] dicts and retain
        original data attributes such as scale_factor, add_offset, etc.

        Parameters
        ----------
        dsname : str or list of str, optional
            full path to the dataset in h5 file. If not present get_data()
            returns dictionary of all valid datasets in the file (default None)
        verbose, bool, optional
            switch to verbose mode (default False)

        Examples
        --------
            >>> h5 = h5_parse('file.h5')
            >>> data = h5.get_data('/dataset/path')

            # Or in one line
            >>> data = h5_parse('file.h5').get_data('/dataset/path')
        """
        self.h5f = h5py.File(self.filename, mode='r')
        dataset = {}
        if dsname is None:
            if verbose:
                print 'Reading all valid datasets in the file..'
            dlist = self.get_dslist()

            for i in dlist:
                if verbose:
                    print i
                # dataset[i[0:]] = h5f[i].value
                dataset[i] = self.h5f[i].value
                # dataset.append({i[1:]: h5f[i].value})
        else:
            if type(dsname) is str:
                dsname = [dsname]
            for ds in dsname:
                if ds in self.h5f:
                    dataset[ds] = self.h5f[ds].value
                else:
                    print "**Cant find '{}' in {}**".format(ds, self.filename)

            # if dsname in h5f:
            #    dataset = h5f[dsname].value
            # else:
            #    h5f.close()
            #    sys.exit("Cant find {} in {}".format(dsname, self.filename))

        self.h5f.close()
        return dataset

    def imshow(self, dsname, fillvalue=RMDI,
               flipx=False, flipy=False, stride=(1, 1), **kw):
        """
        Display a 2-dimensional array dataset.

        Parameters
        ----------
        dsname : str
            dataset name to show
        fillvalue : float, optional
            Real missing data indicator (default -1073741824.0)
        flipx : bool, optional
            flip array in left-right direction (default False)
        flipy : bool, optional
            flip arrau in up-down direction (default False)
        stride : tuple, optional
            stride number of pixels in x and y direction, reducing the
            resulting image array (default (1, 1))
        cmap : str, optional
            colormap name (default 'gray')

        """
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        flx = 1 if flipx is False else -1
        fly = 1 if flipy is False else -1
        cmap = kw.get('cmap', 'gray')
        interp = kw.get('interpolation', 'none')

        ddict = self.get_data(dsname)
        # dataset = ddict[ddict.keys()[0]]
        key, dataset = ddict.popitem()

        if dataset.ndim != 2:
            print 'Dataset must be a 2D array, but dim of {} is {}.'.\
                format(dsname, dataset.ndim)
            raise Exception

        mdata = np.ma.masked_values(dataset, fillvalue)[::fly, ::flx]
        mdata = mdata[::stride[1], ::stride[0]]
        del dataset

        fig = plt.figure()
        ax = fig.add_subplot(111, title=key)

        m = ax.imshow(mdata, cmap,
                      interpolation=interp)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(m, cax)
        nrows, ncols = mdata.shape

        def format_coord(x, y):
            col = int(x + 0.5)
            row = int(y + 0.5)
            if col >= 0 and col < ncols and row >= 0 and row < nrows and \
                    np.ma.is_masked(mdata[row, col]) is False:
                z = mdata[row, col]
                if flipx is False and flipy is False and stride[0] == 1 and \
                        stride[1] == 1:
                    l, t = MSG().pix2geo(x, y)
                else:
                    l, t = x, y
                return 'x=%1.4f, y=%1.4f, z=%1.4f' % (l, t, z)
            else:
                return 'x=%1.4f, y=%1.4f' % (x, y)
        ax.format_coord = format_coord
        plt.show()

    def plot(self, dsname, fillvalue=RMDI):
        """"Line plot 1D array"""

        dataset = self.get_data(dsname)
        mdata = np.ma.masked_values(dataset, fillvalue)
        plt.plot(mdata)
        plt.title(dsname)
        plt.show()

    def __exit__(self, type, value, traceback):
        if self.h5f:
            self.h5f.close()
        return isinstance(value, TypeError)


# ----------------------------------------------------------------------------
# obsolete functions, but kept for now due to dependecy in other scripts.
def get_sd(FILENAME=None, SDSNAME=None, QUIET=False):
    """
    Return SDS object(s) from HDF4 file

    Parameters
    ----------
    FILENAME : str
        Input statfile
    SDSNAME : list, optional
        List of fields to read in
    QUIET : bool, optional
        Silent all information

    Returns
    -------
    SDS object or list of objects where object.get() method returns a
        numpy ndarray
    """
    # Open file(s)
    # hdf_type = [('HDF SD file', '*.hdf')]
    # fi = dialog.pickfile(filetype=hdf_type) if not FILENAME else FILENAME
    fi = dialog.pickfiles(extension='*.hdf') if not FILENAME else FILENAME
    if (QUIET is False):
        print '-' * (len(fi) + 8)
        print "Opening {0}".format(fi)
        print '-' * (len(fi) + 8)

    # Open HDF SD interface
    try:
        hdf4 = SD(fi, SDC.READ)

        if not SDSNAME:
            # List all valid SDS in the file
            datalist = hdf4.datasets().keys()
            datalist.sort()
            for k in datalist:
                print k
            DATAFIELD_NAME = raw_input(
                "\n>> Please type dataset name(s) to read: ")
            DATAFIELD_NAME = map(str, DATAFIELD_NAME.split())
        else:
            if type(SDSNAME) != list:
                SDSNAME = [SDSNAME]
            DATAFIELD_NAME = SDSNAME

        # Read in SDs
        sds = []
        for idx, val in enumerate(DATAFIELD_NAME):
            if (QUIET is False):
                print "Reading {0} {1} ..".format(idx, val)
            sds.append(hdf4.select(val))
            # print sds.info()

        # Return SDS object instead of list if only one dataset is requested
        if len(sds) == 1:
            sds = sds.pop()
        return sds

    except HDF4Error as e:
        raise HDF4Error(e)


def get_h5(statfile, dataset, start=[0, 0], stop=[None, None], stride=[1, 1]):
    '''
    Read a 2D dataset from HDF5 file

    Args:
     * statfile (str) HDF5 statfile
     * dataset (str) HDF5 dataset to read. The assumption here is that the
             dataset is 2 dimensional

    Kwargs:
     * start (list) [x0, y0] starting position if a sub-region is required
     * stop (list) [x0, y0] ending position if a sub-region is required
     * stride (list) [xn, yn] number of pixels to skip while returning the
             dataset

    Returns:
     * 2D numpy array of specified dataset (or a sliced region see Kwargs)
    '''

    print 'Reading', dataset, '..'
    f = h5py.File(statfile, 'r')
    group = f['/']
    #
    # For compound (structured) dataset the slicing approach wont work as is
    # TODO: think of a clever way to deal with that, eg,
    #      data = group[dataset] and then apply slicing?
    #
    if stop[0] is None or stop[1] is None:
        data = np.ma.array(group[dataset][start[0]::stride[0],
                                          start[1]::stride[1]])
    else:
        data = np.ma.array(group[dataset][start[0]:stop[0]:stride[0],
                                          start[1]:stop[1]:stride[1]])

    f.close()

    return data


def get_nc(statfile, dataset, verb=False):
    if verb:
        print 'Reading', dataset, '..'
    with Dataset(statfile, 'r', format='NETCDF3_CLASSIC') as f:
        data = f.variables[dataset]
        return data[:]


def get_nc4(statfile, dataset):
    print 'Reading', dataset, '..'
    with Dataset(statfile, 'r') as f:
        # Dataset is the class behaviour to open the
        # file and create an instance of the ncCDF4
        # class.

        # data attributes
        # for a in f.variables[k].ncattrs():#
        # print a, f.variables[k].getncattr(a)

        data = f.variables[dataset]
        return data[:]


def plot_sd(SD):
    # TODO:
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as ticker

    name = SD.info()[0]
    rank = SD.info()[1]
    dims = SD.info()[2]

    fill = SD._FillValue
    scale = SD.scale_factor
    offset = SD.add_offset
    # unit = SD.units
    data = SD.get()
    # print data.shape, rank
    scl_data = data * scale + offset
    vrng = np.array(SD.valid_range) * scale + offset

    if rank == 3:
        nc = 3  # number of columns in figure
        nr = dims[0] / nc + 1 if dims[0] % nc > 0 else dims[0] / nc
        fig = plt.figure(figsize=(12, nr * 3))
        gs = gridspec.GridSpec(nr, nc)

        for i in range(0, dims[0]):
            if data[i, :, :].max() == fill:
                continue

            ax = fig.add_subplot(gs[i])
            frame = plt.gca()
            frame.axes.get_xaxis().set_ticks([])
            frame.axes.get_yaxis().set_ticks([])
            ax.set_title('{0}:{1}'.format(name, i + 1), fontsize=10)
            mdata = np.ma.masked_where(data[i, :, :] == fill,
                                       scl_data[i, :, :])
            # print 'Reading array', i
            # print mdata.min(), mdata.max(), mdata.ptp()
            vrng = [mdata.min(), mdata.max()]
            plt.imshow(mdata.T, vmin=vrng[0], vmax=vrng[1],
                       cmap='Spectral_r')

            cb = plt.colorbar(orientation='vertical', fraction=.03,
                              format='%0.2f')
            tick_locator = ticker.MaxNLocator(nbins=6)
            cb.locator = tick_locator
            cb.update_ticks()

    elif rank == 2:
        if data.max() == fill:
            # No valid sds in the selected field, skip further processing..
            raise SystemExit('** W: No valid data in {0} for '
                             'plotting **\n'.format(name))
        else:
            # Create masked array of the selected field filtering out
            # fill_values
            mdata = np.ma.masked_where(data == fill, scl_data)
            fig = plt.figure(figsize=(4, 3))
            ax = fig.add_subplot(111)
            ax.set_title(name)
            plt.imshow(mdata.T, vmin=mdata.min(), vmax=mdata.max(),
                       cmap='Spectral_r')  # @UndefinedVariable
            ax.set_aspect('equal')
            plt.colorbar(orientation='vertical', fraction=.03)

    fig.tight_layout()
    plt.show()


h4_parse = h4Parse
h5_parse = h5Parse

if __name__ == '__main__':

    h5f = '/scratch/fra6/MSG/slotstores/MSG_201706280800.h5'
#     h5Parse(h5f).get_data('MSG/Ch06/BT')

#     h5Parse(h5f).imshow('MSG/Ch09/BT', cmap='Spectral')
    h5Parse(h5f).imshow('MSG/Ch09/BT', flipx=True, flipy=True,
                        stride=(4, 4), cmap='Spectral')
#     h5_parse(h5f).plot('RTTOV/UKV/Ch10/ClearSkyBT')

    # d.imshow('MSG/Ch09/BT', flipx=True, flipy=True)  # , cmap='Spectral_r')
    # Testing
#     h4f = '/data/local/fra6/MODIS_c6c5_comparison/MODIS_SCI_C6/' + \
#           'MYD04_L2/2016/162/MYD04_L2.A2016162.1045.006.2016165154520.hdf'
#     h4f = '/scratch/frmo/Lidar/CALIPSO/' + \
#           'CAL_LID_L15_Exp-Beta-V3-40.2017-06-19T21-00-00ZD.hdf'
#     h5f = '/data/local/fra6/sample_data/slotstores/MSG_201202151600_lite.h5'
#     h5f = '/data/local/fra6/sample_data/ncfile/' + \
#           'W_XX-EUMETSAT-Darmstadt,SOUNDING+SATELLITE,' + \
#           'METOPA+GOME_C_EUMC_20170629032654_55485_eps_o_pmap_l2.nc'
#     h5f = '/scratch/frmo/Lidar/CATS/' + \
#           'CATS_L2O_N-M7.2-V2-TE_05kmNRT.2017-06-13T06-20-31T07-05-53UTC.hdf5'

    # print d.items
    # print d.datasets
    # print len(d.items), len(d.groups), len(d.datasets)
    # d.ls()
