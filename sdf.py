#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Scientific Data Format parser

:author: yaswant.pradhan
:copyright: Crown copyright. Met Office
"""

from __future__ import print_function
from future.utils import iteritems
import os
import numpy as np
import h5py
from ypylib.utils import log
try:
    from pyhdf.HDF import HDF
    from pyhdf.SD import SD, SDC
    from pyhdf.VS import VS  # noqa
    from pyhdf.error import HDF4Error
    h4err = None
except ImportError as h4err:
    pass

# from datetime import datetime
# import pdb

__version__ = "1.0"
__author__ = "Yaswant Pradhan"

# Real Missing Data Indicator
RMDI = -1073741824.0


class h4Parse(object):
    """
    A pyhdf interface to parse hdf4 file.

    Examples
    --------
    >>> d = h4_parse('file.hdf')
    >>> print d.items  # print available datasets in hdf file

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
            # self.sds = sorted(h4.datasets().keys())
            self.sds = sorted(list(h4.datasets()))  # 2 & 3
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
                # print(attrs)
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
        print(fieldnames)
        # print temp[fieldnames].keys()
        # print dir(temp)
        # print temp.keys()
        scaled = dict.fromkeys(list(temp), None)
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

    Last update: June 2017 yaswant.pradhan

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
                print(item)

    def _print_items(self, name, obj):
        """
        Callable function to visititems()

        Note: This is the official approach to walk through h5 datatree and
        get attributes. However, failing on nc4 (h5 model) with current
        version of h5py (1.8.11), so _print_h5_dsets() below was written to
        list h5/nc4 data structure.

        """
        if isinstance(obj, h5py.Group):
            print(name)
            self.items.append(obj.name)
        elif isinstance(obj, h5py.Dataset):
            print(name, obj.shape, obj.dtype)
            self.datasets.append(name)
            self.items.append(obj.name)
        # and attributes
        # for key, val in obj.attrs.iteritems():
        for key, val in iteritems(obj.attrs):
            print("    %s: %s" % (key, val))

    def _print_h5_dsets(self, obj, offset=''):
        """Print data structure of a h5/nc4 file."""

        if isinstance(obj, h5py.File):
            if self.verbose is True:
                print(obj.file, '(File)', obj.name)
        elif isinstance(obj, h5py.Group) or \
            isinstance(obj, h5py.SoftLink) or \
            isinstance(obj, h5py.ExternalLink) or \
                isinstance(obj, h5py.HardLink):
            pass
        elif isinstance(obj, h5py.Dataset):
            # self.items.append(obj.name)
            if self.verbose is True:
                print(obj.name, "\t", obj.shape, obj.dtype)
        else:
            print('WARNING: Unknown item in HDF5 file', obj.name)
            # sys.exit("Execution terminated.")
            raise Exception

        if isinstance(obj, h5py.File) or isinstance(obj, h5py.Group):
            # for _key, val in sorted(dict(obj).iteritems()):
            for _key, val in sorted(iteritems(dict(obj))):
                # if self.verbose is True: print offset,
                # self._print_h5_dsets(val, offset + '')
                try:
                    self._print_h5_dsets(val, offset + '')
                except UserWarning:
                    print("** Skipping: {} **".format(_key))
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
            print("VersionError: lsattr() requires h5py version >=2.3 "
                  "but installed version is", h5py.version.version)

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

    def get_data(self, dsname=None, verbose=False, order=False):
        """
        Get specific datasets from hdf5 file.

        TODO: A better model would be to return [list of] dicts and retain
        original data attributes such as scale_factor, add_offset, etc.

        Parameters
        ----------
        dsname : str or sequence of str, optional
            Full path to the hdf5 datasets in file. If not present get_data()
            returns all valid datasets in the file.
        verbose
            Switch to verbose mode (default False)
        order : bool, optional
            Retain variable sequence as requested in the output dictionary
            (ordered).

        Examples
        --------
        >>> h5 = h5_parse('file.h5')
        >>> data = h5.get_data('/dataset/path')

        # Or in one line
        >>> data = h5_parse('file.h5').get_data('/dataset/path')

        Returns
        -------
        dict or OrderedDict
            (Ordered) Dictionary of requested or all variables from the file.

        """
        from collections import OrderedDict
        odict = OrderedDict() if order else {}
        basestring = str
        if dsname is None:
            dsname = self.get_dslist()
        elif isinstance(dsname, basestring):
            dsname = (dsname,)

        with h5py.File(self.filename, mode='r') as h5f:
            for item in dsname:
                try:
                    odict.update({item: h5f[item][:]})
                except KeyError:
                    log.error('%s: No such variable in %s', item,
                              self.filename)
        return odict

    def get_attr(self, dsname=None):
        attr = {}

        if dsname is None:
            dsname = self.get_dslist()
        elif isinstance(dsname, basestring):
            dsname = (dsname,)

        with h5py.File(self.filename, mode='r') as h5f:
            for item in dsname:
                tmp = {}
                try:
                    vatts = h5f[item].attrs
                    for k in vatts:
                        tmp.update({k: vatts[k]})
                    attr.update({item: tmp})
                except KeyError:
                    log.error('%s: No such variable in %s', item,
                              self.filename)
        return attr

    def __exit__(self, type, value, traceback):
        if self.h5f:
            self.h5f.close()
        return isinstance(value, TypeError)


h4_parse = h4Parse
h5_parse = h5Parse

if __name__ == '__main__':
    pass
