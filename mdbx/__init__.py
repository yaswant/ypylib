#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Interrogate MetDB to extract or plot stored observations.

"""
import os
import pwd
import numpy as np
from metdb import obs, subtypes
from csv import writer
from ypylib.utils import XYZ, ll_parse, log, mkdirp
__version__ = "2019.05.1"
__author__ = "Yaswant Pradhan"


class Query(object):
    """MetDB query instance. Some default values are picked from static
    variables.

    Attributes
    ----------
    area : sequence of str, optional
        Cut-out geographic Lat Lon extent (default None). The values in
        area must follow North, South, West, East order and with N|S W|E
        suffixes. For example, ['30N', '30S', '10W', '60E'] or
        ('30N', '10N', '60E', '100E')
    constrain : dict
        Constrain extraction based on dict key/value. For example, to constrain
        extraction for a Aqua MODIS (sat id 784), use
        `constrain={'STLT_IDNY': 784}`
    contact : str
        User email address.
    data : numpy masked array
        Extracted data array when the query is executed via the `extract()`
        method.
    ddict : str, optional
        This can be used to specify a user-supplied data dictionary to be
        used for retrievals. Some specific subtypes and most data on test
        server use the following format
        '\"/var/moods/tests/``subtype``/retrieval_table\"'
    elements : str or list of str or tuple of tuple
        Elements or observation fields to retrieve. Default datetime and
        coordinates from class static variables are automatically added - these
        are ['YEAR', 'MNTH', 'DAY', 'HOUR', 'MINT', 'LTTD', 'LNGD']
    fmt : list
        Format specification for the default elements
    hostname : str, optional
        The hostname of the MetDB RPC server (default None which points to
        'mdbapus-prod'). For retrievals from test storage areas use
        'mdb-test' with unique ddict value.
    keep : bool, optional
        Setting this ``True`` will attach requested data to self.data for
        further use by plot(), to_csv(), to_dump() methods without sending
        duplicate requests to MetDB (default False).
        Note that this option has no effect on extract() method which
        always sends a new request to MetDB at each call.
    keywords : list
        Contains Start/Stop Time, Area, DDICT, MERGED data specifications.
    limit : list of list
        [[West, East], [North, South]] geographic limit based on area keyword.
    mdb_extract : bool
        Flag to indicate whether a mdb extract is required
    merged : bool, optional
        This is required to retrieve merged data.  Use only with subtypes
        offering this option.
    platform : str, optional
        Specify Platform ID when required (e.g. for buoy data)
    start, stop : str, optional
        Cut-off start and stop times (default TODAY-1/0000Z, TODAY-1/0600Z).
        The start and stop values accepts datetime in `%Y%m%d/%H%MZ` format.
    subtype : str
        MetDB subtype name. See
        https://metnet2/content/metdb-subtypes-and-retrieval-details.
    user : str
        MetDB user id. Default is to read from environment variable.


    Returns
    -------
    ypylib.mdbx.Query : instance


    Examples
    --------
    # 1. Simple Example: instantiate mdbx.Query
    >>> from ypylib import mdbx
    >>> from datetime import datetime
    >>> req = mdbx.Query('SATAOD', 'AOD_NM550')

    # a. Extract data
    >>> req.extract()

    # b. Plot data
    >>> req.plot('AOD_NM550')

    # c. Write csv
    >>> req.to_csv('filename.csv')


    # 2. Advanced examples: interrogate SATAOD (MODIS) from metdb test server
    >>> req = mdbx.Query(
            'SATAOD', 'AOD_NM550',
            start='20150807/0300Z', stop='20150807/0600Z',
            ddict='\"/var/moods/tests/SATAOD/retrieval_table\"',
            hostname='mdb-test')

    # a. get requested records as a numpy ndarray
    >>> data = req.extract()

    # b. dump pickle of extracted array
    >>> req.to_dump('sataod')

    # c. load array from dump file
    >>> array = numpy.load('sataod')

    # d. write records to text csv file
    >>> req.to_csv('sataod.csv')

    # e. show a specific element on map
    >>> req.plot('AOD_NM550', cb_on=True, globe=True, show=True)

    # f. save map to file
    >>> my_plt = req.plot('AOD_NM550', cb_on=True, globe=True, vmax=2,
                          figsize=(11, 7), cb_title='AOD_550')
    >>> my_plt.savefig('sataod.png')

    # 3. ASCAT SRFC_SOIL_MSTR from production server for yesterday
    #    between 0z and 6Z (default class values).
    >>> req = mdbx.Query('ASCATHR', 'SRFC_SOIL_MSTR')
    >>> req.plot('SRFC_SOIL_MSTR', cb_on=True, vmin=0, vmax=100, show=True)

    # 4. MODIS AOD from production server for yesterday
    >>> req = mdbx.Query('SATAOD', ['AOD_NM550', 'AOD_NM550_ERRR'],
                    stop='TODAY-1/2359Z')
    >>> req.plot('AOD_NM550_ERRR', show=True, vmin=0, vmax=0.5,
                 cb_on=True, cb_pad='15%')

    # 5. MSG AOD between 20171012/1130Z and 20171012/1230Z
    >>> req = mdbx.Query('MSGAOD', 'ARSL_OPTL_DPTH',
                    start='20171012/1130Z', stop='20171012/1230Z')
    >>> req.plot('ARSL_OPTL_DPTH', show=True, cb_on=True,
                 cb_title='AOD_550')
    >>> req.plot('ARSL_OPTL_DPTH', show=True, cb_on=True,
                 cb_title='AOD_550', delta=(0.12, 0.12))

    # 6. MSGAOD for SDA WAS
    >>> req = mdbx.Query('MSGAOD', 'ARSL_OPTL_DPTH',
                    area=['70N', '0N', '60W', '55E'],
                    start='20171012/1130Z', stop='20171012/1230Z')
    >>> plt = req.plot(
        'ARSL_OPTL_DPTH', delta=(0.12, 0.12),
        limit=[[-60, 55], [0, 70]],
        gspacing=(20, 20), vmax=2,
        title='MSG_201710121200_AOD ({})'.format(r'$\tau_{550}$'),
        cb_on=True, cb_extend='neither',
        cb_title='dust optical thickness',
        # figsize=(10, 7), map_res='h', show=True,
        drawcountries=True)
    >>> plt.text(56, 0, r'$\copyright$ Met Office|EUMETSAT  {}'.format(
            datetime.utcnow().strftime('%F %T')),
            color='#555555', rotation='vertical', family='Ubuntu',
            va='bottom')
    >>> plt.savefig('msgaod.png')

    """

    # class static variables

    user = os.getenv('USER')
    contact = '{}@metoffice.gov.uk'.format(pwd.getpwnam(user).pw_gecos)
    elements = ['YEAR', 'MNTH', 'DAY', 'HOUR', 'MINT', 'LTTD', 'LNGD']
    fmt = ['%d', '%d', '%d', '%d', '%d', '%0.5f', '%0.5f']
    start = 'TODAY-1/0000Z'
    stop = 'TODAY-1/0600Z'

    def __init__(self, subtype, elements,
                 start=None, stop=None,
                 area=None, ddict=None,
                 platform=None, over=None,
                 merged=False, hostname=None,
                 keep=False, constrain=None):

        self.subtype = subtype
        self.constrain = constrain
        self.start = [start, Query.start][start is None]
        self.stop = [stop, Query.stop][stop is None]
        self.keywords = ['START TIME ' + self.start, 'END TIME ' + self.stop]
        self.area = area
        if area:
            nswe = [ll_parse(i) for i in area]
            self.limit = [[nswe[2], nswe[3]], [nswe[1], nswe[0]]]
            self.keywords += ['AREA ' + ' '.join(area)]
        if ddict:
            self.keywords += ['DDICT ' + ddict]
        if over:
            self.keywords += ['OVER ' + over]
        if merged is True:
            self.keywords += ['DATA MERGED']
        if platform:
            self.keywords += ['PLATFORM ' + platform]

        # update elements list
        if isinstance(elements, str):
            elements = elements.split()
        # for multi-dimensional array
        if isinstance(elements, tuple):
            elements = [elements]
        # add constrained fields to elements list
        if self.constrain:
            elements += list(constrain.keys())
        elements += Query.elements
        s = set(elements)
        # re-order the elements to put datetime and lat lon first
        self.elements = Query.elements + list(s.difference(Query.elements))
        self.hostname = hostname
        self.ddict = ddict
        self.keep = keep
        self.data = None
        self.mdb_extract = False if self.keep is True else True

    def set_element_dtype(self, elements, dtypes):
        """
        Set element data types if not already defined in the subtypes dict.

        Parameters
        ----------
        elements : str or sequence of str
            element names to update
        dtypes : str or sequence of str
            corresponding data types
        """
        if isinstance(elements, str):
            elements = (elements,)
        if isinstance(dtypes, str):
            dtypes = (dtypes,)
        for el, dt in list(zip(elements, dtypes)):
            subtypes.DTYPE_MAPS[self.subtype][el] = dt

    def extract(self, verbose=False, fix_unknown=False):
        """Extract obs from MetDB.

        Parameters
        ----------
        verbose : bool
            increase verbosity
        fix_unknown : bool, optional
            Set data types for unknown elements e.g., merged back model fields
            to 32bit float.

        Returns
        -------
        numpy ndarray matching request

        """
        try:
            data = obs(Query.contact, self.subtype, self.keywords,
                       self.elements, hostname=self.hostname)

        except ValueError as verr:
            if fix_unknown:
                elements = self.elements[len(Query.elements):]
                log.warning("Setting dtypes for %s as 'float'", elements)

                for e in elements:
                    self.set_element_dtype(e, 'float')
                data = obs(Query.contact, self.subtype, self.keywords,
                           self.elements, hostname=self.hostname)
            else:
                log.error('%s\n** Hint for unknown elements: use '
                          'extract(fix_unknow=True) OR call '
                          'set_element_dtype(element, dtype) prior '
                          'to extract() to bypass this error **', verr)
                raise(verr)

        except IOError as err:
            if verbose:
                log.info('  %s:%s\n  %s - %s\n  %s', self.hostname, self.ddict,
                         self.start, self.stop, sorted(self.elements))
            raise(err)

        if self.constrain:
            # single field constrain
            # key = list(self.constrain.keys())[0]
            # value = list(self.constrain.values())[0]
            # data = data[data[key] == value]

            # iterate through all constained fields
            for key, val in self.constrain.items():
                data = data[data[key] == val]
            if len(data) == 0:
                log.error('Constrain results with no data.')
                print(self.constrain)
                return

        if self.keep:
            self.data = data

        return data

    def to_dump(self, filename):
        """Dump a pickle of the extracted array to the specified file.
        The array can be read back with `pickle.load()` or `numpy.load()`.

        Parameters
        ----------
        filename : str
            A string for output dump filename.

        """
        if self.mdb_extract or self.data is None:
            data = self.extract()
        else:
            log.info('Skip MDB_EXTRACT.')
            data = self.data

        if data is not None:
            mkdirp(os.path.dirname(filename))
            data.dump(filename)
            log.info('File saved: %s', filename)

    def to_txt(self, filename, delimiter=',', fmt=None):
        """Save data to a text file.

        Preferred over to_csv() for speed and control.

        Parameters
        ----------
        filename : str
            A string for output text filename.
        delimiter : str, optional
            String or character separating columns.
        fmt : str or sequence of strs, optional
            Format specifier(s) for requested elements. This should exclude
            Date/Time and coordinate fields.

        See also
        --------
        to_csv  Save records to comma-separated text file.

        TODO: Auto-generate filename
        """
        el_header = delimiter.join(self.elements)
        add_n = len(self.elements) - len(Query.elements)
        add_fmt = ['%g'] * add_n if fmt is None else list(fmt)
        el_fmt = Query.fmt + add_fmt

        if self.mdb_extract or self.data is None:
            data = self.extract()
        else:
            log.info('Skip MDB_EXTRACT.')
            data = self.data

        if data is not None:
            mkdirp(os.path.dirname(filename))
            np.savetxt(filename, data,
                       delimiter=delimiter,
                       header=el_header,
                       fmt=el_fmt,
                       comments='')
            log.info('File saved: %s', filename)

    def to_csv(self, filename):
        """Save data to comma separated value (csv) file.

        Parameters
        ----------
        filename : str
            A string for output csv filename.

        See also
        --------
        to_txt  Save records to text file and with more control and better
                performance.

        TODO: Auto generate filename
        """
        def fields_view(arr, fields):
            """Create columns views of a structured array.

            Parameters
            ----------
            arr : array
                numpy structured array
            fields : list of str
                list of fields to view in specific order
            """
            # dtype2 = np.dtype(
            #     {name: arr.dtype.fields[name] for name in fields})
            dtype2 = np.dtype(
                dict((name, arr.dtype.fields[name]) for name in fields))
            return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

        if self.mdb_extract or self.data is None:
            data = self.extract()
        else:
            log.info('Skip MDB_EXTRACT.')
            data = self.data

        if data is None:
            return

        # save data to csv
        with open(filename, 'wb') as fh:
            write = writer(fh)
            write.writerow(self.elements)
            write.writerows(fields_view(data, self.elements))
        log.info('File saved: %s', filename)

    def plot(self, element, index=None, fix_unknown=False, show=False,
             valid_min=None, valid_max=None, **kw):
        """Plot requested data on a cylindrical map.

        Parameters
        ----------
        element : str
            A string specifying the field to show on map. Usually one from the
            request elements.
        index : int
            index of second dimension when retrieving multi-dim array
            eg, if extracted field is (N, d), where N is number of records
            and d is number of wavelengths then index=1 will plot the field
            corresponding to the second wavelength, and so on...
        show : bool, optional
            True displays the plot (default is False which return a pyplot
            object)
        valid_min : real
            Valid minimum physical value in the array. Setting this will mask
            any values < valid_min
        valid_max : real
            Valid maximim physical value in the array. Setting this will mask
            any values > valid_max

        Keywords:
            See ypylib.XYZ.mapdata() keyword parameters for adjusting plot
            style.

        Returns
        -------
        matplotlib.pyplot : object when ``show=False``

        """
        # extract data from metdb
        if self.mdb_extract is False and self.data is not None:
            data = self.data
            log.info('Skip MDB_EXTRACT.')
        else:
            data = self.extract(fix_unknown=fix_unknown)
        if data is None:
            return

        # show data on map
        host = ['[{}]'.format(self.hostname), ''][self.hostname is None]

        if isinstance(element, str) is False:
            log.warning('ELEMENT should be a string, but given\n%s', element)
            return

        if index is None or len(data[element].shape) == 1:
            plt_arr = data[element]
            element_str = element
        else:
            plt_arr = data[element][:, index]
            element_str = '{} [{}]'.format(element, str(index))

        # vmin = kw.get('vmin', plt_arr.min())
        # vmax = kw.get('vmax', plt_arr.max())
        vmin = [valid_min, plt_arr.min()][valid_min is None]
        vmax = [valid_max, plt_arr.max()][valid_max is None]

        if np.ma.is_masked(vmin) and np.ma.is_masked(vmax):
            log.warning('No valid data points to plot.')
            return

        if len(plt_arr.shape) > 1:
            log.warning('Cant decide which array to plot,')
            log.info('Please give me an index value (e.g. index=1).')
            return

        w = (plt_arr >= vmin) & (plt_arr <= vmax)
        if sum(w) == 0:
            log.warning('Number of constrained data points = 0. Abort plot.')
            return

        xyz = XYZ(data['LNGD'][w], data['LTTD'][w], plt_arr[w])

        # w = (data[element] >= vmin) & (data[element] <= vmax)
        # xyz = XYZ(data['LNGD'][w], data['LTTD'][w], data[element][w])

        # update plot title
        title = kw.get('title', '{}:{} {}\n{}-{}'.format(
            self.subtype, element_str, host, self.start, self.stop))
        kw.update({'title': title})

        # print(self.extent)
        if self.area:
            kw.update({'limit': self.limit})

        # warnings.simplefilter("ignore", UserWarning)
        res = xyz.mapdata(**kw)

        if show:
            res.show()
        else:
            return res


if __name__ == '__main__':
    pass
