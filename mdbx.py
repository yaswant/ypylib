"""
Created on Wed Oct 11 14:38:17 2017

@author: yaswant.pradhan

"""
import os
import pwd
import numpy as np
from metdb import obs
from csv import writer
from ypylib.utils import XYZ


class Query(object):
    """Interrogate MetDB to extract or plot stored obserbations"""
    # class static variables
    user = os.getenv('USER')
    contact = '{}@metoffice.gov.uk'.format(pwd.getpwnam(user).pw_gecos)
    elements = ['YEAR', 'MNTH', 'DAY', 'HOUR', 'MINT', 'LTTD', 'LNGD']
    start = 'TODAY-1/0000Z'
    stop = 'TODAY-1/0600Z'

    def __init__(self, subtype, elements,
                 start=None, stop=None,
                 area=None, ddict=None, hostname=None,
                 keep=False):
        """
        MetDB query instance. Some default values are picked from static
        variables.

        Parameters
        ----------
        subtype : str
            MetDB subtype name. See
            https://metnet2/content/metdb-subtypes-and-retrieval-details.
        elements : str or list of str
            Elements or observation fields to retrieve (default datetime and
            coordinates from class static variables).
        start, stop : str, optional
            Cut-off start and stop times (default TODAY-1/0000Z,
            TODAY-1/0600Z). The start and stop values accepts datetime in
            %Y%m%d/%H%MZ format.
        area : sequence of str, optional
            Cut-off geographic Lat Lon extent (default None). The values in
            area must follow North, South, West, East order and with N|S W|E
            suffixes. For example, ['30N', '30S', '10W', '60E'] or
            ('30N', '10N', '60E', '100E')
        hostname : str, optional
            The hostname of the MetDB RPC server (default None which points to
            'mdbapus-prod'). For retrievals from test storage areas use
            'mdb-test' with unique ddict value.
        ddict : str, optional
            This can be used to specify a user-supplied data dictionary to be
            used for retrievals. Some specific subtypes and most data on test
            server use the following format
            '\"/var/moods/tests/``subtype``/retrieval_table\"'
        keep : bool, optional
            Setting this ``True`` will attach requested data to self.data for
            further use by plot(), to_csv(), to_dump() methods without sending
            duplicate requests to MetDB (default False).
            Note that this option has no effect on extract() method which
            always sends a new request to MetDB at each call.

        Rrturns
        -------
        ypylib.mdbx.Query : instance

        Examples
        --------
            >>> # Simple Example
            >>> from ypylib import mdbx
            >>> from datetime import datetime
            >>> my_query = mdbx.Query('SATAOD', 'AOD_NM550')

            # Extract data
            >>> my_query.extract()

            # Plot data
            >>> my_query.plot('AOD_NM550')

            # Write csv
            >>> my_query.to_csv('filename.csv')


            >>> # Advanced examples
            >>> # Ex 1. SATAOD (MODIS AOD) from metdb test server
            >>> myq = Query(
                    'SATAOD', 'AOD_NM550',
                    start='20150807/0300Z', stop='20150807/0600Z',
                    ddict='\"/var/moods/tests/SATAOD/retrieval_table\"',
                    hostname='mdb-test')

            >>> # a. get requested records as a numpy ndarray
            >>> data = myq.extract()

            >>> # b. dump pickle of extracted array
            >>> myq.to_dump('sataod')

            >>> # c. load array from dump file
            >>> array = numpy.load('sataod')

            >>> # d. write records to text csv file
            >>> myq.to_csv('sataod.csv')

            >>> # e. show a specific element on map
            >>> myq.plot('AOD_NM550', cb_on=True, globe=True, show=True)

            >>> # f. save map to file
            >>> my_plt = myq.plot('AOD_NM550', cb_on=True, globe=True, vmax=2,
                    figsize=(11, 7), cb_title='AOD_550')
            >>> my_plt.savefig('sataod.png')

            >>> # 2. ASCAT SRFC_SOIL_MSTR from production server for yesterday
                #    between 0z and 6Z (default class values).
            >>> myq = Query('ASCATHR', 'SRFC_SOIL_MSTR')
            >>> myq.plot('SRFC_SOIL_MSTR', cb_on=True, vmin=0, vmax=100,
                         show=True)

            >>> # 3. MODIS AOD from production server for yesterday
            >>> myq = Query('SATAOD', ['AOD_NM550', 'AOD_NM550_ERRR'],
                    stop='TODAY-1/2359Z')
            >>> myq.plot('AOD_NM550_ERRR', show=True, vmin=0, vmax=0.5,
                    cb_on=True, cb_pad='15%')

            >>> # 4. MSG AOD between 20171012/1130Z and 20171012/1230Z
            >>> myq = Query('MSGAOD', 'ARSL_OPTL_DPTH',
                            start='20171012/1130Z', stop='20171012/1230Z')
            >>> myq.plot('ARSL_OPTL_DPTH', show=True, cb_on=True,
                         cb_title='AOD_550')
            >>> myq.plot('ARSL_OPTL_DPTH', show=True, cb_on=True,
                         cb_title='AOD_550', delta=(0.12, 0.12))

            >>> # 5. MSGAOD for SDA WAS
            myq = Query('MSGAOD', 'ARSL_OPTL_DPTH',
                        area=['70N', '0N', '60W', '55E'],
                        start='20171012/1130Z', stop='20171012/1230Z')
            >>> pl = myq.plot(
                    'ARSL_OPTL_DPTH', delta=(0.12, 0.12),
                    limit=[[-60, 55], [0, 70]],
                    gspacing=(20, 20), vmax=2,
                    title='MSG_201710121200_AOD ({})'.format(r'$\tau_{550}$'),
                    cb_on=True, cb_extend='neither',
                    cb_title='dust optical thickness',
                    # figsize=(10, 7),
                    # map_res='h',
                    # show=True,
                    drawcountries=True)
            >>> pl.text(56, 0, r'$\copyright$ Met Office|EUMETSAT  {}'.format(
                    datetime.utcnow().strftime('%F %T')),
                    color='#555555', rotation='vertical', family='Ubuntu',
                    va='bottom')
            >>> pl.savefig('msgaod.png')

        """
        self.subtype = subtype
        self.start = [start, Query.start][start is None]
        self.stop = [stop, Query.stop][stop is None]
        self.keywords = ['START TIME ' + self.start, 'END TIME ' + self.stop]
        if area:
            self.keywords += ['AREA ' + ' '.join(area)]
        if ddict:
            self.keywords += ['DDICT ' + ddict]
        if isinstance(elements, list) is False:
            elements = elements.split()
        elements += Query.elements
        s = set(elements)
        # re-order the elements to put datetime and lat lon first
        self.elements = Query.elements + list(s.difference(Query.elements))
        self.hostname = hostname
        self.ddict = ddict
        self.keep = keep
        self.data = None
        self.mdb_extract = [True, False][self.keep is True]

    def extract(self, verbose=False):
        """
        Extract obs from MetDB

        Parameters
        ----------
        verbose : bool
            increase verbosity

        Returns
        -------
        numpy ndarray matching request

        """
        try:
            data = obs(Query.contact, self.subtype, self.keywords,
                       self.elements, hostname=self.hostname)
            if self.keep:
                self.data = data
            return data
        except IOError as err:
            print "\n", err
            if verbose:
                print '  {}:{}\n  {} - {}\n  {}'.format(
                    self.hostname, self.ddict, self.start, self.stop,
                    sorted(self.elements))

    def to_dump(self, filename):
        """
        Dump a pickle of the extracted array to the specified file.
        The array can be read back with ``pickle.load()`` or ``numpy.load()``.

        Parameters
        ----------
        filename : str
            A string for output dump filename.

        """
        if self.mdb_extract is False and self.data is not None:
            data = self.data
            print('mdbx.Query.to_dump(): Skipping MDB_EXTRACT.')
        else:
            # extract data from metdb, note mdb extract "is" numpy ndarr,
            # so data.dump() should work fine.
            data = self.extract()

        if data is not None:
            data.dump(filename)
            print('mdbx.Query.to_dump(): Saved file {}'.format(filename))

    def to_csv(self, filename):
        """
        Save data to csv text file.

        Parameters
        ----------
        filename : str
            A string for output csv filename.

        """
        def fields_view(arr, fields):
            """
            Create columns views of a structured array.

            Parameters
            ----------
            arr : array
                numpy structured array
            fields : list of str
                list of fields to view in specific order
            """
            dtype2 = np.dtype(
                {name: arr.dtype.fields[name] for name in fields})
            return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

        if self.mdb_extract is False and self.data is not None:
            data = self.data
            print 'mdbx.Query.to_csv(): Skipping MDB_EXTRACT.'
        else:  # extract data from metdb
            data = self.extract()

        if data is None:
            return

        # save data to csv
        with open(filename, 'wb') as fh:
            write = writer(fh)
            write.writerow(self.elements)
            write.writerows(fields_view(data, self.elements))
        print('mdbx.Query.to_csv(): Saved file {}'.format(filename))

    def plot(self, element, show=False, **kw):
        """
        Plot requested data on a cylindrical map.

        Parameters
        ----------
        element : str
            A string specifying the field to show on map. Usually one from the
            request elements.
        show : bool, optional
            Ture displays the plot (default is False which return a pyplot
            object)

        see ypylib.XYZ.mapdata() keyword parameters for adjusting plot style.

        Returns
        -------
        matplotlib.pyplot : object when ``show=False``

        """
        # extract data from metdb
        if self.mdb_extract is False and self.data is not None:
            data = self.data
            print 'mdbx.Query.plot(): Skipping MDB_EXTRACT.'
        else:
            data = self.extract()
        if data is None:
            return

        # show data on map
        host = ['[{}]'.format(self.hostname), ''][self.hostname is None]
        xyz = XYZ(data['LNGD'], data['LTTD'], data[element])

        # update plot title
        title = kw.get('title', '{}:{} {}\n{}-{}'.format(
            self.subtype, element, host, self.start, self.stop))
        kw.update({'title': title})
        plt = xyz.mapdata(**kw)

        if show:
            plt.show()
        else:
            return plt


if __name__ == '__main__':
    pass
