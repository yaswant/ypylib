#!/usr/bin/env python2.7
"""
MODIS Level 2 download, writecsv/netcdf and plotting tool
Created on 05 Feb 2016.

:author: yaswant.pradhan
"""
import os
import socket
import sys
import base64
import csv
import glob
import re
import numpy as np
from contextlib import closing
from datetime import datetime as dt
from datetime import timedelta as td
from shutil import copyfileobj
from urllib2 import urlopen, Request, URLError
from ypylib.hdf import get_sd
from ypylib.stats import bin_xyz
from matplotlib import pyplot as plt
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap
from ImageMetaTag import savefig
if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
    plt.switch_backend('Agg')


# =============================================================================
# Include Lato Font paths for plot titles
# =============================================================================
exfs = r'/usr/lib/rstudio/resources/presentation/revealjs/fonts'
lp = exfs if os.path.exists(exfs) is True else r'$DATADIR/Fonts'
Lato = [lp + '/Lato-' + str(i) for i in ['Bold.ttf', 'Regular.ttf']]
LatoB = fm.FontProperties(fname=Lato[0])
LatoR = fm.FontProperties(fname=Lato[1])


class Level2Files:
    '''
    MODIS aerosol Level-2 hdf read/write/plot module.
    '''

    def __init__(self, collection=6, nrt=False, **kw):
        """
        Kwargs:
        * collection (str) MODIS collection
        * nrt (bool) Option for Near-real-time data
        * date (str) Date in YYYYmmdd form (def: current date)
        * product (str) MODIS level2 product code (def: MYD04_L2 ie aerosol)
        * username (str) username for NRT data access
        * password (str) password for NRT data access
        * local (str) Local path to store level2 files. Only updated via the
          download() method.
       """
        self.files = []
        self.remote = None
        self.local = None
        self.daytitle = None
        self.daybase = None
        self.imagefile = None
        self.collection = collection
        self.status = [0, 0]  # download status, consolidate status
        self.nrt = nrt
        self.date = kw.get('date', dt.today().strftime('%Y%m%d'))
        self.product = kw.get('product', 'MYD04_L2')
        self.username = kw.get('username', None)
        self.password = kw.get('password', None)
        self.c6nrt = r'nrt3.modaps.eosdis.nasa.gov/allData/61'
        self.c6sci = r'ladsftp.nascom.nasa.gov/allData/61'
        self.c5nrt = r'nrt1.modaps.eosdis.nasa.gov/allData/1'
        self.c5sci = r'ladsftp.nascom.nasa.gov/allData/51'
        self.aodfields = kw.get('aodfields', [
            'Longitude', 'Latitude', 'Scan_Start_Time',
            'AOD_550_Dark_Target_Deep_Blue_Combined',
            'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag'])

    def __enter__(self):
        return self

    def dt2yjd(self, date=None):
        """
        Convert date string from YYYYmmdd to YYYY/jjj form.
        This function is kept here because of its relevance to MODIS
        Level2 file structure, but can be used stand-alone.

        Kwargs:
         * date (str) Date in YYYYmmdd. This argument overrides self.date value

        Returns:
         * year/day-of-year (str) in 'YYY/jjj' format.
        """
        if date is None:
            date = self.date
        dtobj = dt.strptime(date, '%Y%m%d')
        return(dtobj.strftime('%Y/%j'))

    def et2dt(self, elapsed_seconds):
        """
        Convert MODIS scan time in elapsed seconds to datetime object

        Returns date and time from MODIS Scan start time which is measured
        in seconds since 1993-1-1 00:00:00.0 0.

        Args:
         * elapsed_seconds (real) Elapsed seconds since 1993-1-1
        """
        base_time = dt(1993, 1, 1, 0, 0, 0)
        return base_time + td(seconds=elapsed_seconds)

    def getCredential(self, profile=None):
        """
        Get MODIS Level-2 NRT data access credential for FTP transfer

        Get eosdis.nasa.gov. user credentials for nrt data access.
        This is read from a text file "$HOME/.eosdis.nasa.gov/profile" where
        the user credentials are stored in "username:password" form as the
        first line.

        Kwargs:
         * profile (str) full path of eosdis user credential file
        """
        if (profile is None):
            profile = os.path.join(os.environ['HOME'],
                                   r'.eosdis.nasa.gov/profile')

        if (self.nrt is True):
            if (self.username is None or self.password is None):
                with open(profile) as f:
                    auth = f.readline().strip().split(':')
                    self.username = base64.b64encode(auth[0])
                    self.password = base64.b64encode(auth[1])

            return [self.username, self.password]

    def getUrl(self):
        """
        Construct URL server path for Level-2 file transfer.

        Get parsed FTP URL for MODIS Level2 nrt or Science product path.
        Note nrt products may not be available on the server for all previous
        dates.
        """
        if (self.nrt is True):
            self.getCredential()

        if (self.collection == 6):
            if (self.nrt is True):
                server = self.c6nrt
                temp = '{0}{1}:{2}@{3}'.format(
                    'ftp://', base64.b64decode(self.username),
                    base64.b64decode(self.password), server)
            else:
                server = self.c6sci
                temp = '{0}{1}'.format('ftp://', server)

        elif (self.collection == 5):
            if (self.nrt is True):
                server = self.c5nrt
                temp = '{0}{1}:{2}@{3}'.format(
                    'ftp://', base64.b64decode(self.username),
                    base64.b64decode(self.password), server)
            else:
                server = self.c5sci
                temp = '{0}{1}'.format('ftp://', server)
        else:
            raise ValueError('Invalid Collection ' + str(self.collection))

        self.remote = '/'.join([server, self.product, self.dt2yjd(self.date)])
        url = '/'.join([temp, self.product, self.dt2yjd(self.date)])
        return url

    def download(self, source=None, destination=None, skip_download=False):
        """
        Download MODIS Level 2 files from NASA server.

        Kwargs:
         * source (str, optional) Full FTP source for MODIS level-2 data to
           override default settings
         * destination (str, optional) Full local path to save the level2 files
        """

        url = self.getUrl() if source is None else source
        if self.nrt is True:
            otag = 'MODIS_NRT_C' + str(self.collection).strip()
        else:
            otag = 'MODIS_SCI_C' + str(self.collection).strip()

        self.local = os.path.join(os.environ['DATADIR']
                                  if destination is None else destination,
                                  otag, self.product, self.dt2yjd())

        if skip_download is True:
            return
        # --------------------------------------------------------- Print info
        print '[Collection={}; NRT={}; Remote={}; Local={}]'.format(
            str(self.collection), str(self.nrt), str(self.remote),
            str(self.local))
        # --------------------------------------------------------- Print info

        # Create destination path if necessary
        if not os.path.exists(self.local):
            try:
                os.makedirs(self.local)
            except IOError as e:
                print e
                sys.exit(1)

        # Begin download l2 files
        cnt = 0
        try:
            # print "URL>>>" + url
            r = urlopen(Request(url))
            regex = r'(' + self.product + '.+?.hdf)'
            # Read l2 filenames from manifest file
            for line in r.readlines():
                if '.hdf' in line:
                    # filter file with matching regex pattern
                    matches = re.findall(regex, line)
                    if (len(matches) > 0):
                        remote = url + '/' + matches[-1]
                        local = self.local + '/' + matches[-1]
                        # Download new files
                        if not os.path.exists(local):
                            cnt += 1
                            # print ' Downloading ' + matches[-1]
                            print matches[-1]
                            with closing(urlopen(remote)) as r:
                                with open(local, 'wb') as f:
                                    copyfileobj(r, f)

            if (cnt > 0):
                print dt.utcnow().strftime('%T') + ' Download complete.\n'
            else:
                print 'All remote files already in: ' + self.local

        except KeyboardInterrupt:
            print 'Interrupted'
            sys.exit(0)
        except URLError, e:
            print ' ** URLError=' + str(e.reason)
            print ' Hint: switching remote server might help'
            if not os.listdir(self.local):
                # remove local directory if empty
                os.rmdir(self.local)
            self.status[0] = -1
            return(self.status[0])

    def consolidateDailyAod(self, filepath=None, valid_time=['0000', '2400'],
                            skip_download=False):
        """
        Concatenate level-2 swath files into single vector array.

        Concatenate AOD data from available 5 minute swath files in to
        single array.

        Kwargs:
         * filepath (str) Location of Level-2 5minute swath files (hdf)

        Returns:
         * consolidated arrays (1D) of: scan_time, longitude, latitude, aod550,
           quality_indicator
        """

        # Download files
        if skip_download is False:
            print dt.utcnow().strftime('%T') + ' Downloading ' + \
                ((self.nrt is True) and 'NRT' or 'SCIENCE') + \
                ' product: ' + self.product

        # assign destination if given
        self.download(destination=filepath, skip_download=skip_download)
        if self.status[0] != 0:
            self.status[1] = -1
            return None, None, None, None, None

        # Search all hdf files
        allfiles = glob.glob('/'.join([self.local, '*.hdf']))
        hdfiles = []
        for k in allfiles:
            if (int(k.split('.')[2]) >= int(valid_time[0]) and
                    int(k.split('.')[2]) <= int(valid_time[1])):
                hdfiles.append(k)

        # hdfiles = glob.glob('/'.join([self.local, '*.hdf']))

        # Initialise consolidated arrays:
        caod, clon, clat, ctime, cqf = [], [], [], [], []

        if (len(hdfiles) > 0):
            print dt.utcnow().strftime('%T') + \
                ' Consolidating swath files...',

            if self.collection == 6:
                for hdfi in hdfiles:
                    try:
                        lon, lat, time, aodc, qf = get_sd(
                            FILENAME=hdfi, SDSNAME=self.aodfields, QUIET=True)
                        aod = aodc.get()
                        fill = aodc.getfillvalue()
                        scale = aodc.scale_factor
                        offset = aodc.add_offset
                        w = np.where(aod > fill)
                        clon.append(lon.get()[w])
                        clat.append(lat.get()[w])
                        ctime.append(time.get()[w])
                        cqf.append(qf.get()[w])
                        caod.append(scale * aod[w] + offset)
                    except:
                        os.remove(hdfi)

            elif self.collection == 5:
                if len(self.aodfields) < 6:
                    self.aodfields[3:5] = [
                        'Optical_Depth_Land_And_Ocean',
                        'Deep_Blue_Aerosol_Optical_Depth_550_Land',
                        'Quality_Assurance_Land']

                for hdfi in hdfiles:
                    lon, lat, time, aodc, aodd, qfl = get_sd(
                        FILENAME=hdfi, SDSNAME=self.aodfields, QUIET=True)

                    qfl5 = qfl.get()[:, :, 4]

                    # Combine DeepBlue and LandOcean AODs
                    aod = aodc.get()
                    fill = aodc.getfillvalue()
                    t = np.where(aod == fill)
                    if len(t) > 0:
                        aod[t] = aodd.get()[t]

                    scale = aodc.scale_factor
                    offset = aodc.add_offset
                    w = np.where(aod > fill)
                    clon.append(lon.get()[w])
                    clat.append(lat.get()[w])
                    ctime.append(time.get()[w])
                    cqf.append(qfl5[w])
                    caod.append(scale * aod[w] + offset)

            else:
                raise ('Invalid Collection ' + str(self.collection))
                sys.exit(0)

        # Flatten list elements
        ctime = [val for sublist in ctime for val in sublist]
        clon = [val for sublist in clon for val in sublist]
        clat = [val for sublist in clat for val in sublist]
        caod = [val for sublist in caod for val in sublist]
        cqf = [val for sublist in cqf for val in sublist]

        z = os.path.basename(hdfi).split('.')
        z[2] = 'daily'
        if self.nrt is not True:
            z[4] = 'SCI'
        self.daytitle = '.'.join([z[i] for i in (0, 1, 3, 4)])
        self.daybase = '.'.join([z[i] for i in (0, 1, 2, 3, 4)])
        print ' done.'
        return ctime, clon, clat, caod, cqf

    def writecsvAod(self, hdfile, csvfile=None, deletehdf=False):
        """
        Read AOD dataset from MODIS c6 hdf and write output to a csv file.

        (lon,lat,time,aod,quality_flag).

        Args:
         * hdfile (str) Level 2 MODIS-C6 AOD hdf filename

        Kwargs:
         * csvfile (str) Output csv filename
         * deletehdf (bool) Remove hdf file after csv write
        """

        if csvfile is None:
                csvfile = os.path.splitext(hdfile)[0] + '.AOD550.csv'

        lon, lat, time, aod, qf = get_sd(FILENAME=hdfile,
                                         SDSNAME=self.aodfields)

        # Get AOD scaling parameters
        scale = aod.scale_factor
        offset = aod.add_offset

        # Get valid sample index
        w = np.where(aod.get() > -9999)

        vlon = lon.get()[w]
        vlat = lat.get()[w]
        vtime = time.get()[w]
        vqf = qf.get()[w]

        # Apply linear scaling to AOD before writing.  "Linear scaling" is
        # assumed here, check the original definition to see if this is
        # still the case
        vaod = scale * aod.get()[w] + offset

        if (vaod.size > 0):
            if (deletehdf is True):
                os.remove(hdfile)

            # Write fields to csv file

            with open(csvfile, 'wb') as fi:
                print ' Writing ' + csvfile
                writer = csv.writer(fi)
                writer.writerow(self.aodfields)
                for i, _ in enumerate(vaod):
                    writer.writerow([vlon[i], vlat[i], self.et2dt(vtime[i]),
                                     vaod[i], vqf[i]])
        else:
            print ' ** Not enough valid AOD retrievals **'
            if (deletehdf is True):
                os.remove(hdfile)

    def writencDailyAod(self, ncfile=None, filepath=None,
                        download=True, **kw):
        """
        Write daily consolidated MODIS AOD data to netCDF file.

        Kwargs:
         * ncfile (str) output netcdf filename
         * filepath (str) override default level 2 file path
         * format (str) netcdf format, default: NETCDF4
        """
        from netCDF4 import Dataset
        ncformat = kw.get('format', 'NETCDF4')
        skip = not download

        time, lon, lat, aod, qf = self.consolidateDailyAod(
            filepath=filepath, skip_download=skip)

        if self.status[1] != 0:
            return -1
        if ncfile is None:
            ncfile = os.path.join(self.local, self.daybase + '.nc')

        print dt.utcnow().strftime('%T') + ' Writing ' + ncfile
        nc = Dataset(ncfile, 'w', clobber=True, format=ncformat)

        # Add Global attributes
        nc.title = self.product
        nc.description = 'MODIS daily aerosol fields accumulated ' + \
            'from {} level 2 swath files'.format(self.product)
        nc.origin = str(self.remote)
        nc.institution = 'Satellite Applications, Met Office, UK.'
        nc.contact = 'yaswant.pradhan@metoffice.gov.uk'
        nc.history = '2016: version 0.1'

        # Define variable dimensions, create variable and add attributes
        # nc.createDimension('x', None) # create sample dimension unlimited
        # Note that an unlimited sample dimension is not efficient at all
        # - creates huge files and take 10x longer for writing
        nc.createDimension('t', np.size(time))
        times = nc.createVariable('scan_time', 'd', 't', zlib=True)
        times.long_name = 'TAI Time at Start of Scan replicated ' + \
            'across the swath'
        times.units = 'Seconds since 1993-1-1 00:00:00.0 0'
        lons = nc.createVariable(
            'lon', 'f', 't', zlib=True, least_significant_digit=3)
        lons.long_name = 'Geodetic Longitude'
        lons.standard_name = 'longitude'
        lons.units = 'degrees_east'
        lats = nc.createVariable(
            'lat', 'f', 't', zlib=True, least_significant_digit=3)
        lats.long_name = 'Geodetic Latitude'
        lats.standard_name = 'latitude'
        lats.units = 'degrees_north'
        aods = nc.createVariable(
            'aod550', 'f', 't', zlib=True, least_significant_digit=4)
        aods.long_name = 'Combined Dark Target, Deep Blue AOT at 0.55 ' + \
            'micron for land and ocean'
        aods.standard_name = 'atmosphere_optical_thickness_due_to_aerosol'
        aods.units = '1'
        aods.valid_range = (-0.1, 5)
        qfs = nc.createVariable(
            'quality_flag', 'h', 't', zlib=True, least_significant_digit=1)
        qfs.long_name = 'Combined Dark Target, Deep Blue Aerosol ' + \
            'Confidence Flag'
        qfs.flag_values = (0, 1, 2, 3)
        qfs.flag_meanings = [
            'No Confidence (or fill), Marginal, ' + 'Good, Very Good']
        qfs.units = '1'
        qfs.valid_range = (0, 3)
        if self.collection == 5:
            qfs.comment = 'C5 Quality_Assurance_Land are encoded 8-bit ' + \
                'data. Decode and combine 1st and 2nd bits to ' + \
                'get AOD retrieval confdence. See http://www-cf/' + \
                '~cfsa/SPS/build/dev/doc/consolidate_mydl2.html'

        # Fill variables
        times[:] = time
        lons[:] = lon
        lats[:] = lat
        aods[:] = aod
        qfs[:] = qf
        nc.close()
        print dt.utcnow().strftime('%T') + ' Done.'

    def writeh5DailyAod(self, h5file=None, filepath=None, download=True):
        """
        Write daily consolidated MODIS AOD data to HDF5 file.
        """
        import h5py
        skip = not download
        time, lon, lat, aod, qf = self.consolidateDailyAod(
            filepath=filepath, skip_download=skip)

        if self.status[1] != 0:
            return -1
        shp = np.shape(time)  # shape

        # print self.local
        # print self.daytitle
        if h5file is None:
            h5file = os.path.join(self.local, self.daybase + '.h5')

        print dt.utcnow().strftime('%T') + ' Writing ' + h5file
        fid = h5py.File(h5file, 'w')

        # Create group and add group attributes
        grp = fid.create_group(self.product)
        grp.attrs['title'] = self.product
        grp.attrs['description'] = 'MODIS Level 2 daily consolidated ' + \
            'aerosol fields from %s' % (self.product)
        grp.attrs['origin'] = str(self.remote)
        grp.attrs['institution'] = 'Satellite Applications, Met Office'
        grp.attrs['contact'] = 'yaswant.pradhan@metoffice.gov.uk'
        grp.attrs['history'] = '2016: version 0.1'

        # Create datasets and add dataset attributes
        times = grp.create_dataset(
            'scan_time', shape=shp, dtype='d', compression='gzip')
        times.attrs['long_name'] = 'TAI Time at Start of Scan replicated ' + \
            'across the swath'
        times.attrs['units'] = 'Seconds since 1993-1-1 00:00:00.0 0'
        lons = grp.create_dataset(
            'lon', shape=shp, dtype='f', compression='gzip')
        lons.attrs['long_name'] = 'Geodetic Longitude'
        lons.attrs['standard_name'] = 'longitude'
        lons.attrs['units'] = 'degrees_east'
        lats = grp.create_dataset(
            'lat', shape=shp, dtype='f', compression='gzip')
        lats.attrs['long_name'] = 'Geodetic Latitude'
        lats.attrs['standard_name'] = 'latitude'
        lats.attrs['units'] = 'degrees_north'
        aods = grp.create_dataset(
            'aod550', shape=shp, dtype='f', compression='gzip')
        aods.attrs['long_name'] = 'Combined Dark Target, Deep Blue AOT ' + \
            'at 0.55 micron for land and ocean'
        aods.attrs['standard_name'] = \
            'atmosphere_optical_thickness_due_to_aerosol'
        aods.attrs['units'] = '1'
        aods.attrs['valid_range'] = [-0.1, 5]
        qfs = grp.create_dataset(
            'quality_flag', shape=shp, dtype='H', compression='gzip')
        qfs.attrs['long_name'] = 'Combined Dark Target, Deep Blue Aerosol ' + \
            'Confidence Flag'
        qfs.attrs['flag_values'] = [0, 1, 2, 3]
        qfs.attrs['flag_meanings'] = [
            'No Confidence (or fill), Marginal, ' + 'Good, Very Good']
        qfs.attrs['units'] = '1'
        qfs.attrs['valid_range'] = [0, 3]
        if self.collection == 5:
            qfs.attrs['comment'] = 'C5 Quality_Assurance_Land are ' + \
                'encoded 8-bit data. Decode and combine 1st ' + \
                'and 2nd bits to get AOD retrieval ' + \
                'confdence. See http://www-cf/~cfsa/' + \
                'SPS/build/dev/doc/consolidate_mydl2.html'

        # Fill datasets
        times[:] = time
        lons[:] = lon
        lats[:] = lat
        aods[:] = aod
        qfs[:] = qf
        fid.close()
        print dt.utcnow().strftime('%T') + ' Done.'

    # -------------------------------------------------------------------------
    # Plot daily consolidated AOD
    # -------------------------------------------------------------------------
    def plotDailyAod(self, rebin=(0.5, 0.5), figsize=(8, 5),
                     valid_time=['0000', '2400'],
                     fieldname='AOD_550_Dark_Target_Deep_Blue_Combined',
                     filepath=None,
                     pngfile=None,
                     download=True):

        skip = not download

        # Consolidate files
        _time, lon, lat, aod, _qf = self.consolidateDailyAod(
            filepath=filepath, valid_time=valid_time, skip_download=skip)

        if self.status[1] != 0:
            return -1

        # Construct png filename
        if pngfile is None:
            pngfile = '/'.join([self.local, self.daybase + '.png'])
        self.imagefile = pngfile

        vdate = dt.strptime(self.date, '%Y%m%d')

        cmap = 'Spectral_r'
        pngtag = {'MODIS product': self.product,
                  'MODIS dataset': fieldname,
                  'MODIS valid date': vdate.strftime('%F'),
                  'MODIS bin resolution': str(rebin[0])}

        # Bin data to regular grid
        print dt.utcnow().strftime('%T') + \
            ' Binning data to ' + str(rebin[0]).strip() + 'deg grid...',
        gaod, glon, glat = bin_xyz(
            lon, lat, aod, delta=rebin, globe=True, order=True)
        print 'done.'

        print dt.utcnow().strftime('%T') + ' Preparing figure...'
        lons, lats = np.meshgrid(glon, glat)

        # create figure, axes instances.
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        gline = (None, None)  # grid line option

        # create Basemap instance : set resolution to None to skip
        # continent processing (this speeds things up a bit)
        m = Basemap(projection='cyl', lon_0=0, resolution='c')
        m.drawmapboundary(fill_color='1')

        # plot data
        im1 = m.pcolormesh(lons, lats, gaod, vmin=0, vmax=2., shading='flat',
                           cmap=cmap, latlon=True)

        # Alternatively plot data as filled contour
        # levels = np.arange(0, 2.01, 0.01)
        # im1 = m.contourf(lons, lats, gaod, levels, cmap=cmap)

        m.drawparallels(np.arange(-90., 99., 30.), labels=[1, 0, 0, 0],
                        fontsize=9, color='gray', linewidth=0.2, dashes=gline)
        m.drawmeridians(np.arange(-180., 180., 30.), labels=[0, 0, 0, 1],
                        fontsize=9, color='gray', linewidth=0.2, dashes=gline)
        m.drawcoastlines(linewidth=0.5)

        # add colour bar
        cb = m.colorbar(im1, location='bottom', size='5%', pad='15%',
                        format='%g')
        cb.ax.xaxis.label.set_font_properties(LatoR)
        cb.set_label('{0} []'.format(fieldname), fontsize='small',
                     labelpad=-42)
        cb.ax.tick_params(labelsize=9, which='both', direction='in')

        # add plot title
        ax.set_title(self.daytitle + ' : ' + vdate.strftime('%d/%m/%Y'),
                     fontproperties=LatoB)

        # add timestamp
        # plt.text(178, -110, dt.utcnow().strftime('%F %T'),
        #         fontsize='xx-small', horizontalalignment='right',
        #         verticalalignment='center',
        #        bbox=dict(facecolor='gray', alpha=0.3))
        plt.figtext(0.945, 0.13, dt.utcnow().strftime('%FZ%T'),
                    fontsize='xx-small', horizontalalignment='right',
                    family='monospace', verticalalignment='baseline',
                    bbox=dict(facecolor='gray', alpha=0.2))
        # print 'done'
        # plt.show()
        # print dt.utcnow().strftime('%T') + ' Saving image ' + pngfile
        # Save figure as pngfile; use imt (a fig wrapper to compress image)
        # fig.savefig(pngfile)
        savefig(pngfile, img_converter=2, img_tags=pngtag)
        print dt.utcnow().strftime('%T') + ' Saved image ' + pngfile
        # print dt.utcnow().strftime('%T') + ' done.'

    # -------------------------------------------------------------------------
    #    Exit
    # -------------------------------------------------------------------------
    def __exit__(self, exc_type, exc_value, traceback):
        for fi in self.files:
            os.unlink(fi)


if __name__ == '__main__':
    pass
