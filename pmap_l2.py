#!/usr/bin/env python2.7
'''
:Module: ypylib.pmap_l2
:File: /net/home/h05/fra6/PyWorkspace/ypy/ypylib/pmap_l2.py
Created on 7 Apr 2016 14:39:14

:author: yaswant.pradhan (fra6)
:copyright: British Crown Copyright 2016, Met Office

'''

import os
import re
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from mpl_toolkits.basemap import Basemap
from ImageMetaTag import savefig
from datetime import datetime as dt
from netCDF4 import Dataset
from ypylib.stats import bin_xyz

# =============================================================================
# Include Lato Font paths for plot titles
# =============================================================================
exfs = r'/usr/lib/rstudio/resources/presentation/revealjs/fonts'
lp = exfs if os.path.exists(exfs) is True else r'/data/users/fra6/Fonts'
Lato = [lp + '/Lato-' + str(i) for i in ['Bold.ttf', 'Regular.ttf']]
LatoB = fm.FontProperties(fname=Lato[0])
LatoR = fm.FontProperties(fname=Lato[1])


class Level2Files(object):
    '''
    classdocs
    '''

    def __init__(self, **kw):
        '''
        Constructor
        '''
        self.files = []
        self.date = kw.get('date', dt.today().strftime('%Y%m%d'))
        self.local = '/scratch/fra6/PMAP/' + self.date
        self.sat = kw.get('sat', 'METOPA')
        self.daytitle = None
        self.daybase = None
        self.imagefile = None
        self.status = 0

    def __enter__(self):
        return self


# ============================================================================
# Read and consolidate data from Level 2 files
# ============================================================================
    def consolidateDailyAod(self, filepath=None, sat=None):
        if filepath is None:
            filepath = self.local
        if sat is None:
            sat = self.sat

        files = glob.glob('/'.join([filepath, 'W_XX-*' + sat + '*.nc']))

        # Initialise consolidated arrays:
        caod, clon, clat = [], [], []

        if len(files) == 0:
            print "No files found in " + filepath
            self.status = -1

        else:
            print dt.utcnow().strftime('%T') + \
                " Consolidating " + self.sat + " files"

            for i, fi in enumerate(files):
                if (i + 1) % 40 == 0:
                    print "."
                else:
                    print ".",

#             for fi in files:

                with Dataset(fi, 'r') as f:
                    try:
                        grp = f.groups['Data'].groups['MeasurementData']
                        aop = grp.groups['ObservationData'].groups['Aerosol']
                        geo = grp.groups['GeoData']
                        aod = aop.variables['aerosol_optical_depth']
                        w = np.where(aod[:] > 0)

                        clon.append(geo.variables['aerosol_center_longitude'][w])
                        clat.append(geo.variables['aerosol_center_latitude'][w])
                        caod.append(aop.variables['aerosol_optical_depth'][w])

                    except KeyError as e:
                        print "Invalid Key:", e

            clon = [val for sublist in clon for val in sublist]
            clat = [val for sublist in clat for val in sublist]
            caod = [val for sublist in caod for val in sublist]
            print "done"

            z = re.split(r"[,\.]", os.path.basename(fi))[2].split('_')
            z[3] = self.date
            self.daytitle = '.'.join([z[i] for i in (0, 1, 2, 3, 5, 6, 7, 8)])
            self.daybase = '.'.join([z[i] for i in (0, 1, 2, 3)])

        return clon, clat, caod

    def plotDailyAod(self, filepath=None, sat='METOPA',
                     rebin=(0.5, 0.5),
                     figsize=[8, 5],
                     pngfile=None):

        self.sat = sat
        lon, lat, aod = self.consolidateDailyAod(filepath=filepath, sat=sat)
        if self.status != 0:
            return -1

        if pngfile is None:
            pngfile = '/'.join([self.local, self.daybase + '.png'])
        self.imagefile = pngfile

        print "\n" + dt.utcnow().strftime('%T') + \
              ' Binning data to ' + str(rebin[0]).strip() + 'deg grid...',
        gaod, glon, glat = bin_xyz(lon, lat, aod, delta=rebin,
                                   globe=True, order=True)
        print 'done'

        cmap = 'Spectral_r'
        print dt.utcnow().strftime('%T') + ' Preparing figure...'
        vdate = dt.strptime(self.date, '%Y%m%d')
        lons, lats = np.meshgrid(glon, glat)

        # create figure, axes instances.
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
        gline = (None, None)  # grid line option
        m = Basemap(projection='cyl', lon_0=0, resolution='c')
        m.drawmapboundary(fill_color='1')
        im1 = m.pcolormesh(lons, lats, gaod, vmin=0, vmax=2., shading='flat',
                           cmap=cmap, latlon=True)

        m.drawparallels(np.arange(-90., 99., 30.), labels=[1, 0, 0, 0],
                        fontsize=9, color='gray', linewidth=0.2, dashes=gline)
        m.drawmeridians(np.arange(-180., 180., 30.), labels=[0, 0, 0, 1],
                        fontsize=9, color='gray', linewidth=0.2, dashes=gline)
        m.drawcoastlines(linewidth=0.5)

        cb = m.colorbar(im1, location='bottom', size='5%', pad='15%',
                        format='%g')
        cb.ax.xaxis.label.set_font_properties(LatoR)
        fieldname = 'AOD_550'
        cb.set_label('{0} []'.format(fieldname), fontsize='small',
                     labelpad=-42)
        cb.ax.tick_params(labelsize=9, which='both', direction='in')

        # add plot title
        ax.set_title(self.daytitle + ' : ' + vdate.strftime('%d/%m/%Y'),
                     fontproperties=LatoB)
        plt.figtext(0.945, 0.13, dt.utcnow().strftime('%FZ%T'),
                    fontsize='xx-small', horizontalalignment='right',
                    family='monospace', verticalalignment='baseline',
                    bbox=dict(facecolor='gray', alpha=0.2))

        pngtag = {'PMAP product': 'Aerosol',
                  'PMAP dataset': fieldname,
                  'PMAP bin resolution': str(rebin[0])}
        savefig(pngfile, img_converter=2, img_tags=pngtag)
        print dt.utcnow().strftime('%T') + ' Saved image ' + pngfile

# ============================================================================
# Exit
# ============================================================================
    def __exit__(self, exc_type, exc_value, traceback):
        for fi in self.files:
            os.unlink(fi)


if __name__ == '__main__':
    # Level2Files().plotDailyAod(
    #     sat='METOPA', filepath='/data/local/fra6/Data/Aerosol/PMAP/20160408')
    pass
