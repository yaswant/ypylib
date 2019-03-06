#!/usr/bin/env python2.7
import os
import glob
import matplotlib
if os.environ.get('DISPLAY') is None:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pyhdf.SD import SD, SDC
from ypylib.geo import wrap_lon
from ypylib.modis_hdf import Level2Files as l2f
from mpl_toolkits.basemap import Basemap
from datetime import datetime as dt
# import pdb
# import sys


def plot_c6_granules(filelist=None, **kw):
    '''pcolormesh is real mess when dealing with discontinuous data
    keywords:
        product
        fieldname
        vmin, vmax
        figsize
        gline
        cmap
    '''
    product = kw.get('product', 'MYD04_L2')
    if filelist is None:
        l2path = kw.get(
            'l2path', os.path.join('$DATADIR/MODIS_NRT_C6/', product,
                                   dt.utcnow().strftime('%Y/%j')))
        filelist = glob.glob(os.path.join(l2path, '*.hdf'))
        # print filelist

    if len(filelist) == 0:
        print "No files found."
        return

    # keywords:
    fieldname = kw.get('fieldname', 'AOD_550_Dark_Target_Deep_Blue_Combined')
    vmin = kw.get('vmin', 0)
    vmax = kw.get('vmax', 2)
    figsize = kw.get('figsize', (8, 5))  # figure size
    gline = kw.get('gline', (None, None))  # grid line option
    cmap = kw.get('cmap', plt.cm.get_cmap('Spectral_r', 20))  # @UndefinedVariable

    filedate = str.split(os.path.basename(filelist[-1]), '.')
    title = '.'.join(filedate[i] for i in (0, 1, 3, 4))
    print title, len(filelist), 'files.'

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])
    m = Basemap(projection='cyl', lon_0=0, llcrnrlon=-180, urcrnrlon=180,
                llcrnrlat=-90, urcrnrlat=90, resolution='c')
    # m.shadedrelief(scale=0.2)

    for mfile in filelist:
        # print os.path.basename(mfile)
        print ".",
        h4 = SD(mfile, SDC.READ)
        lon = h4.select('Longitude')
        lat = h4.select('Latitude')
        field = h4.select(fieldname)
        scl = field.attributes()['scale_factor']

        # pdb.set_trace()
        lon_diff = np.diff(lon[:])

        ww = np.where((lon_diff < -50) | (lon_diff > 50))
        if len(ww[0]) > 0:
            if lat[:].max() > 80:
                # Use pcolor (slow) for high latitude data..
                im = m.pcolor(lon[:], lat[:],
                              np.ma.masked_equal(field[:], -9999) * scl,
                              vmin=vmin, vmax=vmax, cmap=cmap)  # , latlon=True)
            else:
                # split array to get western and eastern parts
                aod_l, lon_l = field[:], lon[:]
                # aod_r, lon_r = field[:], lon[:]
                i = 0
                for r in ww[0]:
                    # western hemisphere; mask right part
                    aod_l[r, ww[1][i]:] = -9999
                    lon_l[r, ww[1][i]:] = -9999
                    # lon_r[r, :ww[1][i]] = -9999
                    # aod_r[r, :ww[1][i]] = -9999
                    i += 1

                # Plot slices: left part up to 180W
                im = m.pcolormesh(lon_l[:], lat[:],
                                  np.ma.masked_equal(aod_l[:], -9999) * scl,
                                  vmin=vmin, vmax=vmax, shading='flat',
                                  cmap=cmap)  # , latlon=True)
                # right part up to 180E
                if lat[:].min() > -60.:
                    im = m.pcolormesh(wrap_lon(lon[:]), lat[:],
                                      np.ma.masked_equal(field[:], -9999) * scl,
                                      vmin=vmin, vmax=vmax, shading='flat',
                                      cmap=cmap)  # , latlon=True)

        else:
            im = m.pcolormesh(lon[:], lat[:],
                              np.ma.masked_equal(field[:], -9999) * scl,
                              vmin=vmin, vmax=vmax, shading='flat',
                              cmap=cmap)  # , latlon=True)
    print 'done'
    m.drawmapboundary(fill_color='1')
    m.drawparallels(np.arange(-90., 99., 30.), labels=[1, 0, 0, 0],
                    fontsize=9, color='gray', linewidth=0.2, dashes=gline)
    m.drawmeridians(np.arange(-180., 180., 30.), labels=[0, 0, 0, 1],
                    fontsize=9, color='gray', linewidth=0.2, dashes=gline)
    m.drawcoastlines(linewidth=0.5)
    ax.set_title(title)  # plot title
    # add colour bar
    try:
        cb = m.colorbar(im, location='bottom', size='5%', pad='15%', format='%g')
        cb.set_label('{0} []'.format(fieldname), fontsize='small', labelpad=-42)
        cb.ax.tick_params(labelsize=9, which='both', direction='in')
    except:
        pass
    plt.show()


def main():
    date = '20160118'  # data date
    prod = 'MOD04_L2'  # product id
    start = '0000'  # start hour
    stop = '2400'  # stop hour
    down = False  # download files?

    # output image file
    pngf = ".".join([prod, date, start + '-' + stop, 'png'])

    l2 = l2f(nrt=False, product=prod, date=date)
    l2.plotDailyAod(valid_time=[start, stop], pngfile=pngf, download=down)
    # l2f(nrt=True, collection=6).plotDailyAod(pngfile='Test.png', download=False)

    # Print current directory where the image is saved by default
    print "in " + os.path.dirname(os.path.realpath(__file__))
    print "L2 file location:" + l2.local

if __name__ == '__main__':
    # plot_c6_granules(fieldname='Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate')
    plot_c6_granules()
