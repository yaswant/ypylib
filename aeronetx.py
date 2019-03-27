#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
AERONET data extractor

"""
import os
import linecache
import tempfile
import requests
from StringIO import StringIO
from datetime import datetime, timedelta
import pandas as pd
from HTMLParser import HTMLParser
__version__ = "2017.01.1"
__author__ = "Yaswant Pradhan"


def date_parse(str_type, fmt="%d:%m:%Y %H:%M:%S"):
    """Return date time parser function.

    Args:
        str_type (str): Date Time string.
        fmt (str, optional): Format specification of str_type.

    Returns:
        TYPE: Description
    """
    return pd.datetime.strptime(str_type, fmt)


class MLStripper(HTMLParser):
    def __init__(self):
        self.reset()
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        return ''.join(self.fed)


def strip_tags(html):
    """Strip HTML tags from string.

    Source: http://stackoverflow.com/questions/753052

    Args:
        html (TYPE): Description

    Returns:
        TYPE: Description
    """
    line = MLStripper()
    line.feed(html)
    return line.get_data()


def parse_v3_site(site, url=None):
    """Validate AERONET v3 site name.

    Compares against https://aeronet.gsfc.nasa.gov/aeronet_locations.txt

    Args:
        site (TYPE): Description
        url (None, optional): Description

    Returns:
        tuple: (boolean, dataframe) - Is site valid, All AERONET sites as a
    pandas dataframe.

    """
    if url is None:
        url = 'https://aeronet.gsfc.nasa.gov/aeronet_locations.txt'
    line = requests.get(url).content
    df = pd.read_csv(StringIO(line.decode('utf-8')), skiprows=1)

    if site in list(df.Site_Name):
        return True, df
    else:
        print "'{}' not found in AERONET Database Site List".format(site)
        return False, df


def show_v3_site():
    """Display AERONET locations on map.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from mpl_toolkits.basemap import Basemap

    loc_df = parse_v3_site(None)
    x = loc_df[1]['Longitude(decimal_degrees)']
    y = loc_df[1]['Latitude(decimal_degrees)']
    z = loc_df[1]['Elevation(meters)']
    # _site = loc_df[1]['Site_Name']

    fig = plt.figure(figsize=(16, 8.2))
    ax = fig.add_subplot(111, title='AERONET Location elevations (m)')
    m = Basemap(resolution='l', projection='cyl', lat_0=0, lon_0=0)
    m.fillcontinents(color='#f2f2f2', lake_color='#46bcec', zorder=0)
    m.drawcoastlines(linewidth=0.5, color='#dcdcdc', zorder=1)
    m.drawcountries(linewidth=0.6, color='#cbcbcb', zorder=2)

    # print z.min(), z.max()
    # norm = mpl.colors.Normalize(vmin=z.min(), vmax=z.max())
    # mc = cm.ScalarMappable(norm=norm, cmap="jet")
    x1, y1 = m(x.values, y.values)
    # m.scatter(x1, y1, c=mc.to_rgba(z), marker="o", zorder=3)
    im = m.scatter(x1, y1, c=z, marker="o", s=15, zorder=3)
    ax.set_xlim([-180, 180])
    ax.set_ylim([-90, 90])

    # colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    plt.colorbar(im, cax=cax)

    # plt.tight_layout()
    plt.show()


def downlaod_v3_site(site, ymd=None, ymd2=None, prd='SDA15', avg='10',
                     hr1=0, hr2=23, verb=False):
    """Download version 3 AERONET data for a specific site.

    Note: Requires curl

    Args:
        site (str): AERONET location name.
        ymd (str, optional): Start date. Defaults to yesterday.
        ymd2 (str, optional): End date. Defaults to yesterday.
        prd (str, optional): AERONET Product code. One from the following:
            AOD10 - Aerosol Optical Depth Level 1.0
            AOD15 - Aerosol Optical Depth Level 1.5
            AOD20 - Aerosol Optical Depth Level 2.0
            SDA10 - SDA Retrieval Level 1.0
            SDA15 - SDA Retrieval Level 1.5 (default)
            SDA20 - SDA Retrieval Level 2.0
            TOT10 - Total Optical Depth based on AOD Level 1.0 (all points)
            TOT15 - Total Optical Depth based on AOD Level 1.5 (all points)
            TOT20 - Total Optical Depth based on AOD Level 2.0 (all points)
        avg (str, optional): AERONET average type:
            10 - All points (default).
            20 - Daily average.
        hr1 (int, optional): Start hour.
        hr2 (int, optional): End hour.
        verb (bool, optional): Verbose mode.

    Returns:
        TYPE: Description

    """
    # Parse AERONET site first
    if parse_v3_site(site)[0] is False:
        return

    host = 'https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3'
    if ymd is None:
        date1 = datetime.utcnow() - timedelta(days=1)
    else:
        date1 = datetime.strptime(ymd, '%Y%m%d')

    y1 = date1.year
    m1 = date1.month
    d1 = date1.day
    if ymd2 is None:
        y2, m2, d2 = y1, m1, d1
        q = "curl -s -X GET '{}?site={}&year={}&month={}&day={}" + \
            "&hour={}&hour2={}&{}=1&AVG={}'"
        cmd = q.format(host, site, y1, m1, d1, hr1, hr2, prd, avg)
    else:
        date2 = datetime.strptime(ymd2, '%Y%m%d')
        y2 = date2.year
        m2 = date2.month
        d2 = date2.day
        q = "curl -s -X GET '{}?site={}&year={}&month={}&day={}&hour={}" + \
            "&year2={}&month2={}&day2={}&hour2={}&{}=1&AVG={}'"
        cmd = q.format(host, site, y1, m1, d1, hr1, y2, m2, d2, hr2, prd, avg)

    if verb:
        print(cmd)
    tmp = next(tempfile._get_candidate_names())
    os.system(cmd + ' > ' + tmp)

    with open(tmp) as f:
        html = f.read()
        os.remove(tmp)
        return strip_tags(html)


def downlaod_v3_region(llyx, uryx, ymd=None, hr1=0, ymd2=None, hr2=23,
                       prd='SDA15', avg='10', strip_html=True, verb=False):
    """Download AERONET data over a rectangular geographic domain.

    Args:
        llyx (sequence): Lower-Left (Latitude, Longitude) limits for region
            extraction.
        uryx (sequence): Upper-Right (Latitude, Longitude) limits for region
            extraction.
        ymd (str, optional): Start date (in YYYYmmdd form) to extract
            (default current date).
        hr1 (int, optional): Start hour to extract (default 0).
        ymd2 (str, optional): End date (in YYYYmmdd form) to extract
            (default current date).
        hr2 (int, optional): End hour to extract (default 23).
        prd (str, optional): Product code to download (default 'SDA15').
            Accepted prd values are:
            AOD10 - Aerosol Optical Depth Level 1.0
            AOD15 - Aerosol Optical Depth Level 1.5
            AOD20 - Aerosol Optical Depth Level 2.0
            SDA10 - SDA Retrieval Level 1.0
            SDA15 - SDA Retrieval Level 1.5 (default)
            SDA20 - SDA Retrieval Level 2.0
            TOT10 - Total Optical Depth based on AOD Level 1.0 (all points)
            TOT15 - Total Optical Depth based on AOD Level 1.5 (all points)
            TOT20 - Total Optical Depth based on AOD Level 2.0 (all points)
        avg (str, optional): Product average indicator (default '10').
            Accepted avg values are:
            10 - All points.
            20 - Daily averages.
        strip_html (bool, optional): Strip HTML tags.
        verb (bool, optional): Verbose mode.

    Returns:
        str: comma separated string of web data that is a very long string
        (use parse_v3_webdata() to convert this data to pandas dataframe).

    """
    host = 'https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3'
    if ymd is None:
        date1 = datetime.utcnow() - timedelta(days=1)
    else:
        date1 = datetime.strptime(ymd, '%Y%m%d')

    no_html = 1 if strip_html else 0
    y1 = date1.year
    m1 = date1.month
    d1 = date1.day
    if ymd2 is None:
        y2, m2, d2 = y1, m1, d1
        q = "curl -s -X GET '{}?lat1={}&lon1={}&lat2={}&lon2={}" +\
            "&year={}&month={}&day={}&hour={}&hour2={}&{}=1&AVG={}" +\
            "&if_no_html={}'"
        cmd = q.format(host, llyx[0], llyx[1], uryx[0], uryx[1],
                       y1, m1, d1, hr1, hr2, prd, avg, no_html)
    else:
        date2 = datetime.strptime(ymd2, '%Y%m%d')
        y2 = date2.year
        m2 = date2.month
        d2 = date2.day
        q = "curl -s -X GET '{}?lat1={}&lon1={}&lat2={}&lon2={}" + \
            "&year={}&month={}&day={}&hour={}" + \
            "&year2={}&month2={}&day2={}&hour2={}&{}=1&AVG={}" +\
            "&if_no_html={}'"
        cmd = q.format(
            host, float(llyx[0]), float(llyx[1]),
            float(uryx[0]), float(uryx[1]),
            y1, m1, d1, hr1, y2, m2, d2, hr2, prd, avg, no_html)

    if verb:
        print cmd

    # download records and return out string with all html tags stripped
    recs = os.popen(cmd).read()
    if strip_html is True:
        # already stripped via if_no_html option
        return recs
    else:
        return strip_tags(recs)


def read_data(filename, version=2):
    """AERONET data reader.

    Read a given AERONET AOT data file, and return it as a pandas dataframe.

    Note: there is a column offset in AERONET Version-3 total AOD files, which
    has been reported to the AERONET web database team, so I wont use any
    hacks to deal with the staggered columns at >=Optical_Air_Mass.

    Args:
        filename (str): AERONET filename.
        version (int, optional): AERONET version. Defaults to 2.

    Returns:
        pandas.DataFrame: DataFrame containing the AERONET data, with the index
        set to the time-stamp of the AERONET observations. Rows or columns
        consisting entirely of missing data are removed. All other columns
        are left as-is.

    Raises:
        ValueError: If version other than 2 or 3.

    """
    # Identify AERONET product name
    file_info = linecache.getlines(filename)[0:7]
    for line in file_info:
        if 'Version' in line:
            prodname = str.strip(line)
            print prodname

    if version == 2:
        skipr = 4
        na = 'N/A'
        renameCol = 'Last_Processing_Date(dd/mm/yyyy)'
        df = pd.read_csv(filename, skiprows=skipr, na_values=[na],
                         parse_dates={'date_time': [0, 1]},
                         date_parser=date_parse)
    elif version == 3:  # read version 3 data
        skipr = 6
        na = -999.0
        renameCol = 'Last_Date_Processed)'
        #
        # read actual header in the Aeronet file
        # add extra column to header so that V3 ragged dataset (ie without
        # headers) can be read correctly as data frame
        #
        hdr = (pd.read_csv(filename, skiprows=skipr, nrows=0)).columns.tolist()
        #
        # update header with dummy wavelength columns
        #
        hdr[-1] = 'w1'
        hdr.extend(['w' + x for x in map(str, range(2, 11))])
        #
        # read values into data frame
        #
        df = pd.read_csv(filename, skiprows=skipr + 1, names=hdr,
                         na_values=[na], parse_dates={'date_time': [0, 1]},
                         date_parser=date_parse)
    else:
        raise ValueError()

    df = df.set_index('date_time')
    # del df['Julian_Day']
    #
    # Drop any rows that are all NaN and any columns that are all NaN and
    # then sort by the index
    #
    an = (df.dropna(axis=1, how='all').dropna(axis=0, how='all').
          rename(columns={renameCol: 'Last_Processing_Date'}).sort_index())
    return an


def parse_v3_web_data(web_data, skip_rows=6):
    """Parse AERONET version 3 web data.

    Args:
        web_data (TYPE): Description
        skip_rows (int, optional): Description

    Returns:
        TYPE: Description
    """
    skipr = skip_rows
    if len(web_data) < 100:
        print '** No records in the web_data'
        return

    hdr = pd.read_csv(StringIO(web_data), skiprows=skipr,
                      nrows=1).columns.tolist()
    hdr[-1] = 'w1'
    hdr.extend(['w' + x for x in map(str, range(2, 13))])

    # read values into data frame
    df = pd.read_csv(StringIO(web_data), skiprows=skipr + 1, names=hdr,
                     na_values=[-999.0], parse_dates={'Date_Time': [1, 2]},
                     date_parser=date_parse)

    df = df.set_index('Date_Time')
    return (df.dropna(axis=1, how='all').dropna(
        axis=0, how='all').sort_index())


def plot_v3_site_sda(site, ymd=None, ymd2=None, hr1=0, hr2=23,
                     prd='SDA15', avg='10', hourly=False, verb=False):
    """Plot AERONET version 3 SDA data.

    Args:
        site (TYPE): Description
        ymd (None, optional): Description
        ymd2 (None, optional): Description
        hr1 (int, optional): Description
        hr2 (int, optional): Description
        prd (str, optional): Description
        avg (str, optional): Description
        hourly (bool, optional): Description
        verb (bool, optional): Description

    Returns:
        TYPE: Description
    """
    import math
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import matplotlib.pyplot as plt
        import matplotlib.dates as dates

    # Parse AERONET site first
    if not parse_v3_site(site)[0]:
        return

    skip_rows = 6
    data = downlaod_v3_site(site, ymd=ymd, ymd2=ymd2, hr1=hr1, hr2=hr2,
                            prd=prd, avg=avg, verb=verb)
    if len(data) < 100:
        print '** No data found for {} on {} **'.format(site, ymd)
        return

    hdr = pd.read_csv(StringIO(data), skiprows=skip_rows,
                      nrows=1).columns.tolist()
    hdr[-1] = 'w1'
    hdr.extend(['w' + x for x in map(str, range(2, 11))])

    # read values into data frame
    df = pd.read_csv(StringIO(data), skiprows=skip_rows + 1, names=hdr,
                     na_values=[-999.0], parse_dates={'date_time': [1, 2]},
                     date_parser=date_parse)
    df = df.set_index('date_time')
    sda = (df.dropna(axis=1, how='all').dropna(
        axis=0, how='all').sort_index())

    # print sda.columns
    lat = sda['Site_Latitude(Degrees)'][0]
    lon = sda['Site_Longitude(Degrees)'][0]
    elvs = 'Alt: {:.0f}m'.format(sda['Site_Elevation(m)'][0])
    instr = 'Id: ' + str(sda['AERONET_Instrument_Number'][0])

    lats = r'{:.3f}$^\circ$N'.format(lat) if lat >= 0 else r'{:.3f}$^\circ$S'.\
        format(math.fabs(lat))
    lons = r'{:.3f}$^\circ$E'.format(lon) if lon >= 0 else r'{:.3f}$^\circ$W'.\
        format(math.fabs(lon))

    # Start SDA plot
    plot_columns = ['Total_AOD_500nm[tau_a]',
                    'Fine_Mode_AOD_500nm[tau_f]',
                    'Coarse_Mode_AOD_500nm[tau_c]']

    if set(plot_columns).issubset(sda.columns):
        pdf = sda[plot_columns]
        del df, sda
    else:
        print "** No SDA data for " + site + ' **'
        return

    # Hourly average?
    pd_version = pd.__version__.split('.')[1]
    if hourly is True:
        if pd_version >= 19:
            dfh = pdf.resample("H", loffset='30Min').mean()
            dfs = pdf.resample("H", loffset='30Min').std()
        else:
            dfh = pdf.resample("H", how='mean', loffset='30Min')
            dfs = pdf.resample("H", how='std', loffset='30Min')
        # print pdfs
        df = dfh
        obs = r'(hourly avg $\pm1\sigma$)'
    else:
        df = pdf
        obs = '(all points)'

    if prd == 'SDA10':
        lev = '1.0'
    elif prd == 'SDA15':
        lev = '1.5'
    elif prd == 'SDA20':
        lev = '2.0'
    else:
        lev = ''

    # -------------------------------------------------- get series statistics
    st = df.describe()
    df.columns = [
        r'Total: {:.3f}; {:.3f}; {:.3f}'.format(st.loc['mean'][0],
                                                st.loc['50%'][0],
                                                st.loc['std'][0]),
        r'Fine: {:.3f}; {:.3f}; {:.3f}'.format(st.loc['mean'][1],
                                               st.loc['50%'][1],
                                               st.loc['std'][1]),
        r'Coarse: {:.3f}; {:.3f}; {:.3f}'.format(st.loc['mean'][2],
                                                 st.loc['50%'][2],
                                                 st.loc['std'][2])]
    if hourly is True:
        dfs.columns = df.columns
    # ------------------------------------------------------------- start plot
    fig, ax = plt.subplots(figsize=(10, 5))  # @UnusedVariable
    # at this point we could use ax = df.plot() for default chart BUT,
    # we want the chart with custom styles for each series, so:
    styles = ['g-', 'bo-', 'rD-']
    sdclr = ['g', 'b', 'r']
    lwd = [1.5, 1.5, 1.5]
    msz = [7, 7, 5]
    mwd = [1.5, 1.5, 0.2]
    for c, st, lw, mw, ms, sdc in zip(
            df.columns, styles, lwd, mwd, msz, sdclr):
        df[c].plot(style=st, lw=lw, ax=ax, ms=ms,
                   markeredgecolor='w', markeredgewidth=mw)

        if hourly is True:
            ax.fill_between(dfs.index, df[c] - dfs[c], df[c] + dfs[c],
                            color=sdc, alpha=0.2)
    # -------------------- format ticks and tick labels based on series length
    ax.set_xticklabels(df.index, rotation=0, ha='center')
    ax.xaxis.set_major_locator(dates.DayLocator())
    ax.xaxis.set_minor_locator(dates.HourLocator(interval=3))
    tsecs = (df.index.max() - df.index.min()).total_seconds()
    if tsecs > 180 * 86400.0:
        ax.xaxis.set_major_locator(dates.YearLocator())
        ax.xaxis.set_minor_locator(dates.MonthLocator())
    elif tsecs > 30 * 86400.0:
        ax.xaxis.set_major_locator(dates.MonthLocator())
        ax.xaxis.set_minor_locator(dates.DayLocator())
    elif tsecs < 3 * 86400.0:
        ax.xaxis.set_major_locator(dates.HourLocator(byhour=range(0, 24, 6)))
        ax.xaxis.set_minor_locator(dates.HourLocator())
    else:
        pass
    ax.xaxis.set_major_formatter(dates.DateFormatter('%H:%M\n%d%b%y'))
    ax.grid(True, which='major', linestyle='-', alpha=0.2)
    ax.grid(True, which='minor', linestyle='-', alpha=0.1)
    # ax.tick_params(axis='both', direction='out')
    # ax.minorticks_on()

    # --------------------------------------------- plot title and axes labels
    ax.set_title(
        'AERONET_V3_L{} SDA {}\n{} ({}, {}, {}, {})'.
        format(lev, obs, site, lats, lons, elvs, instr))
    ax.set_xlabel('')
    ax.set_ylabel(r'AOD_500nm ($\tau_{500}$)')

    plt.tight_layout()
    plt.legend(loc='best', prop={'size': 12})
    # plt.gcf().autofmt_xdate()
    plt.show()


def read_all_sites_times_daily_averages(filename, skip_rows=6,
                                        ymd1=None, ymd2=None,
                                        lat1=-90, lat2=90,
                                        lon1=-180, lon2=180):
    """Read AERONET daily average data.

    Args:
        filename (TYPE): Description
        skip_rows (int, optional): Description
        ymd1 (None, optional): Description
        ymd2 (None, optional): Description
        lat1 (TYPE, optional): Description
        lat2 (int, optional): Description
        lon1 (TYPE, optional): Description
        lon2 (int, optional): Description

    Returns:
        TYPE: Description
    """
    df = pd.read_csv(filename, skiprows=skip_rows, na_values=-999.0)
    return df
    pass


if __name__ == "__main__":
    show_v3_site()

    # d = downlaod_v3_region(
    #    (-5, -35), (35, 5),
    #    ymd='20150807', hr1=0, ymd2='20150807', hr2=23,
    #    prd='SDA15', avg='10')
    # d = parse_v3_web_data(d, skip_rows=5)

    # Data download examples
    # Example 1: download V3 L15 SDA all points for Dakar station
    # x = downlaod_v3_site('Dakar', ymd='20150801',ymd2='20150831')
    # print x

    # Example 2: download V3 L15 SDA all points over a specific geo region
    # web_data = downlaod_v3_region([0, -40], [40, 25],
    #                        ymd='20150812', ymd2='20150812')

    # web_data = downlaod_v3_site('Santa_Cruz_Tenerife', ymd='20150812',
    #                            ymd2='20150812', hr1=8, hr2=9)
    # print parse_v3_web_data(web_data)

    # Example 3: Yesterday's data over India..
    # x = downlaod_v3_region([17.7,72.5],[27.0,90.8])
    # print x

    # Site-specific plotting examples:
    # Example 1: Plot yesterdays observation over Kanpur
    # plot_v3_site_sda('Kanpur')
    # plot_v3_site_sda('Kanpur', ymd='20170303', hr1=5)  # , hourly=True)

    # plot_v3_site_sda('Kanpur', ymd='20170213', ymd2='20170216')
    # plot_v3_site_sda('Bhola', ymd='20170227')

    # Example 2: Plot from a start date to yesterday with hourly averages:
    # plot_v3_site_sda('Kuwait_University', ymd='20170114', hourly=True)

    # plot_v3_site_sda('Capo_Verde', ymd='20150812', ymd2='20150822')

    # locs = ['Praia', 'Calhau', 'Capo_Verde',
    #         'Teide', 'Izana', 'La_Laguna', 'Santa_Cruz_Tenerife']

    # for i in ['Teide', 'Izana', 'La_Laguna', 'Santa_Cruz_Tenerife']:
    #     print i
    #     plot_v3_site_sda(i, ymd='20150812', ymd2='20150813')

    # plot_v3_site_sda('Capo_Verde', ymd='20150820', ymd2='20150821')

    # web_data = downlaod_v3_site('Teide', ymd='20150812', ymd2='20150813',
    #                             verb=True)
    # pdata = parse_v3_web_data(web_data)
    # # print web_data
    # print pdata

    # for i in ['Praia', 'Capo_Verde']:
    #     print i
    #     plot_v3_site_sda(i, ymd='20150820', ymd2='20150821')
