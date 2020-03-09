"""Demonstration of direct plotting of from MetDB using mdbx"""
import os
from ypylib import mdbx


def plot_sataod(ELEMENT='AOD_NM550', AREA=None, START=None, STOP=None,
                constrain=None, **kw):
    """Plot MODIS AOD."""
    q = mdbx.Query('SATAOD', ELEMENT, area=AREA,
                   start=START, stop=STOP,
                   constrain=constrain, keep=True)
    return q.plot(ELEMENT, **kw)


def plot_sataod_test(ELEMENT='AOD_NM550', AREA=None, START=None, STOP=None,
                     constrain=None, **kw):
    """Plot MODIS AOD from mdb-test server."""
    q = mdbx.Query('SATAOD', [ELEMENT, 'STLT_IDNY'], area=AREA,
                   start=START, stop=STOP,
                   ddict='\"/path/to/tests/SATAOD/retrieval_table\"',  # update
                   hostname='mdb-test',
                   constrain=constrain, keep=True)
    return q.plot(ELEMENT, **kw)


def plot_sataod_arabian_peninsula():
    """Plot SATAOD over the Arabian Peninsula."""
    from ypylib.utils import seqdate, log
    save_dir = os.path.expandvars('$SCRATCH') + '/dust/20150820-20150920'
    for date in seqdate('2015-08-20', '2015-09-20',
                        in_fmt='%Y-%m-%d', out_fmt='%Y%m%d'):
        log.info(date)

        plt = plot_sataod(
            AREA=['38N', '12.5N', '34E', '67E'],
            START='/'.join([date, '0000Z']), STOP='/'.join([date, '2359Z']),
            plt_type='scatter', markersize=5,
            vmin=0, vmax=2,
            gspacing=(5, 5),
            map_res='l',
            drawcountries=True,
            title='MODIS aerosol optical depth (DT+DB): ' + date,
            cb_on=True, cb_title='AOD at 550nm []')
        plt.savefig('/'.join([save_dir, 'SATAOD_' + date + '.png']))
        plt.close()


def plot_msgradfd_geo():
    """MSG Cloud Top Temperature Full Disc."""
    import cartopy.crs as ccrs
    SUBTYPE = 'MSGRADFD'
    ELEMENT = 'CLOD_TOP_TMPR'
    req = mdbx.Query(
        SUBTYPE, ELEMENT,
        area=['67N', '67S', '67W', '67E'],  # must specify
        start='TODAY-1/1000Z',
        stop='TODAY-1/1230Z', keep=True)

    return req.plot(
        ELEMENT, use_cartopy=True,
        projection=ccrs.Geostationary(),
        delta=(.5, .5),
        globe=True, cb_on=True,)


def plot_sataod_india_dust():
    """SATAOD India dust case"""
    plt = plot_sataod(  # use_cartopy=True, map_res='h',
        AREA=['40N', '15N', '60E', '90E'],
        START='20180502/0000Z', STOP='20180502/2359Z',
        delta=(0.25, 0.25),
        describe_data=True,
        # vmax=2,
        # plt_type='contour',
        map_res='l',
        drawcountries=True,
        gspacing=(10, 10),
        cb_on=True, cb_title='Aerosol Optical Depth at 550nm')

    # -- Add location markers
    locs = {
        'Delhi': {'lat': 28.7041, 'lon': 77.1025},
        'Agra': {'lat': 27.1767, 'lon': 78.0081},
    }
    arrow = dict(facecolor='black', width=1, shrink=0.1, headwidth=5)
    for i in list(locs.keys()):
        plt.scatter(
            locs[i]['lon'], locs[i]['lat'],
            marker='s', color='k', alpha=0.5)
        plt.annotate(
            i, xy=(locs[i]['lon'], locs[i]['lat']),
            xytext=(locs[i]['lon'] + 3, locs[i]['lat'] + 3.5),
            arrowprops=arrow, fontsize='11')
    return plt


def plot_crete_dust():
    """SATAOD Crete dust case."""
    date = '20180325'
    plt = plot_sataod(
        AREA=['42N', '25N', '12E', '35E'],
        START='/'.join([date, '0000Z']), STOP='/'.join([date, '2359Z']),
        # delta=(0.1, 0.1), plt_type='scatter',
        delta=(.25, .25),
        vmin=0, vmax=3,
        map_res='l',
        drawcountries=True,
        gspacing=(5, 5),
        title='MODIS Aqua+Terra aerosol optical depth (DT+DB) on ' + date,
        cb_on=True, cb_title='AOD at 550nm []')
    return plt


def plot_wow_air_temp():
    """Plot WOW surface air temperature."""
    SUBTYPE = 'WOW'
    ELEMENTS = 'SRFC_AIR_TMPR'
    AREA = ['63N', '49N', '12W', '4E']

    q = mdbx.Query(SUBTYPE, ELEMENTS, area=AREA)
    return q.plot(
        ELEMENTS, delta=(0.2, 0.2), cmap='viridis', cb_on=True,
        use_cartopy=True, map_res='h',
        gspacing=(4, 4),
        describe_data=True,
        vmin=270, vmax=290)


def plot_ascat_mergde_model_field():
    """Plot merged model data in ASCAT subtype."""
    SUBTYPE = 'ASCAT'
    MERGED_FIELD = 'MODL_SRFC_HGHT'
    AREA = ['63N', '49N', '12W', '4E']

    # Merged fields are not defined in the subtype table, so first we do it:
    # 1. manually:
    # >>> subtypes.DTYPE_MAPS[SUBTYPE][MERGED_FIELD] = 'f4'
    #
    # or 2. using set_element_dtype() method:
    # >>> req = mdbx.Query(SUBTYPE, ELEMENTS, ...)
    # >>> req.set_elem      ent_dtype(ELEMENTS, dtypes)
    #
    # or 3. using fix_unknown keyword in extrtact() or plot() methods:
    # >>> req.extract(fix_unknown=True)
    # >>> req.plot(fix_unknown=True)

    # Then send the MetDB request with mergeback option enabled
    req = mdbx.Query(SUBTYPE, MERGED_FIELD, stop='TODAY-1/2345Z', area=AREA,
                     keep=True, merged=True)
    map_title = '{}:{} {}-{}'.format(
        SUBTYPE, MERGED_FIELD, req.start, req.stop)

    plot_kw = dict(
        fix_unknown=True,
        title=map_title,
        cb_on=True, cb_title='model surface height [m]',
        use_cartopy=True,
        gspacing=(4, 4),
        describe_data=True,
        cmap='viridis',
        plt_type='tricontour',
        c_levels=30,
        show=True,
    )
    req.plot(MERGED_FIELD, **plot_kw)


def plot_argo(subtype='ARGOB', elements=(('SALNY',), 100), platform='029',
              start='20180121/0000Z', stop='20180121/1559Z'):
    """
    Extract salinity profiles from ARGO buoy observation network data

    Parameters
    ----------
    subtype : str, optional
        MetDB subtype name:
        * 'ARGO' : for observation-only data from 2008-10-27 to 2018-07-01
        * 'ARGOB' : for observation-only data from 2017-10-16 onwards
    elements : tuple, optional
        Field names and number of records to extract
    platform : str, optional
        Platform id
    start : str, optional
        Start date and time
    stop : str, optional
        End date and time

    """
    req = mdbx.Query(subtype, elements, platform=platform, start=start,
                     stop=stop)
    req.plot('SALNY', index=0, title='salinity',
             plt_type='scatter', show=True, valid_min=10, map_buffer=5)


def plot_atdnet(subtype='ATDNET', elements='LGHN_STRK',
                area=['40N', '0N', '60E', '100E'],
                start='20190112/0000Z', stop='20190210/2300Z'):
    req = mdbx.Query(subtype, elements, start=start, stop=stop, area=area,
                     keep=True)
    req.plot(elements, plt_type='scatter', s_marker='+', s_lw=0.1,
             s_color='r', use_cartopy=True, drawcountries=True,
             drawstates=True,
             map_res='l', show=True)


def plot_snow_depth():
    """Plot Snow Depth show only 0 values."""

    req = mdbx.Query(
        'LNDSYB', 'SNOW_DPTH', area=['90N', '20N', '20W', '80E'],
        start='20171220/0000Z', stop='20171220/2359Z')

    plt_kw = dict(
        plt_type='scatter',
        s_marker='+', s_ms=5, s_color='k',
        cb_on=False,
        use_cartopy=True, map_res='h',
        valid_min=0, valid_max=0,
        show=True)

    req.plot('SNOW_DPTH', cb_title='snow depth [m]', **plt_kw)


def extract_giirs():
    """
    Define subtype/elements and extract GIIRS data from mdb-test.

    GIIRS Radiance stored at 1650 channels
    """
    SUBTYPE = 'GIIRS'
    ELEMENTS = ['LNGD_LW', 'LTTD_LW', (('CHNL_RDNC',), 1650)]
    AREA = ['90N', '90S', '180W', '180E']
    START, STOP = '20200130/0500Z', '20200130/0959Z'

    # GIIRS entries are not in subtype definition, so update manually
    GIIRS_dict = {
        'GIIRS': {
            'BAND_TYPE': 'i4',  # band (set to missing)
            'BAND_TYPE_LW': 'i4',  # band (2=LW)
            'BAND_TYPE_MW': 'i2',  # band (3=MW)
            'CHNL_END_LW': 'f4',  # endChannel
            'CHNL_END_MW': 'f4',  # endChannel
            'CHNL_NMBR': 'i4',  # channelNumber
            'CHNL_RDNC': 'f4',
            'CHNL_RDNC_NSE': 'f4',  # channelRadiance (use for noise, NEdR)
            'CHNL_RPLTN_CONT': 'i4',  # extndDelayedDescriptorReplicationFactor
            'CHNL_STRT_LW': 'f4',  # startChannel
            'CHNL_STRT_MW': 'f4',  # startChannel
            # confidenceFlag (dataset "QF_LWElementExploration").
            # 0=valid, 1=invalid, 15=missing
            'CNFC_FLAG_LW': 'i4',
            # confidenceFlag (dataset "QF_MWElementExploration").
            # 0=valid, 1=invalid, 15=missing
            'CNFC_FLAG_MW': 'i4',
            'DAY': 'f4',
            'FOR_NMBR': 'f4',  # fieldOfRegardNumber (attribute "Dwell number")
            'FOV_NMBR': 'i4',
            'HOUR': 'i4',
            # l1ProcessingFlag (attribute "Data Quality"). 0=OK, 1=not OK
            'L1_PRCSG_FLAG': 'i2',
            'LNGD': 'f4',
            'LNGD_LW': 'f4',
            'LNGD_MW': 'f4',
            'LTTD': 'f4',
            'LTTD_LW': 'f4',
            'LTTD_MW': 'f4',
            'MINT': 'i4',
            'MNTH': 'i4',
            'ORGNG_GNRTG_CNTR': 'S20',  # centre (CCT C-1)
            'ORGNG_GNRTG_CNTR2': 'S20',  # subCentre (CCT C-12)
            'RAD_FLGS': 'i4',  # radianceTypeFlags (4=apodized, 5=unapodized)
            'RCPT_DAY': 'i4',
            'RCPT_HOUR': 'i4',
            'RCPT_MINT': 'i4',
            'RCPT_MNTH': 'i4',
            'RCPT_YEAR': 'i4',
            'SCAN_LINE_NMBR': 'i4',
            'SCND': 'i4',
            'SOLR_AZMH_LW': 'f4',  # solarAzimuth
            'SOLR_AZMH_MW': 'f4',  # solarAzimuth
            'SOLR_ZNTH_ANGL_LW': 'f4',  # solarZenithAngle
            'SOLR_ZNTH_ANGL_MW': 'f4',  # solarZenithAngle
            'STLT_AZMH_LW': 'f4',  # bearingOrAzimuth
            'STLT_AZMH_MW': 'f4',  # bearingOrAzimuth
            'STLT_CLAS': 'i4',  # satelliteClassification (383=FY-4)
            'STLT_IDNY': 'i4',
            'STLT_INST': 'i4',  # satelliteInstruments (962=GIIRS, CCT C-8)
            'STLT_ZNTH_ANGL_LW': 'f4',  # satelliteZenithAngle
            'STLT_ZNTH_ANGL_MW': 'f4',  # satelliteZenithAngle
            # heightOfStation (m, to nearest 100m, geostationary height range)
            'STTN_HGHT': 'f4',
            'WAVE_NMBR_END_LW': 'f4',  # waveNumber (end)
            'WAVE_NMBR_END_MW': 'f4',  # waveNumber (end)
            'WAVE_NMBR_STRT_LW': 'f4',  # waveNumber (start)
            'WAVE_NMBR_STRT_MW': 'f4',  # waveNumber (start)
            'YEAR': 'i4',
        }
    }
    mdbx.subtypes.DTYPE_MAPS.update(GIIRS_dict)

    # -- Create a Query instance (update actual ddict path)
    req = mdbx.Query(
        SUBTYPE, ELEMENTS, area=AREA, start=START, stop=STOP,
        ddict='\"/path/to/tests/GIIRS/tables/retrieval_table_GIIRS\"',
        hostname='mdb-test'
    )

    data = req.extract(fix_unknown=True)
    return data

# req = Query('MWHS', 'STLT_IDNY', start='20190621/1300Z')
# print(req.extract())


# sataod: 783 (Terra), 784 (Aqua)
# data from test server
# plot_sataod_test(
#     constrain={'STLT_IDNY': 783},
#     START='20170928/0000Z',
#     STOP='20170928/2300Z',
#     delta=(0.5, 0.5),
#     cb_on=True, use_cartopy=True).show()

# data from operational server
# import cartopy.crs as ccrs
# plot_sataod(constrain={'SRFC_TYPE': 3},
#             # constrain={'STLT_IDNY': 783, 'SRFC_TYPE': 3},
#             # AREA=['40N', '10N', '30E', '60E'],
#             AREA=['40N', '1N', '30E', '100E'],
#             START='20190314/0600Z',
#             STOP='20190314/0900Z',
#             delta=(0.2, 0.2),
#             # plt_type='hexbin',
#             # projection=ccrs.Robinson(),
#             vmin=0, vmax=2,
#             cb_on=True, use_cartopy=True, show=True)

# plot_sataod(use_cartopy=True, cb_on=True, map_res='h').show()

# AREA = ['60N', '30N', '30W', '50E']

# q = Query('SATAOD', 'AOD_NM550')
# data = q.extract()
# print(data['AOD_NM550'])

# # ----------------------------------------------------------------------
# SUBTYPE = 'SMOS'
# ELEMENTS = ['BRGTS_TMPR_REAL', 'WATR_FRCN', 'PLRZN']
# VMIN, VMAX = 100, 300
# # ----------------------------------------------------------------------

# SUBTYPE = 'MSGRAD'
# # ELEMENT = 'CLOD_TOP_TMPR'
# ELEMENTS = (('CSR_STLT_INST', 'CSR_CNTR_FRQY', 'CSR_SPCL_RDNC',
#              'CSR_RDNC', 'CSR_BRGTS_TMPR'), 12)

# # AREA = ['63N', '49N', '12W', '4E']
# # # VMIN, VMAX = 270, 290
# # ----------------------------------------------------------------------

# # subtypes.DTYPE_MAPS[SUBTYPE][ELEMENT] = 'f4'
# req = Query(SUBTYPE, ELEMENTS,
#             area=AREA,
#             start='TODAY/1000Z',
#             stop='TODAY/1030Z',
#             # merged=False,
#             keep=True)

# data = req.extract()

# # print data['CSR_CNTR_FRQY'][0][7]

# # print data['CSR_BRGTS_TMPR'][0]
# # print data['CSR_RDNC'][:, 6].max()

# req.plot(ELEMENTS[0][3], index=7, delta=0.25,
#          use_cartopy=True, map_res='h',
#          cb_on=True, show=True)

# # req.plot(ELEMENTS,
# #          use_cartopy=True,
# #          # method='max',
# #          # map_res='h',
# #          delta=(0.15, 0.15),
# #          # vmin=VMIN, vmax=VMAX,
# #          # plt_type='contour',
# #          cb_on=True).show()


if __name__ == "__main__":
    # plot_atdnet(start='TODAY-1/0000Z', stop='TODAY/0000Z')
    # plot_argo()
    # plot_ascat_mergde_model_field()
    # plot_wow_air_temp().show()
    # plot_sataod_india_dust().show()
    # plot_crete_dust().show()
    # plot_msgradfd_geo().show()

    pass
