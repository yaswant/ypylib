#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Custom matplotlib configuration

"""
import matplotlib as mpl
from sys import platform
import os
import socket

__version__ = "0.1"
__author__ = "Yaswant Pradhan"
__copyright__ = "(c) Crown copyright 2019, the Met Office."


def load():
    # - automatically switch backend to Agg if using a non-interactive or
    #   system with no DISPLAY environment
    if platform.startswith('linux'):
        if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
            mpl.use('Agg')

    # - avoid backend confilct
    mpl.rcParams["backend_fallback"] = True

    # -- font (only ttf, afm), grid
    def_sans = mpl.rcParamsDefault['font.sans-serif']
    def_mono = mpl.rcParamsDefault['font.monospace']
    font = {
        'family': 'sans-serif',
        'sans-serif': ['Helvetica', 'Arial', 'Fira Sans'] + def_sans,
        'monospace': ['Consolas', 'Courier New', 'Fira Mono'] + def_mono,
    }
    mpl.rc('font', **font)
    mpl.rc('axes', grid=True, axisbelow=True)
    mpl.rc(('xtick', 'ytick'), direction='in')
    mpl.rc('lines', markersize=5)


def reset():
    """Load matplotlib default rc parameters"""
    mpl.rcdefaults()
