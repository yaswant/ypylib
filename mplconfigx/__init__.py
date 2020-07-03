#!/usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from sys import platform
import os
import socket

__version__ = "0.1"
__author__ = "Yaswant Pradhan"
__copyright__ = "(c) Crown copyright 2019, the Met Office."


def load():
    """ Load custom configuration for matplotlib.

    - automatically switch backend

    If available,

    - use Helvetica/Arial/Fira Sans for sans-serif (default) font
    - use Consolas, Courier New, Fira Mono for monospace font
    - show grid lines
    - smaller markersize (5)
    - inward tick direction

    """

    # -- switch backend to Agg if using a non-interactive environment
    if platform.startswith('linux'):
        if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
            plt.switch_backend('Agg')

    # -- avoid backend confilct
    plt.rcParams["backend_fallback"] = True

    # -- custom fonts (ttf, afm), grid, axes, ticks, etc
    def_sans = plt.rcParamsDefault['font.sans-serif']
    def_mono = plt.rcParamsDefault['font.monospace']
    font = {
        'family': 'sans-serif',
        'sans-serif': ['Helvetica', 'Arial', 'Fira Sans'] + def_sans,
        'monospace': ['Consolas', 'Courier New', 'Fira Mono'] + def_mono,
    }
    plt.rc('axes', grid=True)  # axisbelow=True, titleweight='bold'
    plt.rc('font', **font)
    plt.rc('lines', markersize=5)
    plt.rc(('xtick', 'ytick'), direction='in')


def reset():
    """ Reset/unload any custom matplotlib configuration"""
    plt.rcdefaults()
