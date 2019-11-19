#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Custom matplotlib configuration

"""
from sys import platform
import os
import socket
import matplotlib as mpl

__version__ = "0.1"
__author__ = "Yaswant Pradhan"
__copyright__ = "(c) Crown copyright 2019, the Met Office."


def loadrc():
    # - automatically switch backend to Agg if using a non-interactive or
    #   system with no DISPLAY environment
    if platform.startswith('linux'):
        if 'spice' in socket.gethostname() or os.getenv('DISPLAY') is None:
            mpl.use('Agg')
    print('[info] backend: ', mpl.get_backend())

    # - avoid backend confilct
    mpl.rcParams["backend_fallback"] = True

    # - use custom sans and monospace fonts
    #   if not available, fallback to default
    mpl.rcParams["font.sans-serif"].insert(0, "Fira Sans")
    mpl.rcParams["font.monospace"] = [
        "Consolas", "Fira Mono"] + mpl.rcParamsDefault["font.monospace"]


def unloadrc():
    mpl.rcParams["backend_fallback"] = mpl.rcParamsDefault["backend_fallback"]
    mpl.rcParams["font.sans-serif"] = mpl.rcParamsDefault["font.sans-serif"]
    mpl.rcParams["font.monospace"] = mpl.rcParamsDefault["font.monospace"]
