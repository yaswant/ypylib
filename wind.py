#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple scalar wind conversion

"""
import numpy as np


class ScalarWind(object):

    """
    ScalarWind Class method

    Attributes:
        direction (TYPE): Description
        speed (TYPE): Description
    """

    def __init__(self, speed, direction=None):
        super(ScalarWind, self).__init__()
        self.speed = np.array(speed)
        self.direction = np.array(direction)

    def uv(self):
        """
        Return zonal (u) and meridional (v) components of vector wind from
        windspeeds and directions (degrees). u and v have same units as
        windspeed.
        """
        u = -self.speed * np.sin(np.radians(self.direction))
        v = -self.speed * np.cos(np.radians(self.direction))
        return u, v

    def stress(self, Cd=None, rho_air=1.3):
        """
        Return wind stress components (tau_u, tau_v) from wind speeds and
        directions.
        """
        if Cd is None:
            Cd = self.drag_coeff()

        # get u and v components of wind
        u, v = self.uv()

        tau_x = rho_air * Cd * np.abs(self.speed) * u
        tau_y = rho_air * Cd * np.abs(self.speed) * v

        return tau_x, tau_y

    def drag_coeff(self):
        """
        Return drag coefficient (Cd) over ocean from 10m windspeed using
        Yelland and Taylor (1996) method.
        """
        spd = np.atleast_1d(self.speed)
        med = np.where((self.speed > 3.0) & (self.speed <= 6.0))
        high = np.where((self.speed > 6.0) & (self.speed <= 26.0))

        c_drag = np.zeros_like(spd)
        c_drag[:] = 1e-4

        if len(med) > 0:
            c_drag[med] = (0.29 + (3.1 / spd[med]) +
                           (7.7 / (spd[med]**2))) * 1e-3

        if len(high) > 0:
            c_drag[high] = (0.06 + (0.07 / spd[high])) * 1e-3

        return c_drag
