#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Simple scalar wind conversion

"""
import numpy as np


class ScalarWind(object):

    """ScalarWind Class method.

    Attributes:
        direction (TYPE): Description
        speed (TYPE): Description
    """

    def __init__(self, speed, direction=None):
        super(ScalarWind, self).__init__()
        self.speed = np.array(speed)
        self.direction = np.array(direction)

    def uv(self):
        """Return zonal (u) and meridional (v) components of vector wind from
        windspeeds and directions (degrees).

        u and v have same units as windspeed.
        """
        u = -self.speed * np.sin(np.radians(self.direction))
        v = -self.speed * np.cos(np.radians(self.direction))
        return u, v

    def stress(self, Cd=None, rho_air=1.3):
        """Return wind stress components (tau_u, tau_v) from wind speeds and
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
        """Return drag coefficient (Cd) over ocean from 10m windspeed using
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

    def divergence_cfd(self, dx=1, dy=1):
        # TODO: compute divergence via centred finite difference

        # div = Dv/Dy + Du/Dx -(v/a)*tan(phi)
        # div(j,i) = (v(j+1,i)-v(j-1,i))/dy2(j)
        #       + (u(j,i+1)-u(j,i-1))/dx2(j)
        #       - (v(j,i)/a)*tan(phi(j))

        # compare against metpy
        # from metpy import calc
        # from metpy.units import units
        # print(calc.kinematics.divergence(u, v, 1 * units.m, 1 * units.m))

        dudx = np.gradient(self.u, dx, axis=1)
        dvdy = np.gradient(self.v, dy, axis=0)
        return dudx + dvdy

        pass

    def vorticity_cdf(self, dx, dy):
        # TODO: compute vorticity via centred finite difference

        # vor = Dv/Dx - Du/Dy + (u/a)*tan(phi)

        # vor(j,i) = (v(j,i+1)-v(j,i-1))/dx2(j)
        #       - (u(j+1,i)-u(j-1,i))/dy2(j)
        #       + (u(j,i)/a)*tan(phi(j))
        pass


def divergence(u, v, dx=1, dy=1):
    dudx = np.gradient(u, dx, axis=1)
    dvdy = np.gradient(v, dy, axis=0)
    print(u, v)
    return dudx + dvdy


if __name__ == "__main__":

    u = np.array([[1, 2, 6], [3, 4, 5], [3, 4, 5]])
    v = np.array([[4, 5, 7], [1, 4, 5], [3, 4, 5]])

    print(divergence(u, v, 1, 1))
    # print(divergence(u, v, 0.5))
