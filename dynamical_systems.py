'''
Author: Edoardo Caldarelli
Affiliation: Institut de Robòtica i Informàtica Industrial, CSIC-UPC
email: ecaldarelli@iri.upc.edu
July 2024
'''

from abc import ABC, abstractmethod
import numpy as np

class DynamicalSystem(ABC):
    @abstractmethod
    def update_SOM(self, xbef, control):
        pass


class DuffingOscillator(DynamicalSystem):
    """
    Duffing oscillator with Runge-Kutta discretization.
    """
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def _f_u(self, x, u):
        return -np.vstack((-x[1, :],
                          0.5 * x[1, :]+ x[0, :] * (4 * x[0, :] ** 2 - 1) - 0.5 * u))

    def _f_ud(self, x, u):
        return (x + (self.Ts / 6) * (self._k1(x, u) + 2 * self._k2(x, u) + 2 * self._k3(x, u) +
                                     self._k4(x, u)))

    def _k1(self, x, u):
        return self._f_u(x, u)

    def _k2(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts / 2, u)

    def _k3(self, x, u):
        return self._f_u(x + self._k2(x, u) * self.Ts / 2, u)

    def _k4(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts, u)

    def update_SOM(self, xbef, u):
        if len(xbef.shape) == 1:
            xbef = xbef.reshape([-1, 1])
        return self._f_ud(xbef, u)

class DoubleIntegrator(DynamicalSystem):
    """
    Double integrator with Runge-Kutta discretization.
    """
    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def _f_u(self, x, u):
        return np.vstack((x[1, :],
                           u))

    def _f_ud(self, x, u):
        return (x + (self.Ts / 6) * (self._k1(x, u) + 2 * self._k2(x, u) + 2 * self._k3(x, u) +
                                     self._k4(x, u)))

    def _k1(self, x, u):
        return self._f_u(x, u)

    def _k2(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts / 2, u)

    def _k3(self, x, u):
        return self._f_u(x + self._k2(x, u) * self.Ts / 2, u)

    def _k4(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts, u)

    def update_SOM(self, xbef, u):
        if len(xbef.shape) == 1:
            xbef = xbef.reshape([-1, 1])
        return self._f_ud(xbef, u)

class HJB(DynamicalSystem):
    """
    Tutorial system from Guo et al. 2022 with Runge-Kutta discretization.
    """

    def __init__(self, **kwargs):
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def _f_u(self, x, u):
        return -x ** 3 + u

    def _f_ud(self, x, u):
        return (x + (self.Ts / 6) * (self._k1(x, u) + 2 * self._k2(x, u) + 2 * self._k3(x, u) +
                                     self._k4(x, u)))

    def _k1(self, x, u):
        return self._f_u(x, u)

    def _k2(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts / 2, u)

    def _k3(self, x, u):
        return self._f_u(x + self._k2(x, u) * self.Ts / 2, u)

    def _k4(self, x, u):
        return self._f_u(x + self._k1(x, u) * self.Ts, u)

    def update_SOM(self, xbef, u):
        return self._f_ud(xbef, u)