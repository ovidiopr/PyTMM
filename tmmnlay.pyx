# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2014-2015 Pavel Dmitriev <pavel.a.dmitriev@gmail.com>
#    Copyright (C) 2018 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of tmmnlay.
#
#    tmmnlay was forked from PyTMM, which was developed by Pavel Dmitriev
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
from enum import Enum


class Polarization(Enum):
    s = 0
    p = 1


class MultiLayer(object):
    """
        MultiLayer TMM

        Problem geometry:

        | 1 |   |      | |    | |      |     |    | |          |   | T |
        |   | = | I_01 | | L1 | | I_12 | ... | Ln | | I_n(n+1) | = |   |
        | R |   |      | |    | |      |     |    | |          |   | 0 |

        In addition

    """
    def __init__(self, n, d, wvl, aoi=0.0):
        """
        Initialize a MultiLayer class.
        
        Description:
          This class contains the functions needed to solve the wave
          propagation in a multilayer structure.
          
        Inputs:
          n -- A 2D array of possibly complex values for the index of
               refraction of the layers. The first dimension is for the
               wavelengths and the second for the layers
          d -- A 1D array of thickness for all layers. d[0] and d[-1]
               (corresponding to the outer semi-infinite layers) will be
               set to zero, regardless of whatever value they have before you
               pass the array to this class.  If you don't want this done to
               your array, make a copy first. (units must match that of the 
               wavelength)
          wvl -- A 1D array of wavelength (units must match that of the
               thickness)
          aoi -- Angle of incidence (in degrees)
        """
        # Make sure that the list of refractive indexes is a numpy array
        if type(n) is np.ndarray:
            self._n = n
        elif type(n) in (int, float, complex):
            self._n = np.array([n])
        else:
            self._n = np.array(n)

        # Make sure that the list of thicknesses is a numpy array
        if type(d) is np.ndarray:
            self._d = d
        elif type(d) in (int, float, complex):
            self._d = np.array([d])
        else:
            self._d = np.array(d)
        # Enforce this requirement for the outer layers
        self._d[0] = self._d[-1] = 0.0

        # Make sure that the list of wavelengths is a numpy array
        if type(wvl) is np.ndarray:
            self._wvl = wvl
        elif type(wvl) in (int, float, complex):
            self._wvl = np.array([wvl])
        else:
            self._wvl = np.array(wvl)

        # If n is an 1D array then assume that it is constant
        # then, we repeat the values for all wavelengths
        if len(self._n.shape) == 1:
            ln = self._n.shape[0]
            lw = self._wvl.shape[0]
            self._n = np.tile(self._n, lw).reshape((lw, ln))

        self._aoi = aoi
        self._sin2 = np.sin(self._aoi*np.pi/180.0)**2.0

        # Make sure that the problem is well defined
        assert self._d.shape[0] >= 2, "We need at least two layers (%i found)" % (self._d.shape[0])
        assert self._n.shape[0] == self._wvl.shape[0],\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self._n.shape[0], self._wvl.shape[0])
        assert self._n.shape[1] == self._d.shape[0],\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self._n.shape[1], self._d.shape[0])

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffi_TE = self._coeffi_TE = None

    @property
    def n(self):
        """Returns the refractive index values."""
        return self._n

    @n.setter
    def n(self, value):
        """Updates the refractive index values."""
        if type(value) is np.ndarray:
            self._n = value
        elif type(value) in (int, float, complex):
            self._n = np.array([value])
        else:
            self._n = np.array(value)

        # If n is an 1D array then we repeat the values for all wavelengths
        if len(self._n.shape) == 1:
            ln = self._n.shape[0]
            lw = self._wvl.shape[0]
            self._n = np.tile(self._n, lw).reshape((lw, ln))

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffi_TE = self._coeffi_TE = None

    @property
    def d(self):
        """Returns the thickness values."""
        return self._d

    @d.setter
    def d(self, value):
        """Updates the thickness values."""
        if type(value) is np.ndarray:
            self._d = value
        elif type(value) in (int, float, complex):
            self._d = np.array([value])
        else:
            self._d = np.array(value)
        # Enforce this requirement for the outer layers
        self._d[0] = self._d[-1] = 0.0

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffi_TE = self._coeffi_TE = None

    @property
    def wvl(self):
        """Returns the wavelength values."""
        return self._wvl

    @wvl.setter
    def wvl(self, value):
        """Updates wavelength values."""
        if type(value) is np.ndarray:
            self._wvl = value
        elif type(value) in (int, float, complex):
            self._wvl = np.array([value])
        else:
            self._wvl = np.array(value)

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffi_TE = self._coeffi_TE = None

    @property
    def aoi(self):
        """Returns the angle of incidence (degrees)."""
        return self._aoi

    @aoi.setter
    def aoi(self, value):
        """Updates the angle of incidence (degrees)."""
        self._aoi = value
        self._sin2 = np.sin(self._aoi*np.pi/180.0)**2.0

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffi_TE = self._coeffi_TE = None

    @property
    def num_layers(self):
        """Returns the number of layers."""
        return self._d.shape[0]

    @property
    def num_lambda(self):
        """Returns the number of wavelength values."""
        return self._wvl.shape[0]

    @property
    def interface_matrices(self):
        """
        Calculates (only if necessary) and returns the interface matrices.
        To save memory it only saves the rjk and tjk coefficients corresponding
        to the TE (s) and TM (p) polarizations.
        """
        # Make sure that the problem is well defined
        assert self._d.shape[0] >= 2, "We need at least two layers (%i found)" % (self._d.shape[0])
        assert self._n.shape[0] == self._wvl.shape[0],\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self._n.shape[0], self._wvl.shape[0])
        assert self._n.shape[1] == self._d.shape[0],\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self._n.shape[1], self._d.shape[0])

        if self._im is None:
            n2 = self.n**2.0
            n2j = n2[:, :-1]
            n2k = n2[:, 1:]

            q = np.sqrt(n2 - (n2[:, 0]*self._sin2)[:, None])
            qj = q[:, :-1]
            qk = q[:, 1:]

            num1 = qj + qk
            num2 = n2k*qj + n2j*qk

            self._im = np.zeros((self.num_lambda, self.num_layers - 1, 4), dtype=complex)
            # rjk for the TE (s) polarization
            self._im[:, :, 0] = (qj - qk)/num1
            # tjk for the TE (s) polarization
            self._im[:, :, 1] = 2.0*qj/num1
            # rjk for the TM (p) polarization
            self._im[:, :, 2] = (n2k*qj - n2j*qk)/num2
            # tjk for the TM (p) polarization
            self._im[:, :, 3] = 2.0*self.n[:, :-1]*self.n[:, 1:]*qj/num2

        return self._im

    @property
    def layer_matrices(self):
        """
        Calculates (only if necessary) and returns the layer matrices.
        To save memory it only saves the values of exp(-i*beta_j) and
        exp(i*beta_j).
        """
        # Make sure that the problem is well defined
        assert self._d.shape[0] >= 2, "We need at least two layers (%i found)" % (self._d.shape[0])
        assert self._n.shape[0] == self._wvl.shape[0],\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self._n.shape[0], self._wvl.shape[0])
        assert self._n.shape[1] == self._d.shape[0],\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self._n.shape[1], self._d.shape[0])

        if (self._lm is None) and (self.num_layers > 2):
            n2 = self.n**2.0

            q = np.sqrt(n2 - (n2[:, 0]*self._sin2)[:, None])
            qj = q[:, 1:-1]

            Bj = 2.0j*np.pi*qj*self.d[1:-1, None]/self.wvl[:, None]

            self._lm = np.zeros((self.num_lambda, self.num_layers - 2, 2), dtype=complex)
            # exp(-i*beta_j)
            self._lm[:, :, 0] = np.exp(-Bj)
            # exp(i*beta_j)
            self._lm[:, :, 1] = np.exp(Bj)

        return self._lm

    @property
    def matrix_TE(self):
        """
        Calculate total TE matrix for the multilayer
        """
        if self._matrix_TE is None:
            lm = self.layer_matrices
            im = self.interface_matrices

            self._matrix_TE = np.zeros((self.num_lambda, 2, 2), dtype=complex)

            self._matrix_TE[:, 0, 0] = self._matrix_TE[:, 1, 1] = 1.0/im[:, 0, 1]
            self._matrix_TE[:, 0, 1] = self._matrix_TE[:, 1, 0] = im[:, 0, 0]/im[:, 0, 1]

            for j in range(1, self.num_layers - 1):
                B11 = lm[:, j - 1, 0]/im[:, j, 1]
                B12 = lm[:, j - 1, 0]*im[:, j, 0]/im[:, j, 1]
                B21 = lm[:, j - 1, 1]*im[:, j, 0]/im[:, j, 1]
                B22 = lm[:, j - 1, 1]/im[:, j, 1]

                C11 = self._matrix_TE[:, 0, 0]*B11 + self._matrix_TE[:, 0, 1]*B21
                C12 = self._matrix_TE[:, 0, 0]*B12 + self._matrix_TE[:, 0, 1]*B22
                C21 = self._matrix_TE[:, 1, 0]*B11 + self._matrix_TE[:, 1, 1]*B21
                C22 = self._matrix_TE[:, 1, 0]*B12 + self._matrix_TE[:, 1, 1]*B22

                self._matrix_TE[:, 0, 0] = C11
                self._matrix_TE[:, 0, 1] = C12
                self._matrix_TE[:, 1, 0] = C21
                self._matrix_TE[:, 1, 1] = C22

        return self._matrix_TE

    @property
    def matrix_TM(self):
        """
        Calculate total TM matrix for the multilayer
        """
        if self._matrix_TM is None:
            lm = self.layer_matrices
            im = self.interface_matrices

            self._matrix_TM = np.zeros((self.num_lambda, 2, 2), dtype=complex)

            self._matrix_TM[:, 0, 0] = self._matrix_TM[:, 1, 1] = 1.0/im[:, 0, 3]
            self._matrix_TM[:, 0, 1] = self._matrix_TM[:, 1, 0] = im[:, 0, 2]/im[:, 0, 3]

            for j in range(1, self.num_layers - 1):
                B11 = lm[:, j - 1, 0]/im[:, j, 3]
                B12 = lm[:, j - 1, 0]*im[:, j, 2]/im[:, j, 3]
                B21 = lm[:, j - 1, 1]*im[:, j, 2]/im[:, j, 3]
                B22 = lm[:, j - 1, 1]/im[:, j, 3]

                C11 = self._matrix_TM[:, 0, 0]*B11 + self._matrix_TM[:, 0, 1]*B21
                C12 = self._matrix_TM[:, 0, 0]*B12 + self._matrix_TM[:, 0, 1]*B22
                C21 = self._matrix_TM[:, 1, 0]*B11 + self._matrix_TM[:, 1, 1]*B21
                C22 = self._matrix_TM[:, 1, 0]*B12 + self._matrix_TM[:, 1, 1]*B22

                self._matrix_TM[:, 0, 0] = C11
                self._matrix_TM[:, 0, 1] = C12
                self._matrix_TM[:, 1, 0] = C21
                self._matrix_TM[:, 1, 1] = C22

        return self._matrix_TM

    @property
    def rt_TE(self):
        """
        """
        S = self.matrix_TE
        r = S[:, 1, 0]/S[:, 0, 0]
        t = 1.0/S[:, 0, 0]
        #print r, t

        return r, t

    @property
    def rt_TM(self):
        """
        """
        S = self.matrix_TM
        r = S[:, 1, 0]/S[:, 0, 0]
        t = 1.0/S[:, 0, 0]

        return r, t

