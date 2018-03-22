# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2018 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of tmmnlay.
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

(TE, TM) = range(2)

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
        Initialize the MultiLayer class.
        
        Description:
          This class contains the functions needed to solve the wave
          propagation in a multilayer structure.
          
        Inputs:
          n -- A 1D or 2D array of possibly complex values for the index of
               refraction of the layers. The first dimension is for the
               wavelengths and the second for the layers. If the array is
               1D  then assume that it is constantand repeat the values
               for all wavelengths
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
        self.d = d
        self.wvl = wvl
        self.n = n
        self.aoi = aoi

        # Make sure that the problem is well defined
        assert self.num_layers >= 2, "We need at least two layers (%i found)" % (self.num_layers)
        assert self.num_lambda_n == self.num_lambda,\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self.num_lambda_n, self.num_lambda)
        assert self.num_layers_n == self.num_layers,\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self.num_layers_n, self.num_layers)

    @property
    def num_layers(self):
        """Returns the number of layers."""
        return self._d.shape[0]

    @property
    def num_lambda(self):
        """Returns the number of wavelength values."""
        return self._wvl.shape[0]

    @property
    def num_layers_n(self):
        """Returns the number of layers."""
        return self._n.shape[1]

    @property
    def num_lambda_n(self):
        """Returns the number of wavelength values."""
        return self._n.shape[0]

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

        # If n is an 1D array we repeat the values for all wavelengths
        if len(self._n.shape) == 1:
            ln = self._n.shape[0]
            lw = self._wvl.shape[0]
            self._n = np.tile(self._n, lw).reshape((lw, ln))

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffs_TE = self._coeffs_TM = None

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
        self._coeffs_TE = self._coeffs_TM = None

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
        self._coeffs_TE = self._coeffs_TM = None

    @property
    def aoi(self):
        """Returns the angle of incidence (degrees)."""
        return self._aoi

    @aoi.setter
    def aoi(self, value):
        """Updates the angle of incidence (degrees)."""
        self._aoi = value

        self._im = self._lm = None
        self._matrix_TE = self._matrix_TM = None
        self._coeffs_TE = self._coeffs_TM = None

    @property
    def interface_matrices(self):
        """
        Calculates (only if necessary) and returns the interface matrices.
        To save memory it only saves the rjk and tjk coefficients corresponding
        to the TE (s) and TM (p) polarizations.
        """
        # Make sure that the problem is well defined
        assert self.num_layers >= 2, "We need at least two layers (%i found)" % (self.num_layers)
        assert self.num_lambda_n == self.num_lambda,\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self.num_lambda_n, self.num_lambda)
        assert self.num_layers_n == self.num_layers,\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self.num_layers_n, self.num_layers)

        if self._im is None:
            s2 = np.sin(self.aoi*np.pi/180.0)**2.0
            n2 = self.n**2.0
            nj = self.n[:, :-1]
            nk = self.n[:, 1:]

            c = np.sqrt(1.0 - (n2[:, 0]*s2)[:, None]/n2)
            cj = c[:, :-1]
            ck = c[:, 1:]

            num1 = nj*cj + nk*ck
            num2 = nk*cj + nj*ck

            self._im = np.zeros((self.num_lambda, self.num_layers - 1, 4), dtype=complex)
            # rjk for the TE (s) polarization
            self._im[:, :, 0] = (nj*cj - nk*ck)/num1
            # tjk for the TE (s) polarization
            self._im[:, :, 1] = 2.0*nj*cj/num1
            # rjk for the TM (p) polarization
            self._im[:, :, 2] = (nk*cj - nj*ck)/num2
            # tjk for the TM (p) polarization
            self._im[:, :, 3] = 2.0*nj*cj/num2

        return self._im

    @property
    def layer_matrices(self):
        """
        Calculates (only if necessary) and returns the layer matrices.
        To save memory it only saves the values of exp(-i*beta_j) and
        exp(i*beta_j).
        """
        # Make sure that the problem is well defined
        assert self.num_layers >= 2, "We need at least two layers (%i found)" % (self.num_layers)
        assert self.num_lambda_n == self.num_lambda,\
            "The number of refractive index values (%i) does not match the number of wavelengths (%i)"\
            % (self.num_lambda_n, self.num_lambda)
        assert self.num_layers_n == self.num_layers,\
            "The number of refractive index values (%i) does not match the number of layers (%i)"\
            % (self.num_layers_n, self.num_layers)

        if (self._lm is None) and (self.num_layers > 2):
            s2 = np.sin(self.aoi*np.pi/180.0)**2.0
            n2 = self.n**2.0

            kj = np.sqrt(n2[:, 1:-1] - (n2[:, 0]*s2)[:, None])*self.d[1:-1]

            Bj = 2.0j*np.pi*kj/self.wvl[:, None]

            self._lm = np.zeros((self.num_lambda, self.num_layers - 2, 2), dtype=complex)
            # exp(-i*beta_j)
            self._lm[:, :, 0] = np.exp(-Bj)
            # exp(i*beta_j)
            self._lm[:, :, 1] = np.exp(Bj)

        return self._lm

    def matrix(self, pol=TE):
        im = self.interface_matrices

        S = np.zeros((self.num_lambda, 2, 2), dtype=complex)

        S[:, 0, 0] = S[:, 1, 1] = 1.0/im[:, 0, 2*pol + 1]
        S[:, 0, 1] = S[:, 1, 0] = im[:, 0, 2*pol]/im[:, 0, 2*pol + 1]

        if self.num_layers > 2:
            lm = self.layer_matrices

            B11 = lm[:, :, 0]/im[:, 1:, 2*pol + 1]
            B12 = lm[:, :, 0]*im[:, 1:, 2*pol]/im[:, 1:, 2*pol + 1]
            B21 = lm[:, :, 1]*im[:, 1:, 2*pol]/im[:, 1:, 2*pol + 1]
            B22 = lm[:, :, 1]/im[:, 1:, 2*pol + 1]

            for j in range(self.num_layers - 2):
                C11 = S[:, 0, 0]*B11[:, j] + S[:, 0, 1]*B21[:, j]
                C12 = S[:, 0, 0]*B12[:, j] + S[:, 0, 1]*B22[:, j]
                C21 = S[:, 1, 0]*B11[:, j] + S[:, 1, 1]*B21[:, j]
                C22 = S[:, 1, 0]*B12[:, j] + S[:, 1, 1]*B22[:, j]

                S[:, 0, 0] = C11
                S[:, 0, 1] = C12
                S[:, 1, 0] = C21
                S[:, 1, 1] = C22

        return S

    @property
    def matrix_TE(self):
        """
        Calculate total TE matrix for the multilayer
        """
        if self._matrix_TE is None:
            self._matrix_TE = self.matrix(pol=TE)

        return self._matrix_TE

    @property
    def matrix_TM(self):
        """
        Calculate total TM matrix for the multilayer
        """
        if self._matrix_TM is None:
            self._matrix_TM = self.matrix(pol=TM)

        return self._matrix_TM

    @property
    def rt_TE(self):
        """
        """
        S = self.matrix_TE
        r = S[:, 1, 0]/S[:, 0, 0]
        t = 1.0/S[:, 0, 0]

        return r, t

    @property
    def rt_TM(self):
        """
        """
        S = self.matrix_TM
        r = S[:, 1, 0]/S[:, 0, 0]
        t = 1.0/S[:, 0, 0]

        return r, t

    def coeffs_matrix(self, layer, pol=TE):
        im = self.interface_matrices
        lm = self.layer_matrices

        Sl = np.zeros((self.num_lambda, 2, 2), dtype=complex)
        Sr = np.zeros((self.num_lambda, 2, 2), dtype=complex)

        Sl[:, 0, 0] = Sl[:, 1, 1] = 1.0/im[:, 0, 2*pol + 1]
        Sl[:, 0, 1] = Sl[:, 1, 0] = im[:, 0, 2*pol]/im[:, 0, 2*pol + 1]

        Sr[:, 0, 0] = Sr[:, 1, 1] = 1.0/im[:, layer + 1, 2*pol + 1]
        Sr[:, 0, 1] = Sr[:, 1, 0] = im[:, layer + 1, 2*pol]/im[:, layer + 1, 2*pol + 1]

        B11 = lm[:, :, 0]/im[:, 1:, 2*pol + 1]
        B12 = lm[:, :, 0]*im[:, 1:, 2*pol]/im[:, 1:, 2*pol + 1]
        B21 = lm[:, :, 1]*im[:, 1:, 2*pol]/im[:, 1:, 2*pol + 1]
        B22 = lm[:, :, 1]/im[:, 1:, 2*pol + 1]

        for j in range(layer):
            C11 = Sl[:, 0, 0]*B11[:, j] + Sl[:, 0, 1]*B21[:, j]
            C12 = Sl[:, 0, 0]*B12[:, j] + Sl[:, 0, 1]*B22[:, j]
            C21 = Sl[:, 1, 0]*B11[:, j] + Sl[:, 1, 1]*B21[:, j]
            C22 = Sl[:, 1, 0]*B12[:, j] + Sl[:, 1, 1]*B22[:, j]

            Sl[:, 0, 0] = C11
            Sl[:, 0, 1] = C12
            Sl[:, 1, 0] = C21
            Sl[:, 1, 1] = C22

        for j in range(layer + 1, self.num_layers - 2):
            C11 = Sr[:, 0, 0]*B11[:, j] + Sr[:, 0, 1]*B21[:, j]
            C12 = Sr[:, 0, 0]*B12[:, j] + Sr[:, 0, 1]*B22[:, j]
            C21 = Sr[:, 1, 0]*B11[:, j] + Sr[:, 1, 1]*B21[:, j]
            C22 = Sr[:, 1, 0]*B12[:, j] + Sr[:, 1, 1]*B22[:, j]

            Sr[:, 0, 0] = C11
            Sr[:, 0, 1] = C12
            Sr[:, 1, 0] = C21
            Sr[:, 1, 1] = C22

        return Sl, Sr

    def coefficients(self, pol=TE):
        """
        Calculate the field coefficients
        """
        s2 = np.sin(self.aoi*np.pi/180.0)**2.0
        n2 = self.n**2.0
        qj = np.sqrt(n2[:, 1:-1] - (n2[:, 0]*s2)[:, None])
        Bj = 4.0j*np.pi*qj*self.d[1:-1]/self.wvl[:, None]

        coeffs = np.zeros((self.num_lambda, self.num_layers, 2), dtype=complex)

        for j in range(self.num_layers - 2):
            Sp, Spp = self.coeffs_matrix(layer=j, pol=pol)

            coeffs[:, j + 1, 0] = 1.0/(Sp[:, 0, 0] + Sp[:, 0, 1]*Spp[:, 1, 0]*np.exp(Bj[:, j])/Spp[:, 0, 0])
            coeffs[:, j + 1, 1] = coeffs[:, j + 1, 0]*Spp[:, 1, 0]*np.exp(Bj[:, j])/Spp[:, 0, 0]

        return coeffs

    @property
    def coeffs_TE(self):
        """
        Calculate the field coefficients for the TE polarization
        """
        if self._coeffs_TE is None:
            self._coeffs_TE = self.coefficients(pol=TE)

            # Define coefficients in the outer layers
            r, t = self.rt_TE
            self._coeffs_TE[:, 0, 0] = 1.0 + 0.0j
            self._coeffs_TE[:, 0, 1] = r
            self._coeffs_TE[:, -1, 0] = t
            self._coeffs_TE[:, -1, 1] = 0.0 + 0.0j

        return self._coeffs_TE

    @property
    def coeffs_TM(self):
        """
        Calculate the field coefficients for the TM polarization
        """
        if self._coeffs_TM is None:
            self._coeffs_TM = self.coefficients(pol=TM)

            # Define coefficients in the outer layers
            r, t = self.rt_TM
            self._coeffs_TM[:, 0, 0] = 1.0 + 0.0j
            self._coeffs_TM[:, 0, 1] = r
            self._coeffs_TM[:, -1, 0] = t
            self._coeffs_TM[:, -1, 1] = 0.0 + 0.0j

            # Correct the TM coefficients to ensure continuity
            self._coeffs_TM[:, :, 0] *= self.n
            self._coeffs_TM[:, :, 1] *= self.n

        return self._coeffs_TM

    def field(self, x, coeffs):
        """Returns the electric field at specified values of x.

        Inputs:
          coeffs -- Array of field coefficients
          x -- A 1D array of any length specifying the x values for which the field
               should be returned. It must be sorted in increasing order.

        Outputs: (E)
          E -- an array of field values for the given x (complex-valued)
        """    
        xl = np.cumsum(self.d)
        xl -= xl[0] # Just in case self.d[0] was not 0.0
        xl[-1] = np.inf

        s2 = np.sin(self.aoi*np.pi/180.0)**2.0
        n2 = self.n**2.0
        Zeta = 2.0j*np.pi*np.sqrt(n2 - (n2[:, 0]*s2)[:, None])/self.wvl[:, None]

        E = np.zeros((self.num_lambda, len(x)), complex)
        j = 0
        for i, xi in enumerate(x):
            while (xi > xl[j]):
                j += 1
            if (j == 0):
                xj = xi
            else:
                xj = xi - xl[j - 1]

            E[:, i] = coeffs[:, j, 0]*np.exp(Zeta[:, j]*xj) + coeffs[:, j, 1]*np.exp(-Zeta[:, j]*xj)

        return E

    def field_TE(self, x):
        """
        Calculate the electric field coefficients for the TE polarization
        """
        return self.field(x, self.coeffs_TE)

    def field_TM(self, x):
        """
        Calculate the electric field coefficients for the TE polarization
        """
        return self.field(x, self.coeffs_TM)

