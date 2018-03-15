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


class TransferMatrix:
    """
        Dielectric layer TMM

        How the functions eat structure matricies:

        | T |   |        | |        | |     |   | 1 |
        |   | = | Bottom | | Matrix | | Top | = |   |
        | 0 |   |        | |        | |     |   | R |

    """

    @staticmethod
    def structure(*args):
        """
        args - separate structure matricies
        Left to Right = Bottom to Top
        :param args:
        """
        mat = np.identity(2, dtype=np.complex128)
        for m in args:
            mat = np.dot(m.matrix, mat)
        return TransferMatrix(mat)

    @staticmethod
    def layer(n, d, wavelength, theta=0, pol=Polarization.s):
        """
        Creates a Air-DielectricLayer-Air Transfer Matrix
        :param n:
        :param d:
        :param wavelength:
        """
        bottomBoundary = TransferMatrix.boundingLayer(1, n, theta, pol)
        topBoundary = TransferMatrix.boundingLayer(n, 1, theta, pol)
        propagation = TransferMatrix.propagationLayer(n, d, wavelength, theta, pol)

        return TransferMatrix.structure(bottomBoundary, propagation, topBoundary)

    @staticmethod
    def boundingLayer(n1, n2, theta=0, pol=Polarization.s):
        """
        Creates a DielectricLayer-DielectricLayer Boundary Transfer Matrix
        :param n1:
        :param n2:
        """
        # if np.abs((n1/n2)*np.sin(theta)) >= 1.0:
        #     theta2 = np.pi/2*np.sign(np.sin(theta))
        # else:
        theta2 = np.arcsin((n1/n2)*np.sin(theta), dtype=np.complex128)

        # TE
        if pol is Polarization.s:
            _n1 = n1*np.cos(theta)
            _n2 = n2*np.cos(theta2)
            a21 = 1

        # TM
        elif pol is Polarization.p:
            _n1 = n1/np.cos(theta)
            _n2 = n2/np.cos(theta2)
            a21 = np.cos(theta2)/np.cos(theta)

        boundary = 1/(2*a21*_n2)*np.array([[(_n1 + _n2), (_n2 - _n1)],
                                    [(_n2 - _n1), (_n1 + _n2)]], dtype=np.complex128)
        return TransferMatrix(boundary)

    @staticmethod
    def propagationLayer(n, d, wavelength, theta=0, pol=Polarization.s):
        """
        Creates a Propagation Transfer Matrix, width d, refractive index n
        :param n:
        :param d:
        :param wavelength:
        """
        theta2 = np.arcsin((1/n)*np.sin(theta), dtype=np.complex128)

        propagation = np.array([[np.exp((-1j*n*d*2*np.pi/wavelength)*np.cos(theta2)), 0],
                                [0, np.exp((1j*n*d*2*np.pi/wavelength)*np.cos(theta2))]],
                               dtype=np.complex128)
        return TransferMatrix(propagation)

    def __init__(self, matrix):
        self.matrix = matrix

    def invert(self):
        """
        Inverts matrix

        """
        self.matrix = np.linalg.inv(self.matrix)

    def appendLeft(self, matrix):
        """

        :param matrix:
        """
        self.matrix = np.dot(matrix.matrix, self.matrix)

    def appendRight(self, matrix):
        """

        :param matrix:
        """
        self.matrix = np.dot(self.matrix, matrix.matrix)


def solvePropagation(transferMatrix, incidentField=1.0):
    """Calculate reflectance and transmittance
    :param transferMatrix:
    :param incidentField:
    """
    # res[1] = transmittance, res[0] = reflectance
    lhs = np.array([[transferMatrix.matrix[0, 1], -1],
                    [transferMatrix.matrix[1, 1], 0]])
    rhs = np.array([-transferMatrix.matrix[0, 0], -transferMatrix.matrix[1, 0]])
    rhs = np.multiply(rhs, incidentField)
    res = np.linalg.solve(lhs, rhs)
    reflectance = res[0]
    transmittance = res[1]
    return reflectance, transmittance


def findReciprocalTransferMatrix(transmittance, reflectance, bottomMat=TransferMatrix(np.identity(2)),
                                 topMat=TransferMatrix(np.identity(2))):  # , incidentField=1.0
    """

    :param transmittance:
    :param reflectance:
    :param bottomMat:
    :param topMat:
    :return:
    """
    assert transmittance != 0

    matrix = np.array([[1/np.conj(transmittance), reflectance/transmittance],
                       [np.conj(reflectance/transmittance), 1/transmittance]])
    matrix = np.dot(np.linalg.inv(bottomMat.matrix), matrix)
    matrix = np.dot(matrix, np.linalg.inv(topMat.matrix))
    return TransferMatrix(matrix)


def findReciprocalTransferMatrixLegacy(transmittance, reflectance, bottomMat=TransferMatrix(np.identity(2)),
                                       topMat=TransferMatrix(np.identity(2))):  # , incidentField=1.0
    """

    :param transmittance:
    :param reflectance:
    :param bottomMat:
    :param topMat:
    :return:
    """
    a = np.identity(2)
    b = np.array([[np.real(reflectance), np.imag(reflectance)],
                  [np.imag(reflectance), -np.real(reflectance)]])
    lhs = np.vstack((np.hstack((a, b)), np.hstack((b, a))))
    rhs = np.array([np.real(transmittance), np.imag(transmittance), 0, 0])
    res = np.linalg.solve(lhs, rhs)
    matrix = np.array([[res[0] + 1j*res[1], res[2] - 1j*res[3]],
                       [res[2] + 1j*res[3], res[0] - 1j*res[1]]])

    matrix = np.dot(np.linalg.inv(bottomMat.matrix), matrix)
    matrix = np.dot(matrix, np.linalg.inv(topMat.matrix))
    return TransferMatrix(matrix)


def findGeneralizedTransferMatrix(transmitance1, reflectance1, transmitance2, reflectance2,
                                  bottomMat1=TransferMatrix(np.identity(2)),
                                  topMat1=TransferMatrix(np.identity(2)),
                                  bottomMat2=TransferMatrix(np.identity(2)),
                                  topMat2=TransferMatrix(np.identity(2))):
    """

    :param transmitance1:
    :param reflectance1:
    :param transmitance2:
    :param reflectance2:
    :param bottomMat1:
    :param topMat1:
    :param bottomMat2:
    :param topMat2:
    :return:
    """
    a12 = np.dot(np.linalg.inv(bottomMat1.matrix), np.array([[transmitance1], [0]]))
    a34 = np.dot(np.linalg.inv(bottomMat2.matrix), np.array([[transmitance2], [0]]))

    b12 = np.dot(topMat1.matrix, np.array([[1], [reflectance1]]))
    b34 = np.dot(topMat2.matrix, np.array([[1], [reflectance2]]))

    rhs = np.array([a12[0, 0], a34[0, 0], a12[1, 0], a34[1, 0]])

    bmat = np.array([[b12[0, 0], b12[1, 0]],
                     [b34[0, 0], b34[1, 0]]])

    lhs = np.vstack((np.hstack((bmat, np.zeros((2, 2)))),
                     np.hstack((np.zeros((2, 2)), bmat))))
    res = np.linalg.solve(lhs, rhs)

    mat = np.array([[res[0], res[2]],
                    [res[1], res[3]]])
    return TransferMatrix(mat)
