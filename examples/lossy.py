#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2014-2015 Pavel Dmitriev <pavel.a.dmitriev@gmail.com>
#    Copyright (C) 2018 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of tmmnlay
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
import matplotlib.pyplot as plt

from tmmnlay import *

r, t = [], []
r1, t1 = [], []
r2, t2 = [], []
r3, t3 = [], []

wavelengths = np.linspace(300, 1500, 2000)
for i in wavelengths:
    a = TransferMatrix.layer(1.46, 200, i)
    b = TransferMatrix.layer(1.46 - 0.001j, 200, i)
    c = TransferMatrix.layer(1.46 - 0.01j, 200, i)
    d = TransferMatrix.layer(1.46 - 0.1j, 200, i)

    R = solvePropagation(a)[0]
    r.append(np.abs(R) ** 2)
    R = solvePropagation(b)[0]
    r1.append(np.abs(R) ** 2)
    R = solvePropagation(c)[0]
    r2.append(np.abs(R) ** 2)
    R = solvePropagation(d)[0]
    r3.append(np.abs(R) ** 2)

plt.plot(wavelengths, r)
plt.plot(wavelengths, r1)
plt.plot(wavelengths, r2)
plt.plot(wavelengths, r3)
plt.legend(['1.46', '1.46-0.001j', '1.46-0.01j', '1.46-0.1j'], loc='best')
plt.show()
