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

from tmmnlay import MultiLayer

aoi = np.linspace(0, 89.9, 1000)
TE = []
TM = []
# Wavelength is not important if both layers are seminfinite
a = MultiLayer(n=(2.0, 1.0), d=(0.0, 0.0), wvl=1)
for a.aoi in aoi:
    # TE
    r, t = a.rt_TE
    TE.append(r.real*r.real + r.imag*r.imag)

    # TM
    r, t = a.rt_TM
    TM.append(r.real*r.real + r.imag*r.imag)


plt.plot(aoi, TE)
plt.plot(aoi, TM)
plt.xlabel("Angle, deg")
plt.ylabel("Reflectance")
plt.title("Angle dependence of reflectivity")
plt.legend(['TE', 'TM'], loc='best')
plt.show(block=True)
