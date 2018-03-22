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
#    but WITHOUT ANY WARwavelengthTY; without even the implied warwavelengthty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt

from tmmnlay import MultiLayer

n1 = 1.5
n2 = np.sqrt(n1)
d = 700./n2/4.     # quarter-wavelength coating

wavelength = np.linspace(200, 1600, 1401)

# substrate layer
a = MultiLayer(n=(1.0, n1), d=(0.0, 0.0), wvl=wavelength)
r, t = a.rt_TE
refl0 = np.abs(r**2)

# antireflective layer layer "left" of substrate
b = MultiLayer(n=(1.0, n2, n1), d=(0.0, d, 0.0), wvl=wavelength)
r, t = b.rt_TE
refl = np.abs(r**2)

plt.plot(wavelength, refl0)
plt.plot(wavelength, refl)
plt.xlabel("Wavelength, nm")
plt.ylabel("Reflectance")
plt.title("Reflectance of ideal single-layer antireflective coating")
plt.legend(['Substrate', 'Coated substrate'], loc='best')
plt.show(block=True)
