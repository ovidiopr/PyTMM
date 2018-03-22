#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#
#    Copyright (C) 2018 Ovidio Peña Rodríguez <ovidio@bytesfall.com>
#
#    This file is part of tmmnlay
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

n = 1.5
d = 600  # slab thickness, nm
l = np.linspace(400., 1000., 601)  # wavelength, nm
aoi = np.linspace(0, 89.9, 1000)

a = MultiLayer(n=(1.0, n, 1.0), d=(0.0, d, 0.0), wvl=l)

legend = []
for a.aoi in (70., 75.):
    # TE
    rs, ts = a.rt_TE
    # TM
    rp, tp = a.rt_TM

    rho = rp/rs

    psi = np.arctan(np.abs(rho))*180./np.pi
    delta = np.angle(rho)*180./np.pi

    plt.plot(l, psi)
    legend += ['Psi(%.1f deg)' % (a.aoi)]
    plt.plot(l, delta)
    legend += ['Delta(%.1f deg)' % (a.aoi)]

plt.xlabel("Wavelength, nm")
plt.ylabel("Psi & Delta, deg")
plt.title("Ellepsometric spectra")
plt.legend(legend, loc='best')
plt.show(block=True)
