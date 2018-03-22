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

p = 20
n1 = 2.5
n2 = 1.5
# Create multilayer structure
n = np.tile(np.array((n1, n2)), p)
n = np.insert(n, 0, 1.0)
n = np.append(n, 1.0)
d = np.tile(np.array((800.0/4./n1, 800.0/4.0/n2)), p)
d = np.insert(d, 0, 0.0)
d = np.append(d, 0.0)
l = (673., 800.)
x = np.linspace(-1000., 5000., 4001)

a = MultiLayer(n=n, d=d, wvl=l)

# TE
E = a.field_TE(x)
TE = np.abs(E**2)

# TM
E = a.field_TM(x)
TM = np.abs(E**2)

legend = []
for i, wvl in enumerate(l):
    plt.plot(x, TE[i])
    legend += ['TE(%.1fnm)' % (wvl)]
    plt.plot(x, TM[i])
    legend += ['TM(%.1fnm)' % (wvl)]

plt.xlabel("Depth, nm")
plt.ylabel("|E|^2/|Eo|^2")
plt.title("Electric field of Bragg mirror")
plt.legend(legend, loc='best')
plt.show(block=True)
