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

__version__ = '1.0.0'
__title__ = 'Cython implementarion of the Transfer Matrix Method'
__mod__ = 'python-tmmnlay'
__author__ = 'Ovidio Peña Rodríguez'
__email__ = 'ovidio@bytesfall.com'
__url__ = 'https://github.com/ovidiopr/tmmnlay'
__download_url__ = 'https://github.com/ovidiopr/tmmnlay'

from distutils.core import setup
from Cython.Build import cythonize

setup(name = __mod__,
      version = __version__,
      description = __title__,
      long_description="""This extension implements the Transfer Matrix Method, which is \
 used in optics and acoustics to analyze the propagation of \
 electromagnetic or acoustic waves through a multilayer medium.""",
      author = __author__,
      author_email = __email__,
      maintainer = __author__,
      maintainer_email = __email__,
      keywords = ['Transfer Matrix Method', 'Multilayered medium', 'Reflection', 'Transmission', 'Absorption'],
      url = __url__,
      download_url = __download_url__,
      classifiers=['Development Status :: 5 - Production/Stable',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Physics'
                  ],
      license = 'GPL',
      platforms = 'any',
      ext_modules = cythonize("tmmnlay.pyx")
)

