# Copyright (C) 2015  Randy Direen <spherepy@direentech.com>
#
# This file is part of SpherePy.
#
# SpherePy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SpherePy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SpherePy.  If not, see <http://www.gnu.org/licenses/>

from __future__ import division

__all__ = ['pysphi']

import json
import os
from os.path import dirname
import sys

try:
    #Python27
    import pysphi
    import csphi
    import file
    import verify
except ImportError:
    #Python3x
    import spherepy.pysphi as pysphi
    import spherepy.csphi as csphi
    import spherepy.file as file
    import spherepy.verify as verify

with open(dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']
__author__ = _info['author']

try:
    #Python27
    from .spherepy import *
    from .sbessel import *
except ImportError:
    #Python3x
    from spherepy.spherepy import *
    from spherepy.sbessel import *

try:
    #Python27
    from .plot_sphere import *
except ImportError:
    try:
        #Python3x
        from spherepy.plot_sphere import *
    except ImportError:
        pass

