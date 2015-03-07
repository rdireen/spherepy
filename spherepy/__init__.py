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
from __future__ import print_function
from __future__ import absolute_import

import json
import os
from os.path import dirname
import sys

import spherepy.pysphi as pysphi
import spherepy.csphi as csphi
import spherepy.file as file
import spherepy.verify as verify
from spherepy.spherepy import *
from spherepy.sbessel import *

with open(dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']
__author__ = _info['author']
__use_cext__ = _info['use_cext']


#Import matplotlib plotting if it has been installed
try:
    #Python3x
    from spherepy.plot_sphere import *
except ImportError:
    pass

