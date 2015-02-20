
import json
from os.path import dirname
import pysphi
import csphi
import file
import verify


with open(dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']
__author__ = _info['author']

from .spherepy import *
from .sbessel import *

#Set to true if you want to use the c extensions
use_cext = True    

try:
    from .plot_sphere import *
except:
    pass

