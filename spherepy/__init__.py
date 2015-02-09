
import json
from os.path import dirname


with open(dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']
__author__ = _info['author']

from .spherepy import *
from .sbessel import *

try:
    from .plot_sphere import *
except:
    pass


__all__ = ['pysphi','csphi']
