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

import spherepy as sp
import profile

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

c = sp.random_coefs(500, 498)

p = sp.ispht(c, 506, 1010)
# profile.run('p = sp.ispht(c,602,602)')

c2 = sp.spht(p, 500, 498)

#profile.run('c2 = sp.spht(p,400,400)')

print(sp.L2_coef(c - c2) / sp.L2_coef(c))

#sp.plot_sphere.plot_mag_on_sphere(T)


