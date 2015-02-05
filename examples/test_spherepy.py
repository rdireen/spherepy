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
import numpy as np
import plot_sphere
import profile


c = sp.zeros_coefs(48,48)
c[4,3] = 1.0
c[3,0] = 1.0

p = sp.ispht(c,150,150)
#profile.run('p = sp.ispht(c,602,602)')

c2 = sp.spht(p,48,48)

#profile.run('c2 = sp.spht(p,200,200)')

print(sp.L2_coef(c - c2)/sp.L2_coef(c))

T = np.abs(p.array)

plot_sphere.plot_mag_on_sphere(T)


a = raw_input()
