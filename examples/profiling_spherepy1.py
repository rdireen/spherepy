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

Nmax = 400
Nrows = 2000

c = sp.random_coefs(Nmax, Nmax)


p = sp.ispht(c, Nrows, Nrows)
# profile.run('p = sp.ispht(c,602,602)')

#c2 = sp.spht(p, Nmax, Nmax)
c2 = None
profile.run('c2 = sp.spht(p,Nmax,Nmax)',sort=1)

print(sp.L2_coef(c - c2) / sp.L2_coef(c))


