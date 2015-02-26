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

"""***************************************************************************

               plot_sphere: Plotting Routines for SpherePy 

Randy Direen
2/11/2015

These routines are used for plotting results calculated with SpherePy. To use
them, matplotlib must be installed.

***************************************************************************"""

#---------------------------------------------------------------------Built-ins
from __future__ import division
import numpy as np
from functools import wraps

#TODO: Change all xrange instances to range
#and do a 'from six.moves import range' here
from six.moves import xrange

#------------------------------------------------------------------------Custom
import spherepy as sp

#==============================================================================
# Matplotlib Checking Decorator 
#==============================================================================

MPLINSTALLED = True
try:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
except:
    MPLINSTALLED = False

msg = "Plotting requires matplotlib. You can install matplotlib by typing " +\
      "'pip install matplotlib' at the command line. If that doesn't work," +\
      " google it."

def matplotlibensure(func):
    """If matplotlib isn't installed, this decorator alerts the user and 
    suggests how one might obtain the package."""  
    @wraps(func)
    def wrap(*args):
        if MPLINSTALLED == False:
            raise ImportError(msg)
        
        return func(*args)   
        
    return wrap

#==============================================================================
# Plotting Routines
#==============================================================================
    
@matplotlibensure
def simple_plot(n, m):

    c = sp.zeros_coefs(48, 48)
    c[n, m] = 1.0
    p = sp.ispht(c, 100, 100)
    T = np.abs(p.array)

    plot_mag_on_sphere(T)

@matplotlibensure
def plot_mag_on_sphere(T):

    (nrows, ncols) = T.shape

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    phi, theta = np.meshgrid(np.linspace(0, 2 * np.pi, ncols + 1),
                             np.linspace(0, np.pi, nrows))

    v = np.array(T[:,0]).reshape(-1, 1)
    T = np.hstack((T, v))

    X = T * np.cos(phi) * np.sin(theta)
    Y = T * np.sin(phi) * np.sin(theta)
    Z = T * np.cos(theta)
    # ax.plot_surface(X,Y, Z,rstride=1, cstride=1, cmap= 'jet',alpha=.5,
    #                linewidth=0.5)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet', alpha=.5,
                    linewidth=0.5)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X.max() - X.min(),
                          Y.max() - Y.min(),
                          Z.max() - Z.min()]).max()
    Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + \
         0.5 * (X.max() + X.min())
    Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + \
         0.5 * (Y.max() + Y.min())
    Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + \
         0.5 * (Z.max() + Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

@matplotlibensure
def plot_sphere_mag(patt):

    T = sp.mag(patt)

    nrows = patt.nrows
    ncols = patt.ncols

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    delta = 2 * np.pi / ncols + 1
    phi, theta = np.meshgrid(np.linspace(0, 2 * np.pi, ncols + 1),
                             np.linspace(0, np.pi, nrows))

    v = np.array(T[:,0]).reshape(-1, 1)
    T = np.hstack((T, v))

    X = T * np.cos(phi) * np.sin(theta)
    Y = T * np.sin(phi) * np.sin(theta)
    Z = T * np.cos(theta)

    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet', alpha=.5,
                    linewidth=0.5)



    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X.max() - X.min(),
                          Y.max() - Y.min(),
                          Z.max() - Z.min()]).max()
    Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + \
         0.5 * (X.max() + X.min())
    Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + \
         0.5 * (Y.max() + Y.min())
    Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + \
         0.5 * (Z.max() + Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

@matplotlibensure
def pcolor_coefs(coefs):

    A = coefs._array_2d_repr()
    z = np.abs(A)
    el = np.zeros((1, 2 * coefs.mmax + 1), dtype=np.float64)
    z = np.concatenate((z, el), axis=0)
    el = np.zeros((coefs.nmax + 2, 1), dtype=np.float64)
    z = np.concatenate((z, el), axis=1)
    z_min, z_max = z.min(), z.max()
    x, y = np.meshgrid(np.array(range(-coefs.mmax, coefs.mmax + 2)) - 0.5,
                      np.array(range(0, coefs.nmax + 2)) - 0.5)

    plt.pcolormesh(x, y, z, cmap='jet', vmin=z_min, vmax=z_max)
    plt.title('Spherical Coefficients')
    plt.xlabel('m')
    plt.ylabel('n')
    # set the limits of the plot to the limits of the data
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()

    plt.show()

@matplotlibensure
def plot_coefs(coefs):

    z = coefs._array_2d_repr()
    for n in xrange(0, coefs.nmax + 1):
        for m in xrange(1, coefs.mmax + 1):
            if np.abs(m) > n:
                z[n, m + coefs.mmax] = np.inf
                z[n, -m + coefs.mmax] = np.inf

    el = np.zeros((1, 2 * coefs.mmax + 1), dtype=np.float64)
    z = np.concatenate((z, el), axis=0)

    z = np.abs(z).T
    x = np.array(range(0, coefs.nmax + 2))

    for bs in z:
        plt.plot(x, bs)

    plt.xlabel('m')
    plt.show()




