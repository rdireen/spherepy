import numpy as np
import spherepy as sp
from mayavi import mlab


def plot_mag_on_sphere(T):

    (nrows, ncols) = T.shape

  
    phi, theta = np.meshgrid(np.linspace(0, 2 * np.pi, ncols),
                             np.linspace(0, np.pi, nrows))
    X = T * np.cos(phi) * np.sin(theta)
    Y = T * np.sin(phi) * np.sin(theta)
    Z = T * np.cos(theta)
 
    s = mlab.mesh(X, Y, Z,  colormap='prism')
    mlab.show()


c = sp.zeros_coefs(10, 10)
c[5,1] = 1.0


p = sp.ispht(c, 100, 100)

T = np.abs(p.array)

plot_mag_on_sphere(T)
