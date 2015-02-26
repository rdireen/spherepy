import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def random_blanket(mag, nrows, ncols, nrow_modes, ncol_modes):

    if (np.mod(nrows,2) == 1) or (np.mod(ncols,2) == 1):
        raise ValueError("rows and columns must be even")

    MaxPmodes = int(nrows / 2) - 1
    MaxZmodes = int(ncols / 2) - 1

    if (ncol_modes > MaxZmodes) or (nrow_modes > MaxPmodes):
        raise ValueError("Number of modes must be less than Length / 2")


    vec = np.random.normal(0.0, 1.0, nrows * ncols)
    data = vec.reshape((nrows, ncols))

    fdata = np.fft.fft2(data) / (nrows * ncols)

    fdata[:,ncol_modes:-ncol_modes + 1] = 0.0
    fdata[nrow_modes:-nrow_modes + 1,:] = 0.0

    filtered = np.fft.ifft2(fdata)

    return mag*filtered / np.max(filtered)


def cylr(len_z,radius,data):
    
    _ncirc = data.shape[0]
    _delta_t = 2*np.pi / _ncirc
    _delta_z = len_z / data.shape[1]
    N_z = data.shape[1]
 
    verts = []
    for n in range(0, N_z):
        for m in range(0, _ncirc):
            verts.append(((data[m,n] + radius)*np.cos(m*_delta_t),
                         (data[m,n] + radius)*np.sin(m*_delta_t) ,
                         _delta_z * n))
        

    faces = []
    for n in range(0, N_z - 1):
        for m in range(0, _ncirc - 1):
            faces.append((n*_ncirc + m,
                        n*_ncirc + m + 1,
                        (n+1)*_ncirc  + m + 1,
                        (n+1)*_ncirc + m ))
        faces.append((n * _ncirc + _ncirc - 1,
                    n * _ncirc,
                    n * _ncirc + _ncirc,
                    n * _ncirc + _ncirc + _ncirc - 1))
                    
    return (verts, faces)
                 

nrows = 16
ncols = 20

Z = np.abs(random_blanket(1.0,nrows, ncols, 3, 3))
(verts, faces) = cylr(20,2.0,Z)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(np.linspace(0,15, ncols),
                   np.linspace(0,2*np.pi,nrows))

ax.plot_surface(X, Y, Z,cmap= 'jet',
                rstride=1, cstride=1,)

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