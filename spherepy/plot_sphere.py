from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

def plot_mag_on_sphere(T):

    (nrows,ncols) = T.shape

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    phi, theta = np.meshgrid(np.linspace(0,2*np.pi,ncols),
                             np.linspace(0,np.pi,nrows))
    X = T*np.cos(phi)*np.sin(theta)
    Y = T*np.sin(phi)*np.sin(theta)
    Z = T*np.cos(theta)
    #ax.plot_surface(X,Y, Z,rstride=1, cstride=1, cmap= 'jet',alpha=.5,linewidth=0.5)
    ax.plot_surface(X,Y, Z,rstride=1, cstride=1, cmap= 'jet',alpha=.5,linewidth=0.5)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X.max()-X.min(),
                          Y.max()-Y.min(),
                          Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + \
         0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + \
         0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + \
         0.5*(Z.max()+Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
       ax.plot([xb], [yb], [zb], 'w')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')


    plt.show()

def pcolor_coefs(coefs):

    A = coefs._array_2d_repr()
    z = np.abs(A)
    el = np.zeros((1,2*coefs.mmax + 1),dtype=np.float64)
    z = np.concatenate((z, el), axis=0)
    el = np.zeros((coefs.nmax + 2,1),dtype=np.float64)
    z = np.concatenate((z, el), axis=1)
    z_min, z_max = z.min(), z.max()
    x,y = np.meshgrid(np.array(range(-coefs.mmax,coefs.mmax+2)) - 0.5,
                      np.array(range(0,coefs.nmax + 2)) - 0.5)

    plt.pcolormesh(x, y, z, cmap= 'jet', vmin=z_min, vmax=z_max)
    plt.title('Spherical Coefficients')
    plt.xlabel('m')
    plt.ylabel('n')
    # set the limits of the plot to the limits of the data
    plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar()

    plt.show()

def plot_coefs(coefs):

    z = coefs._array_2d_repr()
    for n in xrange(0,coefs.nmax + 1):
        for m in xrange(1,coefs.mmax + 1):
            if np.abs(m) > n:
                z[n,m+coefs.mmax] = np.inf
                z[n,-m+coefs.mmax] = np.inf

    el = np.zeros((1,2*coefs.mmax + 1),dtype=np.float64)
    z = np.concatenate((z, el), axis=0)

    z = np.abs(z).T
    x = np.array(range(0,coefs.nmax + 2))

    for bs in z:
        plt.plot(x,bs)

    plt.xlabel('m')
    plt.show()




