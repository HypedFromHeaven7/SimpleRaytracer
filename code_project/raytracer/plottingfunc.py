from raytracer.mymath import neo_semicircle
import numpy as np
import matplotlib.pyplot as plt


def sphericlensshape(sr, no_points=10, smaller_lens=None):
    rad = 1 / sr.curvature()
    
    aperture = sr.aperture()
    if smaller_lens is None or smaller_lens > aperture:
        pass
    else:
        aperture = smaller_lens
    center = rad + sr.z_0()
    z = neo_semicircle(no_points, rad, center, aperture)
    return z

def planeshape(plane, no_points=10):
    intercept = plane.z_0()
    aperture = plane.aperture()
    y = np.linspace(-aperture, aperture, no_points)
    z = np.full(no_points, intercept)
    return z, y

#probably good idea to remove outputplane call within fucnction
def sphericplotter(rays, sr, outplane, plotera, include_lens=1, no_points=40):
    ymaxx = 0.0
    for igoris in range(0, len(rays)):
        
        sr.propagate_ray(rays[igoris])
        outplane.propagate_ray(rays[igoris])
        points=rays[igoris].vertices()
        # x = np.fromiter((pos[0] for pos in points), dtype=float, count=len(points))
        y = np.fromiter((pos[1] for pos in points), dtype=float, count=len(points))
        z = np.fromiter((pos[2] for pos in points), dtype=float, count=len(points))
        plotera.plot(z, y)
        if np.greater(np.max(np.abs(y)), ymaxx):
            ymaxx = np.max(np.abs(y))

    if include_lens == 1:
        lens = sphericlensshape(no_points=no_points, sr=sr)
        plotera.plot(lens[0], lens[1], color='Black')

    if include_lens == 2:
        lens = sphericlensshape(no_points=no_points, sr=sr, smaller_lens=1.2*ymaxx)
        plotera.plot(lens[0], lens[1], color='Black')
    return #plotera

def rayplotter(rays,plotera):
    # print(f"number of rays is {len(rays)}")
    for igoris in range(0, len(rays)):
        
        points=rays[igoris].vertices()
        # x = np.fromiter((pos[0] for pos in points), dtype=float, count=len(points))
        y = np.fromiter((pos[1] for pos in points), dtype=float, count=len(points))
        z = np.fromiter((pos[2] for pos in points), dtype=float, count=len(points))
        plotera.plot(z, y)
        # print(f"Ray {igoris + 1} plotted, with y values {y}")
    return # plotera




# sr1 = SphericalRefraction(z_0=130, curvature=0.03, n_1=1.0, n_2=1.5, aperture=344)
# sr2 = SphericalRefraction(z_0=110, curvature=0.03, n_1=1.0, n_2=1.5, aperture=15)

# lens1 = sphericlensshape(40,sr1)
# lens2 = sphericlensshape(40,sr2)

# plt.plot(lens1[0],lens1[1], color='Blue')
# plt.plot(lens2[0],lens2[1], color='Black')
# sr = SphericalRefraction(z_0=100, curvature=-0.03, n_1=1.0, n_2=1.5, aperture=3.4)
# lens = sphericlensshape(40,sr)
# plt.plot(lens[0],lens[1], color='Red')
# plt.show()

# # def plotter(rays, )