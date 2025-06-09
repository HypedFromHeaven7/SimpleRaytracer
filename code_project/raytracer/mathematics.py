import numpy as np


def semicircle(R, x, c, d):
    if R > 0.:
        y = - np.sqrt(R * R - (x - c) * (x - c)) + d
        return y
    
    else:
        y = + np.sqrt(R * R - (x - c) * (x - c)) + d
        return y

def circle(R, num):
    x1 = np.linspace(start=-R, stop=R, num=num)
    y1 = semicircle(x=x1, c=0., d=0., R=R)

    x2 = np.linspace(start=R, stop=-R, num=num)
    y2 = semicircle(x=x2, c=0., d=0., R=-R)

    x = np.append(x1, x2)
    y = np.append(y1, y2)

    return x, y

def normalise(vector):
    vector = np.array(vector)
    mult = np.sum(vector * vector)
    return vector/np.sqrt(mult)

def neo_semicircle(num=10, radius=1., z_in=0., aperture=1.):
    if aperture >= np.abs(radius):
        aperture = radius

    y = np.linspace(start=-aperture, stop=aperture, num=num)
    if radius > 0:
        z = - np.sqrt(radius * radius - y * y) + z_in
        return z, y

    z = np.sqrt(radius * radius - y * y) + z_in
    # test print
    # print(f"The aperture of the lens = {aperture}")
    return z, y


def x_y_radius(radius=1., z=0., z_0=0.):
    act_z = z - z_0
    squared_z = act_z * act_z
    squared_rad = radius * radius
    square = squared_rad - squared_z

    if np.less(np.abs(square),np.abs(radius)*1e-9):
        square = 0.0
    return np.sqrt(square)

def rtrings(rmax, nrings, multi,z_0=0):
    radis=rmax/nrings
    yield (np.float64(0.0),np.float64(0.0), np.float64(0.0))
    for bn in range(0,nrings):
        rad=radis*(bn+1)
        arang=2*np.pi/multi/(bn+1)
        for tn in range(0,multi*(bn+1)):
            thet=arang*(tn)
            x = rad * np.cos(thet)
            y = rad * np.sin(thet)
            pos = np.array([x, y, z_0])
            yield pos


    

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