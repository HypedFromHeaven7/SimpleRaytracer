"""
Module to make ray objects.
"""
import numpy as np
import matplotlib.pyplot as plt
# from raytracer.mymath import rtrings
from raytracer.genpolar import rtrings, position
from raytracer.plottingfunc import rayplotter

class Ray:
    """
    A superclass defining rays. 
    
    Arguments:
    3D position.
    3D direction.


    methods:
    pos() returns current point of ray
    direc() returns direction of ray
    """
    
    def __init__(self, pos=[0., 0., 0.], direc=[0., 0., 1.]):
        self.__pos = np.array(pos)
        self.__direc = np.array(direc)
        self.__vertices = [self.__pos]


    def pos(self):
        return self.__pos
    
    def direc(self):
        return self.__direc
    
    def append(self, pos, direc=None):
        if pos is not None:
            self.__pos = np.array(pos)
            self.__vertices.append(np.array(pos.copy()))
        
        if direc is not None:
            self.__direc = np.array(direc)

        return self
    
    def vertices(self):
        return self.__vertices
    
    def __str__(self):
        return str(self.__pos) + str(self.__direc)
    
    def __repr__(self):
        return f"r = {self.__pos} k = {self.__direc}"
    

    

class RayBundle():
    """
    A class defining circular rings of rays,
    
    Inputs:
    maximum radius,
    number of rings,
    multitude of points,
    """
    def __init__(self, rmax=5., nrings=5, multi=6):
        self.__rmax=rmax
        self.__generator = rtrings(rmax, nrings, multi)
        self.__rays = []
        for i, (rad, thet) in enumerate(self.__generator):
            self.__rays.append(Ray(pos=position(rad, thet)))

    def propagate_bundle(self, optics):
        for ray_num in range(0, len(self.__rays)):
            for opt_num in range(0, len(optics)):
                optics[opt_num].propagate_ray(self.__rays[ray_num])
        return self

    def track_plot(self, test_figure = None):
        if test_figure is None:
            test_figure = plt.figure()
            posture = test_figure.add_subplot()
            rayplotter(self.__rays, posture)
            return test_figure
        rayplotter(self.__rays, test_figure)
        return test_figure
    
        

    def test(self):
        return(self.__rays)
    
    def rms(self):
        # method to calculate the RMS size of the bundle
        mean_pos = np.sum([positions.pos() for positions in self.__rays], axis=0)/len(self.__rays)
        squared_disfromcenter = np.sum([np.square(positions.pos() - mean_pos) for positions in self.__rays], axis=0)/len(self.__rays)
        rootmeansquared = np.sqrt((squared_disfromcenter[0] + squared_disfromcenter[1])) #np.sum(squared_disfromcenter))
        # print(mean_pos, squared_disfromcenter, rootmeansquared)
        return rootmeansquared

        
    def spot_plot(self, x_y_plot = None):
        coordinates = np.array([positions.pos() for positions in self.__rays])
        if x_y_plot is None:
            plane_look = plt.figure()
            x_y_plot = plane_look.add_subplot()
            x_y_plot.scatter(coordinates[:, 1], coordinates[:, 0])
            return plane_look
        x_y_plot.scatter(coordinates[:, 1], coordinates[:, 0])
    
    def __repr__(self):
        return f"Max radius = {self.__rmax}"


# up.test()
# generator = rtrings(5,5,6)
# rays = []
# for i, (rad, thet) in enumerate(generator):
#     print(rad, thet, i)
    #rays.append(Ray(pos=position(rad, thet)))
    