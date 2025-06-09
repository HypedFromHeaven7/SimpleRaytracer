"""
This is for defining different types of lenses
"""
from raytracer.elements import SphericalRefraction, SinglePlane
from raytracer.plottingfunc import sphericlensshape, planeshape
from raytracer.mymath import x_y_rad_sq
import numpy as np

class PlanoConvex():
    """
    A lens with a convex shape on one side and a planer shape on the other.
    """
    def __init__(self, z_0, curvature, n_inside, n_outside, thickness, aperture):
        self.__radius = 1. / curvature
        self.__radmag = np.abs(self.__radius)
        self.__z_0= z_0
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__aperture = aperture
        self.__thickness = thickness
        if (self.__thickness <= self.__radmag) & (self.__aperture * self.__aperture > x_y_rad_sq(self.__radmag, self.__radmag - self.__thickness, 0.0)):
                self.__aperture = np.sqrt(x_y_rad_sq(self.__radmag, self.__radmag - self.__thickness, 0.0))
        
        if curvature > 0.:
            self.__tip = z_0
            self.__back = z_0 + thickness 
            self.__n_1 = n_outside
            self.__n_2 = n_inside
            
        
        else:
            self.__tip = z_0 + thickness
            self.__back = z_0
            self.__n_2 = n_outside
            self.__n_1 = n_inside
           

        self.__spheric = SphericalRefraction(self.__tip, self.__aperture, curvature, self.__n_1, self.__n_2)
        self.__plane = SinglePlane(self.__back,aperture=self.__aperture, n_1=self.__n_2, n_2=self.__n_1)

    def propagate_ray(self, ray):
        if self.__radius>0:
            self.__spheric.propagate_ray(ray)
            self.__plane.propagate_ray(ray)
        else:
            self.__plane.propagate_ray(ray)
            self.__spheric.propagate_ray(ray)
        return ray
    
    def lensshape(self, figure, num=10):
        spher = sphericlensshape(self.__spheric, num)
        plane = planeshape(self.__plane, num)
        zvals = np.append(spher[0], plane[0])
        yvals = np.append(spher[1], plane[1])
        figure.plot(zvals, yvals)
        return zvals, yvals
        
    def focal_length(self):
        f = 1 / (self.__n_inside/self.__n_outside - 1) * np.abs(self.__radius)
        return f

    
    def focal_point(self):

        if self.__radius > 0:
            addition =- self.__thickness / self.__n_inside
            return self.__z_0 + self.focal_length() + self.__thickness + addition
        return self.__z_0 + self.focal_length() + self.__thickness


class BiConvex():
    """
    A class describing biconvex lenses
    """
    def __init__(self, z_0, curvature1, curvature2, n_inside, n_outside, thickness, aperture):
        self.__n_in = n_inside
        self.__n_out = n_outside
        self.__radius1 = 1 / curvature1
        self.__radius2 = 1 / curvature2
        self.__thickness = thickness
        self.__z_0 = z_0
        self.__cur1 = curvature1
        self.__cur2 = curvature2
        self.__spher1 = SphericalRefraction(z_0, aperture, curvature1, self.__n_out, self.__n_in)
        self.__spher2 = SphericalRefraction(z_0+thickness, aperture, curvature2, self.__n_in, self.__n_out)


    def focal_length(self):
        f = 1 / (self.__n_in/self.__n_out - 1) / (self.__cur1 - self.__cur2)
        #testing
        # print("test edits in lenses biconvex")
        # print(f"cur1 = {self.__cur1}, cur2 = {self.__cur2}, n = {self.__n_in/self.__n_out}")
        #testing
        return f
    
    def focal_point(self):
        addition =- self.__thickness / self.__n_in
        return self.__z_0 + self.focal_length() + self.__thickness + addition
    
    def propagate_ray(self, ray):
        self.__spher1.propagate_ray(ray)
        self.__spher2.propagate_ray(ray)
        return ray

    def lensshape(self, figure, num=10):
        spher1 = sphericlensshape(self.__spher1, num)
        spher2 = sphericlensshape(self.__spher2, num)
        zvals = np.append(spher1[0], spher2[0])
        yvals = np.append(spher1[1], spher2[1])
        figure.plot(zvals, yvals)
        return zvals, yvals
    

