import numpy as np
from raytracer.mymath import vectors, x_y_radius
from raytracer.physics import refract, reflect, normalise
from raytracer.plottingfunc import sphericlensshape

class OpticalElement:
    def intercept(self, ray):
            raise NotImplementedError('intercept() needs to be implemented in derived classes')
    def propagate_ray(self, ray):
            raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')
    


class SphericalRefraction(OpticalElement):
    """
    A subclass of the super class OpticalElement that describes spherical refraction.
    
    Inputs:
    z_0, the z position of the Element that is centered on z.
    aperature, the lens' radius or the maximum extent of the lens on the x-y plane.
    curvature, the curvature of the lens' spherical element.
    n_1, the refractive index before the lens.
    n_2, the refractive index in the lens.
    """
    def __init__(self, z_0=1., aperture=1., curvature=1., n_1=1., n_2=1.):
        self.__z_0 = z_0
        self.__aperture = aperture
        self.__curvature = curvature
        self.__n_1 = n_1
        self.__n_2 = n_2
        self.__circlcent = np.array([0, 0, self.__z_0 + 1/self.__curvature])

    def lensshape(self, num=10):
        spher = sphericlensshape(self, num)
        return spher[0], spher[1]

    def z_0(self):
        return self.__z_0
    
    def aperture(self):
        return self.__aperture
    
    def curvature(self):
        return self.__curvature
    
    def n_1(self):
        return self.__n_1
    
    def n_2(self):
        return self.__n_2


    def intercept(self, ray):
        """
        A method that return the distance of interception and states if there will be an interception.
        """
        #sets up variables
        rad_curv = 1/self.__curvature
        ray_pos = vectors(ray.pos())
        circlecent = vectors(r=[0, 0, self.__z_0 + rad_curv])
        ray_dir = vectors(ray.direc())

        r_dis = ray_pos - circlecent 
        k_dot_r = ray_dir.normalise() @ r_dis

        discriminant = k_dot_r * k_dot_r - (r_dis.mod_sq() - rad_curv * rad_curv)

        # imaginary solution indicates no interception
        if discriminant < 0.0:
            return None
        
        
        # depending on the curvature, concave/convex, we have different solutions to the quadratic.
        if rad_curv > 0:
            l = - k_dot_r - np.sqrt(discriminant)

        else:
            l = - k_dot_r + np.sqrt(discriminant)
            
        Q = (ray_dir.normalise() * l + ray_pos).array()
        
        
        xy_r = x_y_radius(radius=rad_curv, z=Q[2], z_0=circlecent.array()[2]) 
        if self.__aperture >= xy_r:
            return Q
        
        else:
            return None
        
    def propagate_ray(self, ray):
        Q_point = self.intercept(ray)
        if Q_point is None: # reflection
            return None
        
        circent = np.array([0, 0, self.__z_0 + 1/self.__curvature])
        if self.__curvature > 0:
            normal = Q_point - circent
        else:
            normal = circent - Q_point
        new_direction = refract(ray.direc(), normal=normal, n_1=self.__n_1, n_2=self.__n_2)

        ray.append(pos=Q_point, direc=new_direction)
        return ray
    
    def circlenter(self):
        return self.__circlcent
    

    # def focal_array(self):
    #     focal_length = self.__n_1 / self.__curvature / (self.__n_2 - self.__n_1)
    #     return self.__circlcent + np.array([0.,0.,focal_length])
    
    def focal_length(self):
        return self.__n_2 / self.__curvature / (self.__n_2 - self.__n_1)
    
    def focal_point(self):
        return self.__z_0 + self.focal_length()
    

    def __repr__(self):
        return f"z_0 = {self.__z_0}, curvature = {self.__curvature}, aperature = {self.__aperture}"
    


class SphericalReflection(OpticalElement):
    """
    A subclass defining spherical reflection.
    """
    def __init__(self, z_0 = 100, aperture = 6., curvature = -0.02):
        self.__z_0 = z_0
        self.__curvature = curvature
        self.__aperture = aperture 


    def intercept(self, ray):
        """
        A method that return the distance of interception and states if there will be an interception.
        """
        #sets up variables
        rad_curv = 1/self.__curvature
        ray_pos = vectors(ray.pos())
        circlecent = vectors(r=[0, 0, self.__z_0 + rad_curv])
        ray_dir = vectors(ray.direc())

        r_dis = ray_pos - circlecent 
        k_dot_r = ray_dir.normalise() @ r_dis

        discriminant = k_dot_r * k_dot_r - (r_dis.mod_sq() - rad_curv * rad_curv)

        # imaginary solution indicates no interception
        if discriminant < 0.0:
            return None
        
        
        # depending on the curvature, concave/convex, we have different solutions to the quadratic.
        if rad_curv > 0:
            l = - k_dot_r - np.sqrt(discriminant)

        else:
            l = - k_dot_r + np.sqrt(discriminant)
            
        Q = (ray_dir.normalise() * l + ray_pos).array()
        
        
        xy_r = x_y_radius(radius=rad_curv, z=Q[2], z_0=circlecent.array()[2]) 
        if self.__aperture >= xy_r:
            return Q
        
        else:
            return None
        
    def propagate_ray(self, ray):
        Q_point = self.intercept(ray)
        if Q_point is None: # reflection
            return None
        
        circent = np.array([0, 0, self.__z_0 + 1/self.__curvature])
        if self.__curvature > 0:
            normal = Q_point - circent
        else:
            normal = circent - Q_point
        new_direction = reflect(ray.direc(), normal=normal)

        ray.append(pos=Q_point, direc=new_direction)
        return ray
        
    def z_0(self):
        return self.__z_0
    
    def aperture(self):
        return self.__aperture
    
    def curvature(self):
        return self.__curvature
    
    def lensshape(self, num=10):
        spher = sphericlensshape(self, num)
        return spher[0], spher[1]
    
    def plotlens(self, plot):
        x, y = self.lensshape()
        plot.plot(x, y)
        return self
    
    def focallength(self):
        return 1 / (2 * self.__curvature)
    
    def focalpoint(self):
        return self.__z_0 + self.focallength()


class OutputPlane(OpticalElement):
    """
    A subclass of the super class OpticalElement that describes the point at which the ray terminates.

    inputs:
    the z at which ray's journey ends.
    """
    def __init__(self, z_0, refr=1):
        self.__z_0 = z_0
        self.__refr = refr

    def intercept(self, ray):
        """
        A method that return the point of interception and states if there will be an interception.
        """
        ray_pos = ray.pos()
        k_hat = normalise(ray.direc())
        kdotz = np.dot(k_hat, np.array([0., 0., 1.]))
        if kdotz == 0:
            return None
        cons = (self.__z_0 - ray_pos[2]) / k_hat[2]
        return ray_pos + k_hat * cons
    
    def propagate_ray(self, ray):
        ray.append(self.intercept(ray), ray.direc())
        return ray
    

class SinglePlane(OpticalElement):
    """
    A sub class describing a plane intersection.
    """
    def __init__(self, z_0, aperture, n_1, n_2):
        self.__z_0 = z_0
        self.__n_1 = n_1
        self.__n_2 = n_2
        self.__aperture = aperture
        self.__direc = np.array([0.,0., z_0])

    def intercept(self, ray):
        """
        A method that return the point of interception and states if there will be an interception.
        """
        ray_pos = ray.pos()
        k_hat = normalise(ray.direc())
        kdotz = np.dot(k_hat, self.__direc)
        if kdotz <= 0:
            return None
        cons = (self.__z_0 - ray_pos[2]) / k_hat[2]
        inter = ray_pos + k_hat * cons
        rad_dis = np.sqrt(inter[0] * inter[0] + inter[1] * inter[1])
        if rad_dis > self.__aperture:
            return None
        return inter
    
    def propagate_ray(self, ray):
        new_pos = self.intercept(ray)
        if new_pos is not None:
            new_direc = refract(ray.direc(), -self.__direc, self.__n_1, self.__n_2)
            ray.append(new_pos, new_direc)
        return ray
    
    def aperture(self):
        return self.__aperture
    
    def z_0(self):
        return self.__z_0







        
