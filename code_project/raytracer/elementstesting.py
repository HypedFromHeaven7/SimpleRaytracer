import numpy as np
from raytracer.mymath import vectors, semicircle, circle, normalise
from raytracer.rays import Ray
from raytracer.physics import refract
import matplotlib.pyplot as plt

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
            self.__aperature = aperture
            self.__curvature = curvature
            self.__n_1 = n_1
            self.__n_2 = n_2

        def z_0(self):
            return self.__z_0
        
        def aperture(self):
            return self.__aperature
        
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
            if discriminant < 0.0:
                return None
            
            if rad_curv > 0:
                l = - k_dot_r - np.sqrt(discriminant)

            else:
                l = - k_dot_r + np.sqrt(discriminant)

            Q = (ray_dir.normalise() * l + ray_pos).array()


            P = ray_pos.array()
            O = circlecent.array()
            # create circle
        
            y = np.linspace(start=-rad_curv, stop=rad_curv, num=40)
            z = semicircle(y, 0., O[2], R=rad_curv)

            plt.plot(Q[2], Q[1], 'X', color="Black")
            plt.text(Q[2], Q[1] + 1, 'Q', color="Black")
            plt.plot(P[2], P[1], 'X', color="Blue")
            plt.text(P[2], P[1] + 1, "P", color="Blue")
            plt.plot(O[2], O[1], '.', color='Green', label='O')
            plt.text(O[2], O[1] + 1, "O", color="Green")
            plt.plot(z, y, color="Red")
            plt.show()

            bzz = Q - circlecent.array()
            x_y_radius = np.sqrt(rad_curv * rad_curv - bzz[2] * bzz[2])

            x_y_plot = circle(self.__aperature, 500)
            plt.plot(x_y_plot[0], x_y_plot[1], color="Red")
            plt.plot(O[0], O[1], '.', color='Green', label='O')
            plt.text(O[0], O[1] + 1, "O", color="Green")
            plt.plot(Q[0], Q[1], 'X', color='Black')
            plt.text(Q[0], Q[1] + 1, 'Q', color="Black")
            plt.show()
            if self.__aperature < x_y_radius:
                return f"Fasle since the distance = {x_y_radius} while aperature = {self.__aperature}\nQ = {Q}"
            
            inters = f"Fasle since the distance = {x_y_radius}"
            return Q
            
        def propagate_ray(self, ray):
            intersect = self.intercept(ray=ray)
            if intersect is None:
                return None
            
            circent = np.array([0, 0, self.__z_0 + 1/self.__curvature])
            normal = intersect - circent
            new_direction = refract(ray.direc(), normal=normal, n_1=self.__n_1, n_2=self.__n_2)
            ray.append(pos=intersect, direc=new_direction)
            return ray#, intersect, circent, new_direction, normal

ray1 = Ray(pos=[7., 0., 0.])
ray2 = Ray(pos=[-7., 0., 0.])
ray3 = Ray(pos=[10., 0., 0.])
ray4 = Ray(pos=[-10., 0., 0.])
sr = SphericalRefraction(z_0=10, curvature=0.02, n_1=1., n_2=1.5, aperture=9.)

inter1 = sr.intercept(ray1)
print(inter1)

# prop = sr.propagate_ray(ray1)
# print(prop)

inter2 = sr.intercept(ray2)
print(inter2)
inter3 = sr.intercept(ray3)
print(inter3)
inter4 = sr.intercept(ray4)
print(inter4)

