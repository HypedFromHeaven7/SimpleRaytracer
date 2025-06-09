"""
Module to make math easier.
"""

import numpy as np
import matplotlib.pyplot as plt

class vectors:
    """
    A class of vector operations.
    """
 
    unit_x = np.array([1,0,0])
    def __init__(self, r=[0.,0.,0.]):
        self.__vect = np.array(r)

    def normalise(self, n=1):
        mult = np.sum(self.__vect * self.__vect)
        return vectors(self.__vect/np.sqrt(mult))    

    def mod_sq(self):
        return np.sum(self.__vect * self.__vect)
    
    def mod(self):
        return np.sqrt(np.sum(self.__vect * self.__vect))
    
    def see(self):
        print(self.__vect)

    def array(self):
        return self.__vect
    
    def x_y_radius(self):
        x_y = vectors(r=[self.__vect[0], self.__vect[1], 0.])
        return x_y.mod()

    ###################
    #### Operators ####
    ###################

    def __sub__(self, other):
        return vectors(self.__vect - other.__vect)
    
    def __add__(self, other):
        return vectors(self.__vect + other.__vect)
    
    #scalar multiplication
    def __mul__(self, float=1.0):
        return vectors(self.__vect * float)
    
    # dot product
    def __matmul__(self, other):
        return np.dot(self.__vect, other.__vect) 
    
    def __repr__(self):
        return self.__vect
    
    def __str__(self):
        array = self.__vect
        return str(array)
    

def semicircle(x, c, d, R):
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

def x_y_rad_sq(radius=1., z=0., z_0=0.):
    act_z = z - z_0
    squared_z = act_z * act_z
    squared_rad = radius * radius
    square = squared_rad - squared_z

    if np.less(np.abs(square),np.abs(radius)*1e-9):
        square = 0.0
    if square >= 0.0:
        return square
    
    return square





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



