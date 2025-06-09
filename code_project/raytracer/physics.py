"""
A module to describe the physics of light.
"""

from raytracer.mymath import normalise
import numpy as np
import matplotlib.pyplot as plt

def refract(direc, normal, n_1, n_2):
    k_hat = normalise(direc)
    n_hat = normalise(normal)

    # factor between refractive indices
    mu = n_1/n_2
    
    #n_perp = normalise(n_perp)


    ndotk = np.dot(k_hat, n_hat)
    sin_theta1 = np.sqrt(1 - ndotk * ndotk)

    if sin_theta1 > n_2/n_1:
        return None
    
    # Calculating the refraction angles:
    trans_perpendicular = mu * (k_hat - ndotk * n_hat)
    trans_parallel = - np.sqrt(1 - mu * mu * (1 - ndotk * ndotk)) * n_hat

    transmitted_ray = trans_parallel + trans_perpendicular
    return transmitted_ray


def reflect(direc, normal):
    k1_hat = normalise(direc)
    n_hat = normalise(normal)
    k1dotn = np.dot(k1_hat, n_hat)
    return k1_hat - 2 * k1dotn * n_hat
# print(vectors([1,2,2]).normalise())
# print(normalise([1,2,2]))

# print(refract([1,1,1], [0,1,0], 1, 1))

# direc = np.array([0., 0., 1.])
# norm_lower = np.array([0., -1., -1.])
# a = refract(direc=direc, normal=norm_lower, n_1=1.0, n_2=1.5)
# b = refract(direc=np.array([0., 0.05, 1.]), normal=np.array([0., -1.9, -1.]), n_1=1.0, n_2=1.5)
# print(f"{b}\n\n\n")
# print(f"{a}\n\n\n{direc}\n{norm_lower}")