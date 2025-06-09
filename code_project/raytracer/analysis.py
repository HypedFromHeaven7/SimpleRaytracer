"""Analysis module."""
import matplotlib.pyplot as plt
from raytracer.rays import Ray, RayBundle
from raytracer.elements import SphericalRefraction, SphericalReflection, OutputPlane
from raytracer.plottingfunc import sphericlensshape, sphericplotter, rayplotter
from raytracer.lenses import PlanoConvex, BiConvex
import numpy as np
import scipy.optimize as sp

# testbundle = RayBundle()
great_sr = SphericalRefraction(100., 34., 0.03, 1., 1.5)
great_sr_focal_point = great_sr.focal_point()
greatlenses = [great_sr, OutputPlane(great_sr_focal_point)]
great_sr_focal_length = great_sr.focal_length()



def task8():
    ray1 = Ray(pos=[0., 0., -10.])
    ray2 = Ray(pos=[100., 0., -10.])
    ray3 = Ray(pos=[20., 10., -10.])

    sr = SphericalRefraction(z_0=10, curvature=0.02, n_1=1., n_2=1.5, aperture=50.)
    sr.propagate_ray(ray1)
    sr.propagate_ray(ray3)
    assert np.allclose(ray1.pos(), np.array([0.,0.,10.]))
    return True
    # assert sr.propagate_ray(ray2) is None
 

def task10():
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    sr = SphericalRefraction(100., 34., 0.03, 1., 1.5)
    lens = sphericlensshape(sr, 400)
    outplane = OutputPlane(250)
    
    rays = [Ray([0.,4.,0.]), 
            Ray([0.,1.,0.]),
            Ray([0.,0.2,0.]),
            Ray([0, 0, 0]),
            Ray([0, -0.2, 0]),
            Ray([0, -1, 0]),
            Ray([0, -4, 0])]
    figure_task10 = plt.figure()
    new_plot = figure_task10.add_subplot()
    for i in range(len(rays)):
        sr.propagate_ray(rays[i])
        outplane.propagate_ray(rays[i])
        points=rays[i].vertices()
        z = []
        y = []
        for n in range(len(points)):
            z.append(points[n][2])
            y.append(points[n][1])
        
        new_plot.plot(z, y)
    new_plot.plot(lens[0], lens[1], color='black')
    return figure_task10


def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """
    rays = [Ray([0.1, 0.1, 0]),
            Ray([0, 0, 0]),
            Ray([-0.1, -0.1, 0]),]

    
    focalpoint = great_sr.focal_point()
    figure_task11 = plt.figure()
    paraxials = figure_task11.add_subplot()
    outplane = OutputPlane(focalpoint)

    sphericplotter(rays, great_sr, outplane, paraxials, 2)
    
    return figure_task11, focalpoint


def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.
    """
    up = RayBundle(30, 15, 6)
    lenses = [great_sr, OutputPlane(great_sr_focal_point)]
    up.propagate_bundle(lenses)

    return up.track_plot()


def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """
    t13 = RayBundle()
    t13.propagate_bundle(greatlenses)
    return t13.spot_plot(), t13.rms()


def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """
    number = 100
    r_values = np.linspace(start=0.1, stop=10, num=number)
    bundle = [RayBundle(rmax=value).propagate_bundle(greatlenses) for value in r_values]
    ray25 = RayBundle(2.5).propagate_bundle(greatlenses)
    rms25 = ray25.rms()
    diffraction25 = 588e-6 * great_sr_focal_length / 2.5 / 2 # change diffraction function too
    rms_vals = list([r.rms() for r in bundle])
    
    diffraction = 588e-6 * great_sr_focal_length / r_values / 2
    figure_task14 = plt.figure()
    radius_plots = figure_task14.add_subplot()
    radius_plots.plot(r_values,diffraction, color="Blue")
    radius_plots.plot(r_values, rms_vals, color="Red")
    return figure_task14, rms25, diffraction25


def task15():
    """
    Task 15.

    In this function you will create plano-convex lenses in each orientation and propagate a RayBundle
    through each to their respective focal point. You should then plot the spot plot for each orientation.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spot plots and rms for plano-convex and convex-plano.
    """
    
    pc_lens = PlanoConvex(100, -0.02, 1.5168, 1., thickness=5.,aperture=50.)
    cp_lens = PlanoConvex(100, 0.02, 1.5168, 1., thickness=5.,aperture=50.)
    lenses = [pc_lens, cp_lens]
    
    
    beams = [RayBundle(), RayBundle()]

    foc = []
    for i in range(0,2):
        beam = beams[i]
        lens = lenses[i]

        focal_point = lens.focal_point()
        foc.append(focal_point)
        outplane = OutputPlane(focal_point)

        beam.propagate_bundle([lens, outplane])
        
    return beams[0].spot_plot(), foc[0], beams[1].spot_plot(), foc[1]


def task16():
    """
    Task 16.

    In this function you will be again plotting the radial dependence of the RMS and diffraction values
    for each orientation of your lens.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """
    pc_lens = PlanoConvex(100, -0.02, 1.5168, 1., thickness=5.,aperture=50.)
    cp_lens = PlanoConvex(100, 0.02, 1.5168, 1., thickness=5.,aperture=50.)
    lenses = [pc_lens, cp_lens]
    number = 10
    r_values = np.linspace(start=0.1, stop=10, num=number)
    rms35 =[]
    focus = cp_lens.focal_length()
    rms_vals = []
    
    # end = OutputPlane(focus)

    diffraction = 588e-6 * focus / r_values / 2
    diff35 = (588e-6 * focus / 3.5 / 2)
    
    for i in range(0,2):
        ray35 = RayBundle(3.5)
        end = OutputPlane(lenses[i].focal_point())
        bundle = [RayBundle(rmax=value).propagate_bundle([lenses[i], end]) for value in r_values]
        ray35.propagate_bundle([lenses[i], end])
        rms35.append(ray35.rms())
         # change diffraction function too
    
        rms_vals.append(list([r.rms() for r in bundle]))

    figure_task14 = plt.figure()
    radius_plots = figure_task14.add_subplot()
    radius_plots.plot(r_values,diffraction, color="Blue")
    radius_plots.plot(r_values, rms_vals[0], color="Red")
    radius_plots.plot(r_values, rms_vals[1], color="Green")
    return figure_task14, rms35[0], rms35[1], diff35


def task17():
    """cpl = PlanoConvex(z_0=100., curvature1=0.02, curvature2=0., n_inside=1.5168, n_outside=1., thickness=5., aperture=50.) will be first plotting the spot plot for your PlanoConvex lens with the curved
    side first (at the focal point). You will then be optimising the curvatures of a BiConvex lens
    in order to minimise the RMS spot size at the same focal point. This function should return
    the following items as a tuple in the following order:
    1. The comparison spot plot for both PlanoConvex (curved side first) and BiConvex lenses at PlanoConvex focal point.
    2. The RMS spot size for the PlanoConvex lens at focal point
    3. the RMS spot size for the BiConvex lens at PlanoConvex focal point

    Returns:
        tuple[Figure, float, float]: The combined spot plot, RMS for the PC lens, RMS for the BiConvex lens
    """
    bi = RayBundle()
    cp = RayBundle()
    
    cpl = PlanoConvex(z_0=100., curvature=0.02, n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    cp_foc = cpl.focal_point()
    end = OutputPlane(cp_foc)
    cp.propagate_bundle([cpl, end])

    def rms_calculator(cur):
        cur1, cur2 = cur
        bund = RayBundle()
        lens = BiConvex(z_0=100., curvature1=cur1, curvature2=cur2, n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
        bund.propagate_bundle([lens, end])
        return bund.rms()


    lop = sp.minimize(rms_calculator, [0.0146, -0.00563], options={'maxiter': 100000000000000000000, 'disp': False})
    kursk = lop.x

    ideal_lens = BiConvex(z_0=100., curvature1=kursk[0], curvature2=kursk[1], n_inside=1.5168, n_outside=1., thickness=5., aperture=50.)
    bi.propagate_bundle([ideal_lens, end])

    task17_fig = plt.figure()
    comparison = task17_fig.add_subplot()

    bi.spot_plot(comparison)
    cp.spot_plot(comparison)
    
    return task17_fig, cp.rms(), bi.rms()




def task18():
    """
    Task 18.

    In this function you will be testing your reflection modelling. Create a new SphericalReflecting surface
    and trace a RayBundle through it to the OutputPlane.This function should return
    the following items as a tuple in the following order:
    1. The track plot showing reflecting ray bundle off SphericalReflection surface.
    2. The focal point of the SphericalReflection surface.

    Returns:
        tuple[Figure, float]: The track plot and the focal point.

    """
    task18_fig = plt.figure()
    plot = task18_fig.add_subplot()

    reflector = SphericalReflection(aperture=6.)
    beam = RayBundle()

    final = OutputPlane(50)
    rinal = OutputPlane(200)

    beam.propagate_bundle([reflector, final])#, rinal])
    reflector.plotlens(plot)
    beam.track_plot(plot)
    return task18_fig, reflector.focalpoint()


if __name__ == "__main__":

    # Run task 8 function
    # task8()

    # Run task 10 function
    # FIG10 = task10()

    # Run task 11 function
    # FIG11, FOCAL_POINT = task11()

    # Run task 12 function
    # FIG12 = task12()

    # Run task 13 function
    # FIG13, TASK13_RMS = task13()

    # Run task 14 function
    # FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    # FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP = task15()

    # Run task 16 function
    # FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    # Run task 17 function
    # FIG17, CP_RMS, BICONVEX_RMS = task17()

    # Run task 18 function
    # FIG18, FOCAL_POINT = task18()

    plt.show()

