import numpy as np

def rtrings(rmax=5.,nrings=5,multi=6,wing=0, z_0=0.):
    radis=rmax/nrings
    yield (0.0,0.0)
    
    for bn in range(0,nrings):
        rad=radis*(bn+1)
        arang=2*np.pi/multi/(bn+1)
        for tn in range(0,multi*(bn+1)):
            thet=arang*(tn)
            # if wing==0:
            yield rad, thet
            
            # else:
            #     x = rad * np.cos(thet)
            #     y = rad * np.sin(thet)
            #     pos = np.array([x, y, z_0])
            #     yield pos

def position(rad, thet, z_0=0.0):
    x = rad * np.cos(thet)
    y = rad * np.sin(thet)
    return np.array([x, y, z_0])



