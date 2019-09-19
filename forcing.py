#file forcing.py
import math

def melt(x,y,t,thck,topg,gl_proximity=0.0,gl_proximity_scale=1.0):
    d = - math.log(gl_proximity + 1.0e-10) * gl_proximity_scale
    return - (d/1.0e+3)**2 * gl_proximity

def meltforce(x,y,t,thck,topg,gl_proximity=0.0,gl_proximity_scale=1.0):
    mmean = 0.0
    L = 1.0e5
    if t > 20.0 and t < 100.0:
        ft = t - 20.0
        mmean = - 0.1875 * ft
    elif t >= 100.0:
        mmean = - 15.0
    mo = ( mmean * L ) / gl_proximity_scale
    m = mo * gl_proximity
    return m

