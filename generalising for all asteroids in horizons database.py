# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 03:23:19 2020

@author: shari
"""

import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from astroquery.jplhorizons import Horizons

nhao = {"lon" : 52.4862,"lat" : 1.8904,"elevation" : 140}
AU2M = 149597870700
v_0 = -26.762
sph = 3600
mpkm = 1000
mc = -55.87
I = []
sun_dist = []
Type =[]
d = 500
n = 1
while n <= d:
    asteroid = Horizons(n, epochs={'start':'2019-01-01', 'stop':'2019-02-01', 'step':'1d'})
    eph = asteroid.ephemerides()
    H = eph["H"][0]
    sun_distance = np.average(eph["r"])
    obs_dist = np.average(eph["delta"])
    ang_width = np.average(np.deg2rad(eph["ang_width"]/sph))
    diam = ang_width*obs_dist*AU2M/1000  #diameter in km
    S = np.pi*(diam*mpkm/2)**2
    a = 10**((H-v_0-2.5*np.log10(np.pi/S)+mc)/-2.5)
    if a<1 and 1.5<sun_distance<4:
        I.append(a)
        sun_dist.append(sun_distance)
    if I[0] <= 0.1:
        T = 1 #C_type, carbonaceous
    if 0.1 < I[0] <= 0.18:
        T = 2  #"S_type, silicaceous or M_typem metallic"
    if I[0] > 0.18:
        T = 3   #"S_type, silicaceous"
    Type.append(T)
    n = n+1
    
fig, axs = plt.subplots(1, 1, figsize=(11, 6), sharex=False, sharey=False, gridspec_kw=None)
plt.tight_layout(pad=6)
    
axs.plot(sun_dist, I, "+", color="red")
axs.set(ylabel= "Albedo", xlabel= "Distance from the sun [AU]", title="Albedo of asteroids against their average distance from the sun")

z = np.polyfit(sun_dist, I, 1)
p = np.poly1d(z)
py.plot(sun_dist,p(sun_dist),"r")
plt.annotate("y=%.6fx+%.6f"%(z[0],z[1]), xy=(0.05, 0.95), xycoords='axes fraction')
#-0.045610x+0.268640

#y = -0.051142x+0.280705
#-0.059202x+0.302258
#0.063169