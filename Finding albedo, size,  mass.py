# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from astroquery.jplhorizons import Horizons

#DEFINITIONS BELOW!!!!!!!!!!
#https://ssd.jpl.nasa.gov/?horizons_doc#specific_quantities
#https://arxiv.org/pdf/1104.4227.pdf
"""
Constants
AU = Astronomical unit (m)
v_0 = Magnitude of the sun in the V band (mag)
spy = seconds per yes (s)
sph = seconds per hour (s)
mpkm = meters per kilometer (m)
"""

AU = 149597870700
v_0 = -26.762
spy = 3600*24*365
sph = 3600
mpkm = 1000
mc = -55.87

"""
Quoted values for Eckard 694 from https://ssd.jpl.nasa.gov/sbdb.cgi?sstr=694
e= eccentricity
sma = Semi-major axis (m)
i = Inclination (deg)
T = Period (s)
Trot = Rotational Period (s)
H = Absolute Magnitude (mag)
V = Apparent V band magnitude
dg = distance to the sun
dh distance to earth
D = Diameter (km)
Ag = Geometric albedo
G = Slope parameter (https://in-the-sky.org/data/object.php?id=A694#source_00)

Observed Values
dh = Heliocentric distance on night of observation
dg = geocentric distance on night of observation
V = Observed apparrent magnitude
"""

e = 0.3222011376356694
e_err = 3.1947*10**-8
sma = 2.671961685807055*AU
sma_err = AU*6.8881*10**-9
i = 15.88864056439503
i_err = 4.3093*10**-6
T = 4.37*spy
T_err = spy*1.689*10**-8
Trot = 5.925*sph
#H = 9.17
D = 121.891
D_err = 0.721
#Ag = 0.026                #Expected albedo
#Ag_err = 0.003
G = 0.15
dh = 1.076*AU
dg = 2.017*AU
d_err = 0.0005
V = 11.281199999999998
V_err = 0.16



"""
Plotting the phase of Ekard over 3500 days
"""
ekard = Horizons(694, epochs={'start':'2010-01-01', 'stop':'2020-01-01', 'step':'1d'})
eph = ekard.ephemerides()

fig, axs = plt.subplots(1, 1, figsize=(6, 4), sharex=False, sharey=False, gridspec_kw=None)
axs.plot(eph["datetime_jd"], eph["alpha"])
plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
axs.set_xlabel("JD")
axs.tick_params(rotation=45)
axs.set_ylabel(r"Phase angle ($^\circ$)")
axs.set(title = "Phase angle of Ekard over 3500 days")
plt.show()
minimum_alpha = min(eph["alpha"])
maximum_alpha = max(eph["alpha"])
#print("Min phase angle = {} deg" .format(minimum_alpha))
#print("Max phase angle = {} deg" .format(maximum_alpha))



"""
Estimating the geometic albedo
alpha = Phase angle of Eckard on night of observation using cosine formula
G1, = Gerometric albedo estimate (expect ~0.026)
Ha = Reduced magnitude at this phase angle
H0 = Reduced magnitude at opposition
S = averaage projected are of Eckard using quoted diameter
I = Radiance factor
"""
theta = np.linspace(0, maximum_alpha, num = 1000)         #Producing a list of angles to model the Radience factor
rads = np.deg2rad(theta)
alpha = np.arccos((dh**2+dg**2-AU**2)/(2*dh*dg))   #Calculating the Phase angle on the night of observation
Ha = V-5*np.log10(dh*dg/AU**2)                     #Calculating the reduced magnitude at this phase angle

#finding H0
P1 = np.exp(-3.332*(np.tan(alpha/2))**0.631)
P2 = np.exp(-1.862*(np.tan(alpha/2))**1.218)
H0 = Ha + 2.5*np.log10((1-G)*P1+G*P2)              #Absolute magnitude when alpha=0
#print("H0 = {}".format(H0))
#print("Ha = {}".format(Ha))

#finding how H varies with phase angle
P11 = np.exp(-3.332*(np.tan(rads/2))**0.631)
P22 = np.exp(-1.862*(np.tan(rads/2))**1.218)
Halpha = H0 - 2.5*np.log10((1-G)*P11+G*P22)        #How the absolute magnitude varies from opposition
S = np.pi*(D*mpkm/2)**2                            #Average projected area of the asteroid
I = 10**((Halpha-v_0-2.5*np.log10(np.pi/S)+mc)/-2.5)
I_err = I*np.sqrt((np.log(10)*V_err)**2+(2/(D/2)*np.log(10)*D_err/2)**2)

fig, axs = plt.subplots(1, 1, figsize=(6, 4), sharex=False, sharey=False, gridspec_kw=None)
ax = axs
ax2 = ax.twiny()
ax.plot(rads, I)
axs.tick_params(rotation=45)
ax.errorbar(rads,I, yerr = I_err, fmt = "none")
ax.set(ylabel='I/F', xlabel='Phase(Rads)', title='Radiance factor over the range of phase angles for Ekard')
ax2.plot(theta, I, ls='')
ax2.set(xlabel='Phase ($^\circ$)')
ax2.xaxis.set_major_locator(MultipleLocator(10))
ax2.xaxis.set_minor_locator(MultipleLocator(5))
plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
plt.show()
print("Geometric albedo estimate = {}±{}" .format(I[0], I_err[0]))



"""
Finding how THe intensity varies over time
"""
P111 = np.exp(-3.332*(np.tan(np.deg2rad(eph["alpha"])/2))**0.631)
P222 = np.exp(-1.862*(np.tan(np.deg2rad(eph["alpha"])/2))**1.218)
Halpha1 = H0 - 2.5*np.log10((1-G)*P111+G*P222)        #How the absolute magnitude varies from opposition
S = np.pi*(D*mpkm/2)**2                            #Average projected area of the asteroid
I1 = 10**((Halpha1-v_0-2.5*np.log10(np.pi/S)+mc)/-2.5)
#I_err = I*np.sqrt((np.log(10)*V_err)**2+(2/(D/2)*np.log(10)*D_err/2)**2)

fig, axs = plt.subplots(1, 1, figsize=(6, 4), sharex=False, sharey=False, gridspec_kw=None)
ax = axs
ax.plot(eph["datetime_jd"], Halpha1)
plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
axs.set_xlabel("JD")
axs.tick_params(rotation=45)
axs.set_ylabel("I/F")
axs.set(title = "Reduced magnitude of Eckard over 3500 days")
plt.show()


"""
Calculating the average diameter of the asteroid using the "1329" relation.
"""
D1 = 1329*I[0]**(-0.5)*10**(-0.2*Halpha1)
#print("Projected diameter on the night of observation= {}" .format(D1))
fig, axs = plt.subplots(1, 1, figsize=(6, 4), sharex=False, sharey=False, gridspec_kw=None)
ax = axs
ax.plot(eph["datetime_jd"], D1)
plt.tight_layout()
fig.align_ylabels(axs)
fig.align_xlabels(axs)
axs.set_xlabel("JD")
axs.tick_params(rotation=45)
axs.set_ylabel("Projected diameter (km)")
axs.set(title = "Projected diameter of Eckard over 3500 days")
plt.show()
"""
Calculating the three dimentions of the asteroid
assuming a triaxial ellipsoid shape, D is the average extent.
a is the longest extent
c is along the rotational axis
b is perpendicular to a and c
"""
a = round(D*(1.2**2*1.1)**(1/3), 3)
a_err = round(D_err*(1.2**2*1.1)**(1/3), 3)
b = round(a/1.2, 3)
b_err = round(a_err/1.2, 3)
c = round(b/1.1, 3)
c_err = round(a_err/(1.1*1.2), 3)

print("a = ({}±{}) km" .format(a, a_err))
print("b = ({}±{}) km" .format(b, b_err))
print("c = ({}±{}) km" .format(c, c_err))

"""
assuming a density of 2000kg/m^3, the mass can be found
"""
M = round(np.pi*a*b*c*2000/(6*10**9), 3)
M_err = round(M*np.sqrt((a_err/a)**2+(b_err/b)**2+(c_err/c)**2), 3)
print("Mass = ({}±{})E9 kg" .format(M, M_err))

