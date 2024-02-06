from scipy.integrate import solve_ivp
import numpy as np
import math as m
import matplotlib.pyplot as plt
import PythonFunctions as PF

from astropy.coordinates import get_body
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
# Set your observer location (example: Greenwich Observatory)
observer_location = EarthLocation(lat=51.477*u.deg, lon=0*u.deg, height=0*u.m)
# Get the current UTC time
time_utc = Time.now()
# Get positions at the specified time and observer location
sun_position = get_body('sun',time_utc).transform_to(AltAz(obstime=time_utc, location=observer_location))
jupiter_position = get_body('jupiter', time_utc).transform_to(AltAz(obstime=time_utc, location=observer_location))
moon_position = get_body('moon',time_utc).transform_to(AltAz(obstime=time_utc, location=observer_location))
# Get the distance 
distance_to_sun = sun_position.distance
distance_to_jupiter = jupiter_position.distance
distance_to_moon = moon_position.distance
# These masses represent the Earth-Moon-Sun-Jupiter system
m_1 = 5.974E24  # kg
m_2 = 7.348E22 # kg
m_j = 1.8982E27 # kg
m_s = 1.9891E30 # kg
pi_2 = m_2/(m_1 + m_2)
pi_s = m_s/(m_1 + m_2)
pi_j = m_j/(m_1 + m_2)

print(pi_2, pi_s, pi_j)
print(distance_to_sun/distance_to_moon,distance_to_jupiter/distance_to_moon)
X0 = np.array([0.87, 0, 0, 0, -1.48270, 0])
# X0 = np.array([1-pi_2,0.0455,0,-0.5,0.5,0.0])
x_0 = X0[0]
y_0 = X0[1]
z_0 = X0[2]
vx_0 = X0[3]
vy_0 = X0[4]
vz_0 = X0[5]

# Then stack everything together into the state vector
r_0 = np.array((x_0, y_0, z_0))
v_0 = np.array((vx_0, vy_0, vz_0))
Y_0 = np.hstack((r_0, v_0))

PF.nondim_cr3bp(0, Y_0, pi_2, m.pi/4, m.pi/4, distance_to_sun/distance_to_moon,distance_to_jupiter/distance_to_moon, pi_s, pi_j)