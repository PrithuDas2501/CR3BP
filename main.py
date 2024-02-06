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

#print(pi_2, pi_s, pi_j)
#print(distance_to_sun/distance_to_moon,distance_to_jupiter/distance_to_moon)
# Resonant EM Orbits
# X0 = [0.16686 0, 0, 0, 2.99389, 0]; # 1:1
# X0 = [0.17166 0, 0, 0, 3.00312, 0]; # 1:2
# X0 = [3.34922 0, 0, 0, -3.01527, 0]; # 1:3
# X0 = [0.33752 0, 0, 0, 1.86475, 0]; # 3:4
# X0 =  [0.44416, 0, 0, 0, 1.29134, 0]; # 3:2
# X0 = [0.10252, 0, 0, 0, 3.72303, 0]; # 4:1
# X0 = [1.18945, 0, 0, 0, -0.69908, 0]; # 4:3

X0 = np.array([0.87, 0.1, 0, 0, -1.48270, 0])
#X0 = np.array([0.87, 0.1, 0, 0, 0.2, 0])
#X0 = np.array([1-pi_2,0.0455,0,-0.5,0.5,0.0])
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

t_0 = 0  # nondimensional time
t_f = 20  # nondimensional time
t_points = np.linspace(t_0, t_f, 1000*(t_f-t_0))
sol = solve_ivp(PF.nondim_cr3bp, [t_0, t_f], Y_0, args=(pi_2,0,0,distance_to_sun/distance_to_moon,distance_to_jupiter/distance_to_moon,pi_s,pi_j), atol=1e-9, rtol=1e-6, t_eval=t_points)

Y = sol.y.T
r = Y[:, :3]  # nondimensional distance
v = Y[:, 3:]  # nondimensional velocity
t = sol.t

sol = solve_ivp(PF.nondim_cr3bp_without_pert, [t_0, t_f], Y_0, args=(pi_2,0,0,distance_to_sun/distance_to_moon,distance_to_jupiter/distance_to_moon,pi_s,pi_j), atol=1e-9, rtol=1e-6, t_eval=t_points)
Y = sol.y.T
r2 = Y[:, :3]  # nondimensional distance
v2 = Y[:, 3:]  # nondimensional velocity
################################################################3

x_2 = (1 - pi_2) * np.cos(np.linspace(0, np.pi, 100))
y_2 = (1 - pi_2) * np.sin(np.linspace(0, np.pi, 100))
fig, ax = plt.subplots(figsize=(5,5), dpi=96)

# Plot the orbits
plt.figure(1)
ax.plot(r[:, 0], r[:, 1], 'r', label="With Perturbation")
ax.plot(r2[:, 0], r2[:, 1], 'b', label="Without Perturbation")
plt.legend(["With Perturbation", "Without Perturbation"])
ax.axhline(0, color='k')
ax.plot(np.hstack((x_2, x_2[::-1])), np.hstack((y_2, -y_2[::-1])))
ax.plot(-pi_2, 0, 'bo', label="$m_1$")
ax.plot(1 - pi_2, 0, 'go', label="$m_2$")
ax.plot(x_0, y_0, 'ro')
ax.set_aspect("equal")
plt.show()  #ggg

plt.figure(2)
ax.plot(t, np.sum((r[:, :3] - r2[:, :3])**2, axis = 0), 'r', label="Deviation")


