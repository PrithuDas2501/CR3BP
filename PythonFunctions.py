import numpy as np
import math as m

## Taken from the Reference provided
def nondim_cr3bp(t, Y, pi_2, t0s, t0j, rs, rj, pi_s, pi_j): ## Every input is a non dimensional entity
    ang_pos_s = +t0s
    ang_pos_j = +t0j
    X_s = rs*m.cos(ang_pos_s)
    Y_s = rs*m.sin(ang_pos_s)
    X_j = rj*m.cos(ang_pos_j)
    Y_j = rj*m.sin(ang_pos_j)


    # Get the position and velocity from the solution vector
    x, y, z = Y[:3]
    xdot, ydot, zdot = Y[3:]

    # Define the derivative vector
    Ydot = np.zeros_like(Y)
    Ydot[:3] = Y[3:]

    sigma = np.sqrt(np.sum(np.square([x + pi_2, y, z])))
    psi = np.sqrt(np.sum(np.square([x - 1 + pi_2, y, z])))

    ## NEED TO ADD SUN AND JUPITER EFFECT TO THESE 3 LINES OF CODE
    ## Astropy Library has some provisions for getting positional data for celestial bodies.
    Ydot[3] = 2 * ydot + x - (1 - pi_2) * (x + pi_2) / sigma**3 - pi_2 * (x - 1 + pi_2) / psi**3 - pi_s*X_s/rs**3 - pi_j*X_j/rj**3
    Ydot[4] = -2 * xdot + y - (1 - pi_2) * y / sigma**3 - pi_2 * y / psi**3 - pi_s*Y_s/rs**3 - pi_j*Y_j/rj**3
    Ydot[5] = -(1 - pi_2)/sigma**3 * z - pi_2/psi**3 * z ## Assuming only planar motion, no Perturbation here
    return Ydot


'''
from astropy.coordinates import get_body
from astropy.coordinates import EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u


# Set your observer location (example: Greenwich Observatory)
observer_location = EarthLocation(lat=51.477*u.deg, lon=0*u.deg, height=0*u.m)

# Get the current UTC time
time_utc = Time.now()

# Get the Sun's position at the specified time and observer location
sun_position = get_body('sun',time_utc).transform_to(AltAz(obstime=time_utc, location=observer_location))
jupiter_position = get_body('jupiter', time_utc).transform_to(AltAz(obstime=time_utc, location=observer_location))

# Get the distance from the Earth to the Sun
distance_to_sun = sun_position.distance
distance_to_jupiter = jupiter_position.distance

# Print the distance (1 AU = 149,597,870.7 kilometers.)
print("Distance from Earth to Sun:", distance_to_sun)
print(sun_position)
'''