
#######################################################
################### Library Imports ###################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpmath import mp
from rocket_parameters import n, r, t_skin, t_stiff, l_total, f_x, m_z, m_y

#######################################################
################## Decimal Presicion ##################

mp.dps = 20 # Decimal  Precision

#######################################################
################ Calculation Functions ################

### Angles Array Function ###

def calculate_angles(n):
    theta = [mp.mpf(0)] * n
    delta_theta = mp.mpf(2 * np.pi / n)
    for i in range(1, n):
        theta[i] = theta[i-1] + delta_theta
    return theta

### Z - Coordinate Array Function ###

def calculate_z_coordinates(n, theta, r):
    z = [mp.mpf(0)] * n

    for i in range(len(theta)):
        z[i] = r * mp.sin(theta[i])
    
    return z

### Y - Coordinate Array Function ###

def calculate_y_coordinates(n, theta, r):
    y = [mp.mpf(0)] * n

    for i in range(len(theta)):
        y[i] = r * mp.cos(theta[i])
    
    return y

### Lumped Area Array Function ###

def calculate_lumped_areas(n, t_stiff, t_skin, r, l_total):
    delta_theta = mp.mpf(2 * np.pi / n)
    area = [r * delta_theta * t_skin + l_total * t_stiff] * n
    
    return area

### Centroids Function ###

def calculate_centroids(z, y, area):
    total_area = mp.fsum(area)

    z_g = mp.fsum([a * z_i for a, z_i in zip(area, z)]) / total_area
    y_g = mp.fsum([a * y_i for a, y_i in zip(area, y)]) / total_area

    return z_g, y_g

### Moments of Inertia Function ###

def calculate_moments_of_inertia(z, y, area, z_g, y_g):
    i_gy = mp.fsum([a * (z_i - z_g)**2 for a, z_i in zip(area, z)])
    i_gz = mp.fsum([a * (y_i - y_g)**2 for a, y_i in zip(area, y)])
    i_gyz = mp.fsum([a * (y_i - y_g) * (z_i - z_g) for a, y_i, z_i in zip(area, y, z)])

    return i_gy, i_gz, i_gyz

### Flexure Formula Array Function ###

def apply_flexure_formula(i_gy, i_gz, i_gyz, m_z, m_y, y_g, z_g, area, y, z, f_x):
    sigma_flexure = [(f_x / mp.fsum(area) + 1/(i_gy * i_gz - i_gyz**2) * 
               (-(i_gy * m_z + i_gyz * m_y) * (y_i - y_g) + 
                (i_gz * m_y + i_gyz * m_z) * (z_i - z_g)))
                for y_i, z_i in zip(y, z)]

    return sigma_flexure

### Method Calling ###

theta = calculate_angles(n)
z = calculate_z_coordinates(n, theta, r)
y = calculate_y_coordinates(n, theta, r)
area = calculate_lumped_areas(n, t_stiff, t_skin, r, l_total)
z_g, y_g = calculate_centroids(z, y, area)
i_gy, i_gz, i_gyz = calculate_moments_of_inertia(z, y, area, z_g, y_g)
sigma_flexure = apply_flexure_formula(i_gy, i_gz, i_gyz, m_z, m_y, y_g, z_g, area, y, z, f_x)

### Results Print Statement ###

data = {
    'sigma_flexure (Pa)': [float(sf) for sf in sigma_flexure]
}

coordinates = {
    'y_coordinates (m)': y,
    'z_coordinates (m)': z
}

results = pd.DataFrame(data)
coord = pd.DataFrame(coordinates)

print()
print()
print(results)
print()
print()
print(coord)
print()
print()