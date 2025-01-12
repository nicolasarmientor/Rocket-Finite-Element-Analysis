#######################################################
################### Library Imports ###################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpmath import mp
from lumped_area import sigma_flexure, y, z
from rocket_parameters import l_bottom, l_top, l_web, t_stiff, sigma_yield, E, v

#######################################################
################## Decimal Presicion ##################

mp.dps = 20 # Decimal Places

#######################################################
################ Calculation Functions ################

### Warning: Hat Stiffeners Only! ###

### Stiffener Areas Function ###

def calculate_stiffener_areas(l_bottom, l_top, l_web, t_stiff):
    stiff_areas = [mp.mpf(0) for _ in range(5)]
    
    stiff_areas[0] = (l_bottom - t_stiff) * t_stiff
    stiff_areas[1] = l_web * t_stiff
    stiff_areas[2] = (l_top - 2 * t_stiff) * t_stiff
    stiff_areas[3] = stiff_areas[1]
    stiff_areas[4] = stiff_areas[0]

    return stiff_areas

### Stiffener Area Centroid Function ###

def calculate_stiffener_area_centroids(l_web, t_stiff):
    area_centroids = [mp.mpf(0) for _ in range(5)]

    area_centroids[0] = t_stiff * 0.5
    area_centroids[1] = l_web * 0.5
    area_centroids[2] = (l_web - t_stiff * 0.5)
    area_centroids[3] = area_centroids[1]
    area_centroids[4] = area_centroids[0]

    return area_centroids

### Stiffener Centroid Function ###

def calculate_stiffener_centroid(stiff_areas, area_centroids):
    y_g_stiff = mp.fsum([a * y for a, y in zip(stiff_areas, area_centroids)]) / mp.fsum(stiff_areas)

    return y_g_stiff

### Crippling Function ###

def calculate_crippling(t_stiff, l_bottom, l_top, l_web, sigma_yield, E, v):
    b_crippling = [mp.mpf(0) for _ in range(5)]
    sigma_co = sigma_yield
    C = [mp.mpf(0.425), mp.mpf(0.425), mp.mpf(4.0), mp.mpf(0.425), mp.mpf(0.425)]

    b_crippling[0] = l_bottom - t_stiff
    b_crippling[1] = b_crippling[0]
    b_crippling[2] = l_top - 2 * t_stiff
    b_crippling[3] = l_web - 2 * t_stiff
    b_crippling[4] = b_crippling[3]

    lam = [mp.sqrt(sigma_yield / E) * (b / t_stiff) for b in b_crippling]

    n = mp.mpf(0.6)
    alpha = mp.mpf(0.8)
    
    sigma_crippling = [sigma_yield * alpha * ((np.pi**2 * c) / (12 * (1 - v**2) * l))**(1 - n) for c, l in zip(C, lam)]

    for i in range(len(sigma_crippling)):
        if sigma_crippling[i] > sigma_co:
            sigma_crippling[i] = sigma_co

    sigma_cc = mp.fsum([b * t_stiff * s for b, s in zip(b_crippling, sigma_crippling)]) / mp.fsum([b * t_stiff for b in b_crippling])

    ms_cc = np.abs([(sigma_cc - flexure) / flexure for flexure in sigma_flexure])

    return sigma_crippling, sigma_cc, ms_cc

### Buckling Function ###

def calculate_buckling(sigma_cc, v, E, l_web, t_stiff, sigma_flexure):
    sigma_cr = sigma_cc

    c = mp.mpf(1.5)
    K = mp.mpf(0.85)
    L = mp.mpf(14.5)


    b_buckling = [mp.mpf(0) for _ in range(5)]
    h_buckling  = [mp.mpf(0) for _ in range (5)]
    y_c = [mp.mpf(0) for _ in range (5)]
    C = [mp.mpf(0.425), mp.mpf(0.425), mp.mpf(4.0), mp.mpf(0.425), mp.mpf(0.425)]

    b_buckling[0] = l_bottom + t_stiff
    b_buckling[1] = l_bottom + t_stiff
    b_buckling[3] = l_bottom + t_stiff
    b_buckling[4] = l_bottom + t_stiff
    b_buckling[2] = t_stiff

    h_buckling[0] = t_stiff
    h_buckling[1] = t_stiff
    h_buckling[3] = t_stiff
    h_buckling[4] = t_stiff
    h_buckling[2] = l_top

    y_c[0] = t_stiff/2
    y_c[1] = t_stiff/2
    y_c[2] = l_web - t_stiff/2
    y_c[3] = l_web
    y_c[4] = l_web

    for i in range(1000):
        area = mp.fsum([b * h for b, h in zip(b_buckling, h_buckling)])
        y_g = mp.fsum([b * h * yc for b, h, yc in zip(b_buckling, h_buckling, y_c)])/area
        i_z = mp.fsum([(1/12) * b * h**3 + b * h * (yc - y_g)**2 for b, h, yc in zip(b_buckling, h_buckling, y_c)])
        rho_buckling = mp.sqrt(i_z/area)
        l_e = L/mp.sqrt(c)
        slend = l_e/rho_buckling
        slend_cr = np.pi * mp.sqrt(2/(sigma_cr/E)) 

        if slend <= slend_cr:
            sigma_cr = sigma_cc * (1 - 1/(4 * np.pi**2) * sigma_cr/E*slend**2) # Johnson's Formula
        else:
            sigma_cr = np.pi**2 * E/slend**2
        
        w_e = t_stiff * K * E / mp.sqrt(E) * mp.sqrt(1/sigma_cr)

        # Iterate back

    ms_cr = np.abs([(sigma_cr - flexure) / flexure for flexure in sigma_flexure])

    return sigma_cr, ms_cr



stiff_areas = calculate_stiffener_areas(l_bottom, l_top, l_web, t_stiff)
area_centroids = calculate_stiffener_area_centroids(l_web, t_stiff)
y_g_stiff = calculate_stiffener_centroid(stiff_areas, area_centroids)
sigma_crippling, sigma_cc, ms_cc = calculate_crippling(t_stiff, l_bottom, l_top, l_web, sigma_yield, E, v)
sigma_cr, ms_cr = calculate_buckling(sigma_cc, v, E, l_web, t_stiff, sigma_flexure)

data = {
    'ms_cc': [float(ms) for ms in ms_cc],
    'ms_cr': [float(ms) for ms in ms_cr]
}

results = pd.DataFrame(data)

print(results)
print()
print()

### Stiffener Position Plot ###

plt.figure(figsize=(8, 8))
plt.scatter([float(yi) for yi in y], [float(zi) for zi in z], c='blue', label='Stiffeners')

plt.axhline(0, color='black', linewidth=0.5, linestyle='--')  # Horizontal axis
plt.axvline(0, color='black', linewidth=0.5, linestyle='--')  # Vertical axis
plt.xlabel('y (Coordinate)')
plt.ylabel('z (Coordinate)')
plt.title('Stiffener Placement')
plt.grid(True)

plt.legend(loc='center')

plt.gca().set_aspect('equal', adjustable='box')
plt.show()