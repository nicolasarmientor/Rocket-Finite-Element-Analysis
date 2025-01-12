#######################################################
################### Library Imports ###################

import numpy as np
from mpmath import mp

#######################################################
################## Decimal Presicion ##################

mp.dps = 20 # Decimal  Precision

#######################################################
############## Gravitational Properties ###############

g_moon = mp.mpf(1.625)  # Gravity of Moon (m/s^2)
g_mars = mp.mpf(3.71)   # Gravity of Mars (m/s^2)

#######################################################
################ Rocket Specifications ################

m_wet = mp.mpf(557000.0) # Wet Mass (kg)
t_skin = mp.mpf(0.0020320) # Thickness of Skin (m)
r = mp.mpf(4) # Radius of Rocket (m)

#######################################################
################# Material Properties #################

### 304 Stainless Steel

E = mp.mpf(193e9)       # Young's Modulus (Pa)
G = mp.mpf(86e9)        # Shear Modulus (Pa)
v = mp.mpf(0.29)        # Poisson's Ratio
rho = mp.mpf(7900)      # Density (kg/m^3)
sigma_yield = mp.mpf(215e6)  # Yielding Stress (Pa)

#######################################################
################ Stiffener Properties #################

l_bottom = mp.mpf(0.0254) # Bottom Flange Length (m)
l_top = mp.mpf(0.0635) # Top Flange Length (m)
l_web = mp.mpf(0.0254) # Web Flange Length (m) 
l_total = mp.mpf(0.1651) # Total Length of Stiffenner (m)
t_stiff = mp.mpf(0.002032) # Thickness of Stiffener (m)

#######################################################
################ Variable Intialization ###############

n = int(10) # Number of Stiffeners
theta = [mp.mpf(0)] * n # Area Sector Angle
inclination = mp.mpf(35 * np.pi / 180)

#######################################################
################# Boundary Conditions #################

f_x = m_wet * g_mars * 3.0 + m_wet * g_mars * mp.cos(inclination) # Force X Applied (N)
f_y = -m_wet * g_mars * mp.sin(inclination) # Force Y Applied (N)
f_z = 0 # Force Z Applied (N)

x_f = mp.mpf(20.0) # Position of Cut (m)
y_f = 0 # Location of Force (m)
z_f = 0 # Location of Force (m)

v_y = f_y # Shear Force Y (N)
v_z = 0 # Shear Force Z (N)

m_z = f_y * x_f # Moment due to Shear Y (N-m)
m_y = f_z * x_f # Moment due to Shear Z (N-m)

#######################################################