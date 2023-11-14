# -*- coding: utf-8 -*-
"""
Created on Thu Nov 9 15:24:03 2023

@author: nickl

this code takes the second and third columns from rings_final2 (the Radius(arcsec)
and VROT(KM/S) columns) and plots them against a v prop r^(-0.5) rotation curve to
estimate the mass of a central body inhibiting keplarian motion on nearby gas/dust
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# User-defined global scale factor
pc = 3.086e+16
G = 6.6743e-11
SunMass = 1.98847e30
scale_factor_x =  8*10**3 * np.pi * pc / (3600 * 180) ## turns arcsec into m
scale_factor_y = 1000 #turns km/s into m/s

FOLDER = 'e_sma1'
FILEPATH = 'C:/Users/nickl/linux/prog/bbarolo/output/' #NICK
# FILEPATH = 'C:/Users/christdick/linux/prog/bbarolo/output/' #NICK
# Define the model function
def model(x, A):
    return A * x**(-0.5)

# Load the data from rings_final2.txt
data = np.genfromtxt(FILEPATH+FOLDER +'/rings_final2.txt', usecols=(1, 2), skip_header=1, unpack=True)

# Unpack the radius and rotation velocity
rad, vrot = data

# Filter out zero (or close to zero) and negative radius values as they will cause issues in the model
valid_indices = (rad > 0) #| (rad > 0.028)
rad = rad[valid_indices]
print(rad)
vrot = vrot[valid_indices]

# Convert vrot to m/s for the second plot
vrot_m_s = vrot * scale_factor_y

# Provide an initial guess for A, for example, the mean of vrot
initial_guess = [np.mean(vrot)]

# Fit the model to the data
params, covariance = curve_fit(model, rad, vrot, p0=initial_guess, maxfev=10000)

# Unpack the fitting parameters
A = params[0]
#print(A)

# Adjust A for the scaled x-axis
A_scaled = A / scale_factor_x**(-0.5) * scale_factor_y

print(f'The value for A in the scaled plot, equal to sqrt(GM), is {A_scaled:.2e}')
Mass = A_scaled**2 / G
SMass = A_scaled**2 / (G * SunMass)
print(f'\nThis gives a value M, equal to A^2/G as {Mass:.3e} or {SMass:.3e} Solar Masses')
# Create a finer array for plotting the smooth curve
rad_fine = np.linspace(np.min(rad), np.max(rad), 1000)

# Create a figure with two subplots, sharing the y-axis
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 5))

# First subplot: original data
ax1.plot(rad, vrot, 'o', color='blue', label='Rotation Velocity')
ax1.plot(rad_fine, model(rad_fine, A), '-', color='red', label=f'Fit: ${A:.3f} \cdot r^{{ -0.5}}$')
ax1.set_xlabel('Radius (arcsec)')
ax1.set_ylabel('Rotation Velocity (km/s)')
ax1.set_title('3DBbarolo Rotation Curve')
ax1.grid(alpha=0.3)
ax1.legend()

# Second subplot: scaled radius by the user-defined scale factor, velocity in m/s
scaled_rad = scale_factor_x * rad  # Multiply the radius by the scale factor
ax2.plot(scaled_rad, vrot_m_s, 'o', color='green', label='Rotation Velocity (Scaled)')
ax2.plot(scale_factor_x * rad_fine, model(rad_fine, A) * scale_factor_y, '-', color='orange',
         label=f'$Fit: {A_scaled:.2e} \cdot r^{{ -0.5}}$',)
ax2.plot([],[], label = f'M = {SMass:.2f} M$_\odot$')
ax2.set_xlabel('Radius (m)')
ax2.set_ylabel('Rotation Velocity (m/s)')
ax2.set_title('SI Scaled Rotation Curve')
ax2.grid(alpha=0.3)
ax2.legend()

plt.savefig(fname = 'bbarolo_curve_fit.png', bbox_inches = 'tight', dpi = 866)
# Display the plot
plt.show()

