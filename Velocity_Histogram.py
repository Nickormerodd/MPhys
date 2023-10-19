# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:28:01 2023

@author: nickl

The code takes the velocities that the user puts in the first
bit and then takes a selected range (the standard deviation)
around those velocity points and plots 100 points in the range.
Then it adds them to a histogram and plots a curve fit gaussian
around it.

It's quite nice cause it's automatically a weighted mean, so
if you point in anomolous points then it will almost ignore
them.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define your 6 velocities
velocities = [29.96785795, 30.03606668, 29.61066637, 30.47205602,
              30.32752666, 29.89211559, 29.14614509, 28.67315548,
              27.37202032]

# Define the range (+/- 0.25 around each velocity)
range_width = 0.411

# Initialize an empty list to store data points
data_points = []

# Generate 100 data points for each velocity
num_points = 100
for velocity in velocities:
    data_points += list(np.linspace(velocity - range_width, velocity + range_width, num_points))

# Create a histogram
hist, bins, _ = plt.hist(data_points, bins=20, edgecolor='black', alpha=0.7)
bin_centers = (bins[1:] + bins[:-1]) / 2

# Define the Gaussian function
def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

# Perform the curve fit with better initial guesses
mean_estimate = np.mean(data_points)
stddev_estimate = np.std(data_points)
params, _ = curve_fit(gaussian, bin_centers, hist, p0=[max(hist), mean_estimate, stddev_estimate])

# Create the fitted Gaussian curve
fit_curve = gaussian(bin_centers, *params)

# Plot the histogram and the fitted Gaussian curve
plt.plot(bin_centers, fit_curve, 'r-', label='Gaussian Fit')
plt.plot([],[], label = 'v = {:.3f}'.format(params[1]), alpha = 0)
plt.plot([],[], label = '$\sigma$ = {:.3f}'.format(params[2]), alpha = 0)
plt.title('Velocity Distribution with Gaussian Fit')
plt.xlabel('Velocity')
plt.ylabel('Frequency')
plt.grid(True, alpha = 0.5)
plt.legend()

plt.savefig('Velocity_Histogram.png', dpi = 1000)

# Display the Gaussian parameters (Amplitude, Mean, and Standard Deviation)
print("Amplitude (A):", params[0])
print("Mean (mu):", params[1])
print("Standard Deviation (sigma):", params[2])
print("\nThe final velocity (with a weighted mean) is",
      " {:.3f} +/- {:.3f} km/s".format(params[1], params[2]))

# Show the plot
plt.show()
