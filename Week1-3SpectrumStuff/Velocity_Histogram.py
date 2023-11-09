# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:19:45 2023

@author: Nickl
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

velocities = [29.967857949650206, 30.0360666791222, 29.6106663694946, 
              30.47205602161188, 30.327526656494776, 29.892115588352816, 
              29.14614508529477, 30.133613156710926, 30.159453300100875,
              30.202760899874747, 30.28307076697336, 30.276524612408984]
              #28.67315548311474, 27.37202032616293]

#k=0-k=6
#13CO (25)
#HC3N (27)
#CH3OH (29)
#H2CO (29)
#H2CO (31)
#k=7,8 but are # out as outliers

def velocity_graph(velocities):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    range_width = np.std(velocities)
    
    # Initialize an empty list to store data points
    data_points = []
    
    # Generate 100 data points for each velocity
    num_points = 100
    for velocity in velocities:
        data_points += list(np.linspace(velocity - range_width, velocity + range_width, num_points))
    
    # Create a histogram
    hist, bins, _ = plt.hist(data_points, bins=30, edgecolor='black', alpha=0.7, 
                             color = 'indigo')
    bin_centers = (bins[1:] + bins[:-1]) / 2
    
    # Perform the curve fit with better initial guesses
    mean_estimate = np.mean(data_points)
    stddev_estimate = np.std(data_points)
    params, _ = curve_fit(gaussian, bin_centers, hist, p0=[max(hist), mean_estimate, stddev_estimate])
    
    # Create the fitted Gaussian curve
    fit_curve = gaussian(bin_centers, *params)
    
    # Plot the histogram and the fitted Gaussian curve
    ax.plot(bin_centers, fit_curve, label='Gaussian Fit', color='gold')
    ax.plot([],[],label="v = {:.3f} km/s".format(params[1]), alpha =0)
    ax.plot([],[],label="$\sigma$ = {:.3f} km/s".format(params[2]), alpha=0)
    ax.set_title('Weighted Velocity Distribution')
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Frequency')
    ax.grid(True, alpha=0.2)
    ax.legend()
    ax.annotate(r'Mean = ' + str('{:.3f}'.format(np.mean(velocities))) + ' $\pm$ ' + str('{:.3f}'.format(range_width)) + ' km/s',
            xy=(0.5, 0.5), xytext=(0, -145), 
            xycoords='axes fraction', textcoords='offset points', fontsize=11,
            ha='center', va='top')
    ax.annotate(r'Weighted mean = ' + str('{:.3f}'.format(params[1])) + ' $\pm$ ' + str('{:.3f}'.format(params[2])) + ' km/s',
            xy=(0.5, 0.5), xytext=(0, -160), 
            xycoords='axes fraction', textcoords='offset points', fontsize=11,
            ha='center', va='top')
    
    # Display the Gaussian parameters (Amplitude, Mean, and Standard Deviation)
    print("\nWeighted final velocity = ",
          "{} +/- {} km/s".format(params[1], params[2]))
        
    plt.savefig('Weighted_velocity_distribution.png', dpi=1000, bbox_inches='tight')
    
    
    return

def gaussian(x, amplitude, mean, std_dev):
    
    return amplitude * np.exp(-(x - mean) ** 2 / (2 * std_dev ** 2))

velocity_graph(velocities)
