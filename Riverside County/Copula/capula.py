# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 01:30:53 2024

@author: shosseini
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma, norm, multivariate_normal, chi2
from scipy.optimize import minimize

# Given vectors
BR = [0.3391839172095941, 0.338794435197385, 0.34007302165611863, 0.3362855268773844, 0.2789834765504536, 
      0.3379585059683711, 0.33978216281466317, 0.33660871057147745, 0.33996129374371004, 0.33429016649186233, 
      0.34023816320000005, 0.33754852562390825, 0.3331133708766275, 0.1285472208145385]
N = [884592676.96119, 116317161.89268737, 434580939.82169616, 178390924.86499807, 248369696.0773816, 
     251197317.63641196, 451884624.0382529, 175645339.1375478, 681558678.3087615, 359591223.9482173, 
     257343486.11129397, 1142627080.9261448, 251226280.73109725, 207758585.71434098]


BR_array = np.array(BR)
N_array = np.array(N)
def fit_gamma(data):
    def negative_log_likelihood(params):
        shape, scale = params
        return -np.sum(gamma.logpdf(data, shape, scale=scale))
    
    initial_params = [1.0, np.mean(data)]
    result = minimize(negative_log_likelihood, initial_params, bounds=[(0, None), (0, None)])
    return result.x
params_N = fit_gamma(N_array)
shape_N, scale_N = params_N

params_BR = fit_gamma(BR_array)
shape_BR, scale_BR = params_BR

F_N = gamma.cdf(N_array, shape_N, scale=scale_N)
F_BR = gamma.cdf(BR_array, shape_BR, scale=scale_BR)

U = norm.ppf(F_N)
V = norm.ppf(F_BR)

correlation_coefficient = np.corrcoef(BR_array, N_array)[0, 1]
mean = [0, 0]
cov = [[1, correlation_coefficient], [correlation_coefficient, 1]]
rv = multivariate_normal(mean, cov)

n = np.linspace(0, 2*max(N_array), 100)
br = np.linspace(0, 2*max(BR_array), 100)
N_grid, BR_grid = np.meshgrid(n, br)

U_grid = norm.ppf(gamma.cdf(N_grid, shape_N, scale=scale_N))
V_grid = norm.ppf(gamma.cdf(BR_grid, shape_BR, scale=scale_BR))
pos = np.dstack((U_grid, V_grid))
Z = rv.pdf(pos)
probability_levels = [0.2, 0.5, 0.8]
chi2_values = [chi2.ppf(level, df=2) for level in probability_levels]
contour_levels = [rv.pdf([0, 0]) * np.exp(-0.5 * chi2_val) for chi2_val in chi2_values]
fig, ax = plt.subplots(figsize=(12, 8))
contour_filled = ax.contourf(N_grid, BR_grid, Z, levels=50, cmap='viridis')
plt.colorbar(contour_filled, ax=ax, label='Joint PDF')
contour_lines = plt.contour(N_grid, BR_grid, Z, levels=sorted(contour_levels), colors=['red', 'orange', 'black'])
ax.clabel(contour_lines,  fmt={contour_levels[0]: '0.2', contour_levels[1]: '0.5', contour_levels[2]: '0.8'}, inline=True, fontsize=10)
ax.set_xlabel('N')
ax.set_ylabel('BR')
ax.set_title('Joint PDF with Contour Levels')
plt.grid(True)



x0 = 600000000  # Example value for N
y0 = 0.24 
plt.plot(x0, y0, 'r*')
plt.show()
i=0

A = contour_lines.collections[0].get_paths()
for path in A:  # Loop through each path
    if path.contains_point((x0, y0)):  # Check if the point is within the path
        i=i+1
        break  # If the point is found in any path, you may want to break the loop
A = contour_lines.collections[1].get_paths()
for path in A:  # Loop through each path
    if path.contains_point((x0, y0)):  # Check if the point is within the path
        i=i+1
        break  # If the point is found in any path, you may want to break the loop
A = contour_lines.collections[2].get_paths()
for path in A:  # Loop through each path
    if path.contains_point((x0, y0)):  # Check if the point is within the path
        i=i+1
        break  # If the point is found in any path, you may want to break the loop
print(i)
if i==1:print("yellow") 
elif i==2:print("orange")
elif i==3: print("red")
elif i==0:print("green") 

