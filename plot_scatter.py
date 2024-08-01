import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

from itertools import combinations
from cavity_data import cavity, key_to_text

df = pd.DataFrame(cavity)

params = list(df.columns)
params.remove('name')
params.remove('rc')
params.remove('sigmac')
# params.remove('age')
# params.remove('deltad')

param_combinations = list(combinations(params, 2))
num_combinations = len(param_combinations)

if num_combinations % 4 == 0:
    if num_combinations % 3 == 0:
        fig, axes = plt.subplots(6, num_combinations // 6 + num_combinations % 6, figsize=(40, 40))
    else: 
        fig, axes = plt.subplots(4, num_combinations // 4 + num_combinations % 4, figsize=(48, 30))
else:
    fig, axes = plt.subplots(3, num_combinations // 3 + num_combinations % 3, figsize=(45, 20))
axes = axes.flatten()

# Plot each combination
def linear_model(x, a, b):
    return a * x + b

for i, (param1, param2) in enumerate(param_combinations):
    ax = axes[i]
    # color_value = np.array(df['age'])/(np.array(df['Md'])/np.array(df['Mdot']))
    color_value = np.array(df['Mdot'])/np.array(df['Md'])
    ax.scatter(df[param1], df[param2], s=125, c=np.log10(color_value), cmap='rainbow')
    
    params, covariance = curve_fit(linear_model, np.log10(df[param1]), np.log10(df[param2]))
    a, b = params
    fitted_line = linear_model(np.log10(df[param1]), a, b)
    ax.plot(df[param1], 10**(fitted_line))
    
    ax.set_xlabel(key_to_text[param1], fontsize=25)
    ax.set_ylabel(key_to_text[param2], fontsize=25)
    ax.set_title(f'{key_to_text[param1]} vs {key_to_text[param2]}', fontsize=25)
    ax.set_xscale('log')
    ax.set_yscale('log')

# Remove any unused subplots
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

# Adjust layout
plt.tight_layout()

# Show the plot
plt.savefig('scatter.pdf', transparent = True)
plt.close()

