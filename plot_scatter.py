import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

from itertools import combinations
from cavity_data import cavity

# Creating a DataFrame
df = pd.DataFrame(cavity)

# Create a list of parameter names
params = list(df.columns)
params.remove(r'$r_{c}$')
params.remove(r'$\Sigma_{c}$')
params.remove(r'$Myr$')
# params.remove(r'$\delta_{gas}$')
# Generate all combinations of pairs of parameters
param_combinations = list(combinations(params, 2))

# Calculate the number of combinations
num_combinations = len(param_combinations)

# Create subplots
fig, axes = plt.subplots(3, num_combinations // 3 + num_combinations % 3, figsize=(45, 20))

# Flatten the axes array for easy iteration
axes = axes.flatten()

# Plot each combination
def linear_model(x, a, b):
    return a * x + b

for i, (param1, param2) in enumerate(param_combinations):
    ax = axes[i]
    ratio = np.array(df[r'$Myr$'])/(np.array(df[r'$M_{d}$'])/np.array(df[r'$\dot{M}$'])/333000)
    ax.scatter(df[param1], df[param2], s=125, c=np.log10(ratio), cmap='rainbow')
    
    params, covariance = curve_fit(linear_model, np.log10(df[param1]), np.log10(df[param2]))
    a, b = params
    fitted_line = linear_model(np.log10(df[param1]), a, b)
    ax.plot(df[param1], 10**(fitted_line))
    
    ax.set_xlabel(param1, fontsize=25)
    ax.set_ylabel(param2, fontsize=25)
    ax.set_title(f'{param1} vs {param2}', fontsize=25)
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

