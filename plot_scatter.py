import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from scipy.optimize import curve_fit
import corner

from itertools import combinations
from cavity_data import *
###################################################################################################
def plot_all_scatter(cavity_dict, fname='scatter.pdf'):
    df = pd.DataFrame(cavity_dict)

    params = list(df.columns)
    params.remove('name')
    params.remove('rc')
    params.remove('sigmac')
    params.remove('age')
    params.remove('existence')
    params.remove('companion')
    params.remove('rcavd')
    params.remove('deltad')
    params.remove('gamma')
    params.remove('psi')
    # params.remove('vr')
    params.remove('hc')
    
    param_combinations = list(combinations(params, 2))
    num_combinations = len(param_combinations)

    if num_combinations % 4 == 0:
        if num_combinations % 3 == 0:
            fig, axes = plt.subplots(6, num_combinations // 6 + num_combinations % 6, figsize=(40, 40))
        else: 
            fig, axes = plt.subplots(4, num_combinations // 4 + num_combinations % 4, figsize=(60, 30))
    else:
        fig, axes = plt.subplots(3, num_combinations // 3 + num_combinations % 3, figsize=(45, 20))
    axes = axes.flatten()

    # Plot each combination
    def linear_model(x, a, b):
        return a * x + b

    for i, (param1, param2) in enumerate(param_combinations):
        ax = axes[i]
        
        color_value = df['age']/(df['age']+(df['Md']/(df['Mdot']*333000)/1e6))
        
        marker_dict = {
            (True, False)  : 'o',
            (True, True)   : '^',
            (False, False) : 'x'
        }
        for j in range(len(df['name'])):
            status = (df['existence'][j], df['companion'][j])
            scatter = ax.scatter(df[param1][j], df[param2][j], s=150, c=color_value[j],
                                cmap='rainbow', marker = marker_dict[status],
                                vmin=color_value.min(), vmax=color_value.max())
        
        # parameters, covariance = curve_fit(linear_model, np.log10(df[param1]), np.log10(df[param2]))
        # a, b = parameters
        # fitted_line = linear_model(np.log10(df[param1]), a, b)
        # ax.plot(df[param1], 10**(fitted_line))
        try:
            ax.set_xlabel(key_to_text[param1]+key_to_unit[param1], fontsize=25)
        except:
            ax.set_xlabel(key_to_text[param1], fontsize=25)
        try:
            ax.set_ylabel(key_to_text[param2]+key_to_unit[param2], fontsize=25)
        except:
            ax.set_ylabel(key_to_text[param2], fontsize=25)
        
        ax.set_title(f'{key_to_text[param1]} vs {key_to_text[param2]}', fontsize=25)
        ax.set_xscale('log')
        ax.set_yscale('log')

    # Remove any unused subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    cbar_ax = fig.add_axes([0.95, 0.05, 0.01, 0.9])
    cbar = fig.colorbar(scatter, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Proportion of life', fontsize=30)
    cbar.ax.tick_params(labelsize=20)
    
    legend_elements = [
        mlines.Line2D([], [], color='black', marker='o', linestyle='None', markersize=10, label='w/ cav + w/o companion'),
        mlines.Line2D([], [], color='black', marker='^', linestyle='None', markersize=10, label='w/ cav + w/ companion'),
        mlines.Line2D([], [], color='black', marker='x', linestyle='None', markersize=10, label='Full disks')
    ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=20)
    fig.tight_layout(rect=[0, 0, 0.95, 1])
    
    plt.savefig(fname, transparent = True)
    plt.close()
    return

def plot_corner(cavity_dict, fname='corner.pdf'):
    df = pd.DataFrame(cavity_dict)

    params = list(df.columns)
    params.remove('name')
    params.remove('rc')
    params.remove('sigmac')
    # params.remove('age')
    # params.remove('deltad')  

    df_selected = df[params].apply(np.log10)
    figure = corner.corner(df_selected, labels=[key_to_text[param] for param in params],
                        plot_contours=False,
                        show_titles=True, title_kwargs={"fontsize": 12})
    figure.savefig(fname, transparent=True)
    plt.close()

def plot_scatter(data_set=None, x_axis=None , y_axis=None,
                              x_scale=None, y_scale=None,
                              x_label=None, y_label=None,
                              title=None  , fname=None):
    if title is None: title = ''
    if fname is None: fname = 'test.pdf'
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    

    try:
        x = data_set[x_axis]
        y = data_set[y_axis]
        print('Using the given data dictionary')
    except:
        x = x_axis
        y = y_axis
        print('Using the given data arrays')
        if len(x) != len(y):
            print('The two data arrays have different lengths')
            print(f'The length of x-axis is {len(x)}')
            print(f'The length of y-axis is {len(y)}')
            return
        
    ax.scatter(x, y)
    
    try:
        xlabel = key_to_text[x_axis] + key_to_unit[x_axis]
        ylabel = key_to_text[y_axis] + key_to_unit[y_axis]
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    except:
        try:
            xlabel = key_to_text[x_axis]
            ylabel = key_to_text[y_axis]
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        except:
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
        pass
    
    try:
        ax.set_xscale(x_scale)
    except:
        ax.set_xscale('log')
        print('x-axis is automatically set to be log scaling')
    
    try:
        ax.set_yscale(y_scale)
    except:
        ax.set_yscale('log')
        print('y-axis is automatically set to be log scaling')
        
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(fname, transparent = True)
    return
    

###################################################################################################
# cavity = make_dict(target=with_companion+without_companion)
cavity = make_dict(target=with_cav+without_cav)
plot_all_scatter(cavity_dict=cavity, fname='scatter_raw.pdf')

cavity = make_dict(target=with_companion)
plot_all_scatter(cavity_dict=cavity, fname='scatter_companion.pdf')

cavity = make_dict(target=without_companion)
plot_all_scatter(cavity_dict=cavity, fname='scatter_nocompanion.pdf')



# cavity = make_dict(target=with_cav+without_cav)
# plot_scatter(data_set=cavity, x_axis='Mdot', y_axis='Lstar', title='Mdot vs Lstar')

