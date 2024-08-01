import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cavity_data import target, cavity

au = 14959787069100

r_grid = np.logspace(np.log10(0.1), np.log10(200), 1000)


def surface_density(rc, sigmac, delta_gas, rcav, gamma = 1):
    density = np.where(
        r_grid < rcav,
        delta_gas * sigmac * ((r_grid/rc)**(-gamma)) * np.exp(-(r_grid/rc)**(2-gamma)),
        sigmac * ((r_grid/rc)**(-gamma)) * np.exp(-(r_grid/rc)**(2-gamma))
    )
    return density


# Parameters
name_list = list(target.keys())
rc_list = cavity[r'$r_{c}$']
sigmac_list = cavity[r'$\Sigma_{c}$']
deltagas_list = cavity[r'$\delta_{gas}$']
rcav_list = cavity[r'$r_{cav}$']



fig, ax = plt.subplots(5, 4, figsize=(20, 20), sharex=True, sharey=True)
ax = ax.flatten()

for i in range(len(rc_list)):
    if i == 5:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltagas_list[i],
                                rcav=rcav_list[i],
                                gamma=0.3)
    elif i == 16:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltagas_list[i],
                                rcav=rcav_list[i],
                                gamma=1.1)
    else:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltagas_list[i],
                                rcav=rcav_list[i],
                                gamma=1)
    ax[i].plot(r_grid, sigma)
    ax[i].set_title(name_list[i], fontsize = 16)
    ax[i].set_xlim((0.1, 200))
    ax[i].set_ylim((1e-6, 1e6))
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    if i == 16:
        ax[i].set_xlabel('R [AU]', fontsize = 16)
        ax[i].set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    

plt.tight_layout()
plt.savefig('profile.pdf', transparent = True)
plt.close()

    