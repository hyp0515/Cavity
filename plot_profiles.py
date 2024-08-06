import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from cavity_data import *

cavity = make_dict(target=with_cav+without_cav)


r_grid = np.logspace(np.log10(0.1), np.log10(200), 1000)


def surface_density(rc, sigmac, delta_gas, rcav, gamma = 1):
    if delta_gas is None and rcav is None:
        delta_gas, rcav = 1, 0
    density = np.where(
        r_grid < rcav,
        delta_gas * sigmac * ((r_grid/rc)**(-gamma)) * np.exp(-(r_grid/rc)**(2-gamma)),
        sigmac * ((r_grid/rc)**(-gamma)) * np.exp(-(r_grid/rc)**(2-gamma))
    )
    return density

def surface_density_viscous(Lstar, Mstar, Mdot, phi=0.02, alpha=1e-3):
    T_r = (phi*Lstar*Lsun/(8*np.pi*SB*(r_grid*au)**2))**(1/4)
    sigma_vis = Mdot*(Msun/yr)*2*mp/(3*np.pi*alpha*kB*T_r)*np.sqrt(G*Mstar*Msun/(r_grid*au)**3)
    return sigma_vis

# Parameters
name_list   = cavity['name']
rc_list     = cavity['rc']
sigmac_list = cavity['sigmac']
deltag_list = cavity['deltag']
rcavg_list  = cavity['rcavg']
Lstar_list  = cavity['Lstar']
Mstar_list  = cavity['Mstar']
Mdot_list   = cavity['Mdot']
gamma_list  = cavity['gamma']

fig, ax = plt.subplots(7, 4, figsize=(20, 30), sharex=True, sharey=True)
ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    sigma_vis = surface_density_viscous(Lstar=Lstar_list[i],
                                        Mstar=Mstar_list[i],
                                        Mdot=Mdot_list[i])
    ax[i].plot(r_grid, sigma, label = 'Andrews+11')
    ax[i].plot(r_grid, sigma_vis, label = 'Manara+14')
    ax[i].set_title(name_list[i], fontsize = 16)
    ax[i].set_xlim((0.1, 200))
    ax[i].set_ylim((1e-6, 1e6))
    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    if i == 24:
        ax[i].set_xlabel('R [AU]', fontsize = 16)
        ax[i].set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax[i].legend(fontsize = 16)
    

plt.tight_layout()
plt.savefig('profile.pdf', transparent = True)
plt.close()

    