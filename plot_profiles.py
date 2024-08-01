import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from cavity_data import cavity

au   = 1.496e13
SB   = 5.6704e-5
mp   = 1.6726e-24
kB   = 1.3807e-16
Msun = 1.9884e33
Lsun = 3.828e33
G    = 6.6743e-8
yr   = 3.156e7

r_grid = np.logspace(np.log10(0.1), np.log10(200), 1000)


def surface_density(rc, sigmac, delta_gas, rcav, gamma = 1):
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



fig, ax = plt.subplots(5, 4, figsize=(20, 20), sharex=True, sharey=True)
ax = ax.flatten()

for i in range(len(rc_list)):
    if i == 5:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltag_list[i],
                                rcav=rcavg_list[i],
                                gamma=0.3)
    elif i == 16:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltag_list[i],
                                rcav=rcavg_list[i],
                                gamma=1.1)
    else:
        sigma = surface_density(rc = rc_list[i],
                                sigmac=sigmac_list[i],
                                delta_gas=deltag_list[i],
                                rcav=rcavg_list[i],
                                gamma=1)
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
    if i == 16:
        ax[i].set_xlabel('R [AU]', fontsize = 16)
        ax[i].set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax[i].legend(fontsize = 16)
    

plt.tight_layout()
plt.savefig('profile.pdf', transparent = True)
plt.close()

    