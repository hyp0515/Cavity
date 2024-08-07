import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from cavity_data import *

cavity = make_dict(target=with_cav+without_cav)


r_grid = np.logspace(np.log10(0.1), np.log10(200), 200, endpoint=True)

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
    sigma_vis = Mdot*(Msun/yr)*2.3*mp/(3*np.pi*alpha*kB*T_r)*np.sqrt(G*Mstar*Msun/(r_grid*au)**3)
    return sigma_vis

def temperature(Mstar, Mdot):
    n = 3*G*Mstar*Msun*Mdot*(Msun/yr)
    d = 8*np.pi*SB*(r_grid*au)**3
    return (n/d)**(1/4)

def vr_profile(Mstar, Mdot, Sigma):
    vr_cms  = Mdot*(Msun/yr)/(2*np.pi*Sigma*au*r_grid)
    kep_cms = np.sqrt(G*Mstar*Msun/(au*r_grid))
    return vr_cms/kep_cms

def alpha_profile(Mstar, Mdot, Sigma):
    T_disk = ((3*G*Mdot*(Msun/yr)*Mstar*Msun)/(8*np.pi*SB*(r_grid*au)**3))**(1/4)
    a = Mdot*(Msun/yr)*2.3*mp*np.sqrt(G*Mstar*Msun/((r_grid*au)**3))/(3*np.pi*Sigma*kB*T_disk)
    return a

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
    # sigma_vis = surface_density_viscous(Lstar=Lstar_list[i],
    #                                     Mstar=Mstar_list[i],
    #                                     Mdot=Mdot_list[i])
    # ax[i].plot(r_grid, sigma_vis, label = 'Manara+14')
    T = temperature(Mstar=Mstar_list[i],
                    Mdot=Mdot_list[i])
    
    v = vr_profile(Mstar=Mstar_list[i],
                   Mdot=Mdot_list[i],
                   Sigma=sigma)
    print(v.max(), v.min())  # HD139614 is weird
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, sigma, label = r'$\Sigma_{g}$')
    ax2.plot(r_grid, T, label = r'T', color = 'r')
    # ax1.scatter(r_grid, sigma)
    ax1.quiver(r_grid, sigma, -v, np.zeros(v.shape), v,
               cmap='rainbow')
    
    ax1.set_title(name_list[i], fontsize = 16)
    ax1.set_xlim((0.1, 200))
    ax1.set_ylim((1e-6, 1e6))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax1.set_xlabel('R [AU]', fontsize = 16)
    ax1.set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax2.set_ylabel('T [K]', fontsize = 16)
    ax1.legend(fontsize = 16, loc='upper left')
    ax2.legend(fontsize = 16, loc='upper right')
    
plt.tight_layout()
plt.savefig('profile.pdf', transparent = True)
plt.close()


fig, ax = plt.subplots(7, 4, figsize=(20, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    alpha_plot = alpha_profile(Mstar=Mstar_list[i],
                               Mdot=Mdot_list[i],
                               Sigma=sigma)
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, sigma, label = r'$\Sigma_{g}$')
    ax2.plot(r_grid, alpha_plot, label = r'$\alpha$', color = 'r')
    ax2.plot(r_grid, 1e-2*np.ones(r_grid.shape),'r--')
    
    ax1.set_title(name_list[i], fontsize = 16)
    ax1.set_xlim((0.1, 200))
    ax1.set_ylim((1e-6, 1e6))
    ax2.set_ylim((1e-4, 1e4))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
    
    ax1.set_xlabel('R [AU]', fontsize = 16)
    ax1.set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax2.set_ylabel(r'$\alpha$', fontsize = 16)
    ax1.legend(fontsize = 16, loc='upper left')
    ax2.legend(fontsize = 16, loc='upper right')
    
plt.tight_layout()
plt.savefig('alpha_profile.pdf', transparent = True)
plt.close()