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

def temperature(psi, Lstar):
    T_r = (psi*Lstar*Lsun/(8*np.pi*SB*(r_grid*au)**2))**(1/4)
    return T_r

def vr_profile(Mstar, Mdot, Sigma, T):
    vr_cms  = Mdot*(Msun/yr)/(2*np.pi*Sigma*au*r_grid)
    kep_cms = np.sqrt(G*Mstar*Msun/(au*r_grid))
    cs_cms  = (kB*T/(2.3*mp))**0.5
    return vr_cms/kep_cms, vr_cms/cs_cms

def alpha_profile(Mstar, Mdot, Sigma, psi, Lstar):
    T_disk = (psi*Lstar*Lsun/(8*np.pi*SB*(r_grid*au)**2))**(1/4)
    a = Mdot*(Msun/yr)*2.3*mp*np.sqrt(G*Mstar*Msun/((r_grid*au)**3))/(3*np.pi*Sigma*kB*T_disk)
    return a

def h(rc, hc, psi): # H/R
    return hc*(r_grid/rc)**psi

def Bz2(H, Mstar, Sigma, deltag, rcav):
    Omega = np.sqrt(G*Mstar*Msun/((r_grid*au)**3))
    
    if deltag is None and rcav is None:
        return 4*np.sqrt(2*np.pi)*H*(Omega**2)*Sigma*(1e4**(-1))
    
    beta_array = np.where(r_grid < rcav,
                          1, deltag**(-1))
    Bz_2 = 4*np.sqrt(2*np.pi)*H*(Omega**2)*Sigma*(beta_array**(-2))
    return Bz_2

def BzBphi(Mstar, Mdot):
    Omega = np.sqrt(G*Mstar*Msun/((r_grid*au)**3))
    return Mdot*(Msun/yr)*Omega/(2*r_grid*au)

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
psi_list    = cavity['psi']
hc_list     = cavity['hc']

###############################################################################
# T plot
fig, ax = plt.subplots(6, 4, figsize=(25, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    
    T = temperature(psi=psi_list[i],
                    Lstar=Lstar_list[i])
    
    _, v = vr_profile(Mstar=Mstar_list[i],
                      Mdot=Mdot_list[i],
                      Sigma=sigma,
                      T=T)
    
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

###############################################################################
# alpha plot
fig, ax = plt.subplots(6, 4, figsize=(25, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    alpha_plot = alpha_profile(Mstar=Mstar_list[i],
                               Mdot=Mdot_list[i],
                               Sigma=sigma,
                               psi=psi_list[i],
                               Lstar=Lstar_list[i])
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, sigma, label = r'$\Sigma_{g}$')
    ax2.plot(r_grid, alpha_plot, label = r'$\alpha$', color = 'r')
    ax2.plot(r_grid, 1e-2*np.ones(r_grid.shape),'r:')
    ax2.plot(r_grid, np.ones(r_grid.shape),'r--')
    
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

###############################################################################
# vr plot
fig, ax = plt.subplots(6, 4, figsize=(25, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    T = temperature(psi=psi_list[i],
                    Lstar=Lstar_list[i])
    _, vr_plot = vr_profile(Mstar=Mstar_list[i],
                            Mdot=Mdot_list[i],
                            Sigma=sigma,
                            T=T)
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, sigma, label = r'$\Sigma_{g}$')
    ax2.plot(r_grid, vr_plot, label = r'$v_{r}$', color = 'r')

    ax2.plot(r_grid, np.ones(r_grid.shape),'r:')
    
    ax1.set_title(name_list[i], fontsize = 16)
    ax1.set_xlim((0.1, 200))
    ax1.set_ylim((1e-6, 1e6))
    ax2.set_ylim((1e-6, 1e2))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
    
    ax1.set_xlabel('R [AU]', fontsize = 16)
    ax1.set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax2.set_ylabel(r'$v_{r} [c_{s}]$', fontsize = 16)
    ax1.legend(fontsize = 16, loc='upper left')
    ax2.legend(fontsize = 16, loc='upper right')
    
plt.tight_layout()
plt.savefig('vr_profile.pdf', transparent = True)
plt.close()

###############################################################################
# Bz plot
fig, ax = plt.subplots(6, 4, figsize=(25, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])

    if name_list[i] in without_cav:
        rc_list[i] = 100
        
    h_list = h(rc=rc_list[i],
               hc=hc_list[i],
               psi=psi_list[i])
    H_list = h_list*r_grid*au
    Bz_squre = Bz2(H=H_list,
                   Mstar=Mstar_list[i],
                   Sigma=sigma,
                   deltag=deltag_list[i],
                   rcav=rcavg_list[i])
    BzBphi_profile = BzBphi(Mstar=Mstar_list[i],
                            Mdot=Mdot_list[i])
    
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, sigma, label = r'$\Sigma_{g}$')
    ax2.plot(r_grid, np.sqrt(Bz_squre), 'r--', label = r'$\sqrt{B_{z}^{2}}$')
    ax2.plot(r_grid, np.sqrt(BzBphi_profile), 'r:', label = r'$\sqrt{B_{z}B_{\phi}}$')
    
    ax1.set_title(name_list[i], fontsize = 16)
    ax1.set_xlim((0.1, 200))
    ax1.set_ylim((1e-6, 1e6))
    ax2.set_ylim((1e-8, 1e4))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
    
    ax1.set_xlabel('R [AU]', fontsize = 16)
    ax1.set_ylabel(r'$\Sigma_{g} [g cm^{-2}]$', fontsize = 16)
    ax2.set_ylabel(r'$[G]$', fontsize = 16)
    ax1.legend(fontsize = 16, loc='upper left')
    ax2.legend(fontsize = 16, loc='upper right')
    
plt.tight_layout()
plt.savefig('B_profile.pdf', transparent = True)
plt.close()

###############################################################################
# alpha vs B
fig, ax = plt.subplots(6, 4, figsize=(25, 30), sharex=True, sharey=True)

ax = ax.flatten()

for i in range(len(rc_list)):
    sigma = surface_density(rc = rc_list[i],
                            sigmac=sigmac_list[i],
                            delta_gas=deltag_list[i],
                            rcav=rcavg_list[i],
                            gamma=gamma_list[i])
    alpha_plot = alpha_profile(Mstar=Mstar_list[i],
                               Mdot=Mdot_list[i],
                               Sigma=sigma,
                               psi=psi_list[i],
                               Lstar=Lstar_list[i])
    
    
    if name_list[i] in without_cav:
        rc_list[i] = 100
        
    h_list = h(rc=rc_list[i],
               hc=hc_list[i],
               psi=psi_list[i])
    H_list = h_list*r_grid*au
    Bz_squre = Bz2(H=H_list,
                   Mstar=Mstar_list[i],
                   Sigma=sigma,
                   deltag=deltag_list[i],
                   rcav=rcavg_list[i])
    BzBphi_profile = BzBphi(Mstar=Mstar_list[i],
                            Mdot=Mdot_list[i])
    
    ax1 = ax[i]
    ax2 = ax1.twinx()
    
    ax1.plot(r_grid, alpha_plot, label = r'$\alpha$')
    ax2.plot(r_grid, np.sqrt(Bz_squre), 'r--', label = r'$\sqrt{B_{z}^{2}}$')
    ax2.plot(r_grid, np.sqrt(BzBphi_profile), 'r:', label = r'$\sqrt{B_{z}B_{\phi}}$')
    
    ax1.set_title(name_list[i], fontsize = 16)
    ax1.set_xlim((0.1, 200))
    ax1.set_ylim((1e-5, 1e4))
    ax2.set_ylim((1e-8, 1e4))
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    
    
    ax1.set_xlabel('R [AU]', fontsize = 16)
    ax1.set_ylabel(r'$\alpha$', fontsize = 16)
    ax2.set_ylabel(r'$[G]$', fontsize = 16)
    ax1.legend(fontsize = 16, loc='upper left')
    ax2.legend(fontsize = 16, loc='upper right')
    
plt.tight_layout()
plt.savefig('alpha_vs_B_profile.pdf', transparent = True)
plt.close()