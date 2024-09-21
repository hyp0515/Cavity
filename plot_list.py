import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
from scipy.optimize import curve_fit

from cavity_data import *
###############################################################################
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

def Bz2(Mdot, Mstar, bphibz_ratio=10):
    Omega = np.sqrt(G*Mstar*Msun/((r_grid*au)**3))
    return  Mdot*(Msun/yr)*Omega/bphibz_ratio/(2*r_grid*au)

def beta(Sigma, Mstar, H, Bz2):
    Omega = np.sqrt(G*Mstar*Msun/((r_grid*au)**3))
    return (4*np.sqrt(2*np.pi)*H*(Omega**2)*Sigma)/Bz2
###############################################################################

cavity = make_dict(target=with_cav+without_cav)

name_list   = np.array(cavity['name'])
rcavg_list  = np.array(cavity['rcavg'])
Mdot_list   = np.array(cavity['Mdot'])
alpha_list  = np.array(cavity['alpha'])
rc_list     = np.array(cavity['rc'])
sigmac_list = np.array(cavity['sigmac'])
deltag_list = np.array(cavity['deltag'])
Lstar_list  = np.array(cavity['Lstar'])
Mstar_list  = np.array(cavity['Mstar'])
gamma_list  = np.array(cavity['gamma'])
psi_list    = np.array(cavity['psi'])
hc_list     = np.array(cavity['hc'])
td_list     = np.array(cavity['existence'])

sigma_list        = np.zeros((len(name_list), len(r_grid)))
alpha_list        = np.zeros((len(name_list), len(r_grid)))
beta_list         = np.zeros((len(name_list), len(r_grid)))
rcav_index        = np.zeros((len(name_list)))
alpha_rcav_list   = np.zeros((len(name_list)))
beta_rcav_list    = np.zeros((len(name_list)))
alpha_upper_error = np.zeros((len(name_list)))
alpha_lower_error = np.zeros((len(name_list)))
beta_upper_error  = np.zeros((len(name_list)))
beta_lower_error  = np.zeros((len(name_list)))

for i in range(len(name_list)):

    sigma_list[i, :]  = surface_density(
        rc=rc_list[i],
        sigmac=sigmac_list[i],
        delta_gas=deltag_list[i],
        rcav=rcavg_list[i],
        gamma=gamma_list[i]
    )
    
    try:
        rcav_index[i] = np.searchsorted(r_grid, rcavg_list[i]) - 1
    except:
        rcav_index[i] = np.searchsorted(r_grid, 10) - 1
    
    alpha_list[i, :] = alpha_profile(
        Mstar=Mstar_list[i],
        Mdot=Mdot_list[i],
        Sigma=sigma_list[i, :],
        psi=psi_list[i],
        Lstar=Lstar_list[i]
    )
    alpha_rcav_list[i] = alpha_list[i, int(rcav_index[i])]
    alpha_upper_error[i] = np.max(alpha_list[i, :int(rcav_index[i])]) - alpha_rcav_list[i]
    alpha_lower_error[i] = alpha_rcav_list[i] - np.min(alpha_list[i, :int(rcav_index[i])])
    if alpha_upper_error[i] < 0: alpha_upper_error[i] = 0
    if alpha_lower_error[i] < 0: alpha_lower_error[i] = 0
    
    
    if td_list[i] is False:
        rc_list[i] = 100
    
        
    h_value = h(rc=rc_list[i],
               hc=hc_list[i],
               psi=psi_list[i])
    H_value = h_value*r_grid*au
    Bz_squre = Bz2(Mdot=Mdot_list[i],
                   Mstar=Mstar_list[i])
    beta_list[i, :] = beta(
        Sigma=sigma_list[i, :],
        Mstar=Mstar_list[i],
        H=H_value,
        Bz2=Bz_squre
    )
    beta_rcav_list[i] = beta_list[i, int(rcav_index[i])]
    beta_upper_error[i] = np.max(beta_list[i, :int(rcav_index[i])]) - beta_rcav_list[i]
    beta_lower_error[i] = beta_rcav_list[i] - np.min(beta_list[i, :int(rcav_index[i])])
    if beta_upper_error[i] < 0: beta_upper_error[i] = 0
    if beta_lower_error[i] < 0: beta_lower_error[i] = 0
    



sorted_indices = np.argsort(alpha_rcav_list)

td_list            = td_list[sorted_indices]
Mdot_list          = Mdot_list[sorted_indices]
rcavg_list         = rcavg_list[sorted_indices]
name_list          = name_list[sorted_indices]
alpha_rcav_list    = alpha_rcav_list[sorted_indices]
alpha_upper_error  = alpha_upper_error[sorted_indices]
alpha_lower_error  = alpha_lower_error[sorted_indices]
beta_rcav_list     = beta_rcav_list[sorted_indices]
beta_upper_error   = beta_upper_error[sorted_indices]
beta_lower_error   = beta_lower_error[sorted_indices]

alpha_yerr         = [alpha_lower_error, alpha_upper_error]
beta_yerr          = [beta_lower_error, beta_upper_error]

fig, ax = plt.subplots(4, 1, figsize=(12, 12), sharex=True)

# for i in td_list:
#     if i is True:
#         ax[0].scatter(range(len(name_list)), Mdot_list, marker='o')
#     else:
#         ax[0].scatter(range(len(name_list)), Mdot_list, marker='v')
ax[0].scatter(range(len(name_list)), Mdot_list, marker='o')
# ax[0].set_xticks(range(len(name_list)))
# ax[0].set_xticklabels(name_list, rotation=90, fontsize=10)
ax[0].set_yscale('log')
ax[0].set_ylabel('Accretion Rate', fontsize=12)
ax[0].grid(True, which="major", ls=":", linewidth=0.75)


# for i in td_list:
#     if i is True:
#         ax[1].scatter(range(len(name_list)), rcavg_list, marker='o')
#     else:
#         ax[1].scatter(range(len(name_list)), rcavg_list, marker='v')
ax[1].scatter(range(len(name_list)), rcavg_list)
ax[1].set_xticks(range(len(name_list)))
ax[1].set_xticklabels(name_list, rotation=90, fontsize=10)
ax[1].set_yscale('linear')
ax[1].set_ylabel('Cavity Size', fontsize=12)
ax[1].grid(True, which="both", ls=":", linewidth=0.75)


ax[2].errorbar(range(len(name_list)), alpha_rcav_list,
                yerr=alpha_yerr, fmt='o', ecolor='blue', capsize=5, elinewidth=1)

ax[2].set_xticks(range(len(name_list)))
ax[2].set_xticklabels(name_list, rotation=90, fontsize=10)
ax[2].set_yscale('log')
ax[2].set_ylabel(r'$\alpha$', fontsize=12)
ax[2].grid(True, which="major", ls=":", linewidth=0.75)
ax[2].axhline(y=1, color='red', linestyle='-')
ax[2].axhline(y=1e-2, color='red', linestyle='--')



ax[3].errorbar(range(len(name_list)), beta_rcav_list,
                yerr=beta_yerr, fmt='o', ecolor='blue', capsize=5, elinewidth=1)

ax[3].set_xticks(range(len(name_list)))
ax[3].set_xticklabels(name_list, rotation=90, fontsize=10)
ax[3].set_yscale('log')
ax[3].set_ylabel(r'$\beta$', fontsize=12)
ax[3].grid(True, which="major", ls=":", linewidth=0.75)
ax[3].axhline(y=1e2, color='red', linestyle='-')
# ax[3].axhline(y=1e-2, color='red', linestyle='--')







plt.tight_layout()
plt.savefig('test.pdf', transparent=True)
plt.show()
