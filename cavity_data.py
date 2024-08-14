import numpy as np
au   = 1.496e13     # cm/au
SB   = 5.6704e-5    # g/s^3/K^4
mp   = 1.6726e-24   # g/mp
kB   = 1.3807e-16   # cm^2 g/s^2/K
Msun = 1.9884e33    # g/Msun
Lsun = 3.839e33     # erg/s
G    = 6.6743e-8    # cm^3/g/s^2
yr   = 3.156e7      # s/yr

#                  rc,  Σc,   rcg,    δg,        Mdot,  M*,  Age,     Md,   L*,  rcd,   δd,    γ,    ψ,   hc
all_target = {
    'SR21'      :[ 15,  400,   25,  1e-2,  10**(-7.9),   2,  1.5,   4000,   11,   56,  1e-3, 1.0, 0.15, 0.07], # spectroscopic binary/substellar companion
    'HD135344B' :[ 25,  120,   30,  2e-4,  10**(-7.4), 1.6,  5.7,   5000,  6.7,   52,  2e-4, 1.0, 0.05, 0.15],
    'LkCa15'    :[ 85,   34,   45,  1e-1,  10**(-8.4), 0.8,  1.5,  34600,  1.3,   76,  1e-4, 1.0, 0.04, 0.06], # substellar companion
    'RXJ1615'   :[115,   50,   20,  1e-4,  10**(-8.5), 0.8,  1.5, 157675,  0.9,   17,  1e-5, 1.0, 0.20, 0.04], # substellar companion
    'J1604'     :[ 60,   12,   30,  1e-5, 10**(-10.5), 1.2,   10,   6700,  0.7,   87,  1e-6, 1.0, 0.68, 0.07],
    'CQTau'     :[ 56,  2.5,   20,  1e-2,    10**(-7), 1.6,  8.9,   1158,   10,   20,  1e-2, 0.3, 0.05, 0.13],
    'DoAr44'    :[ 25,   60,   16,  1e-4,  10**(-8.3),   1,  1.5,    850,  1.9,   40,  1e-2, 1.0, 0.10, 0.10], # substellar companion
    'IRS48'     :[ 60,  0.5,   25,  1e-3,  10**(-8.4),   2,  1.5,    188, 17.8,   70,  1e-3, 1.0, 0.22, 0.14],
    'TWHya'     :[ 35,   30,    4,  1e-2,  10**(-8.9), 0.6,    9,  16650, 0.28,    3,  1e-2, 1.0, 0.30, 0.10],
    'HD169142'  :[100,  6.5,   56, 0.025,  10**(-8.7),   2,    9,   2700,    8,   24,  0.27, 1.0, 0.01, 0.07], # ψ=0
    'Sz91'      :[ 75,    7,   50,  1e-5,  10**(-8.7), 0.5,  2.5,    300,  0.2,   67, 1e-20, 1.0, 0.30, 0.10],
    'J16083070' :[ 50,   39,   60,  1e-4,  10**(-9.1), 1.5,  2.5,   1000,  1.8,   62, 1e-20, 1.0, 0.10, 0.10],
    'Sz111'     :[ 50, 1500,   45,  1e-2,  10**(-9.1), 0.5,  2.5,   1600,  0.2,   55, 1e-20, 1.0, 0.10, 0.08],
    'RYLup'     :[ 25,  200,   50,  1e-1,  10**(-8.2), 1.5,  2.5,   2300,   15,   27,  1e-1, 1.0, 0.15, 0.15],
    'Sz118'     :[ 25,  100,   40,  1e-3,  10**(-9.2),   1,  2.5,    700,  0.7,   64,  1e-1, 1.0, 0.10, 0.07],
    'Sz123A'    :[ 25,   36,   30,  1e-3,  10**(-9.2), 0.6,  2.5,    400,  0.1,   39,  1e-2, 1.0, 0.15, 0.15],
    'HD100546'  :[ 30,  250,   15,  1e-5,    10**(-7), 2.2,  5.5,  16650,   25,   25, 1e-10, 1.1, 0.01, 0.10], # ψ=0
    # 'HD139614'  :[  6, 0.17,    6,  1e-2,    10**(-8), 1.8,  8.8,   3330,   11,    6,  1e-4, 1.0, None],
    'HD142527'  :[200,   22,   90,  2e-2,  10**(-7.5), 2.3,  6.6,  33300,  9.9,  185,  3e-5, 1.0, 0.01, 0.10], # binary/ ψ=0
    'DMTau'     :[124, 0.65,   21,  0.15,  10**(-8.3), 0.3,  1.5,   2000,  0.2,   18,  6e-9, 1.0, 0.20, 0.06], # tentative binary/substellar companion
    'IMLup'     :[100, 28.4, None,  None,  10**(-7.9), 1.1,  2.5,  66600, 2.57, None,  None, 1.0, 0.17, 0.10],
    'GMAur'     :[176,  9.4, None,  None,  10**(-8.1), 1.1,  1.5,  66600,    1, None,  None, 1.0, 0.35, 0.08],
    'AS209'     :[ 80,  1.0, None,  None,  10**(-7.3), 1.2,  1.5,   1500, 1.41, None,  None, 1.0, 0.25, 0.06],
    'HD163296'  :[165,  8.8, None,  None,  10**(-7.4), 2.0,    6,  46620,   17, None,  None, 0.8, 0.08, 0.08],
    'MWC480'    :[200,  5.8, None,  None,  10**(-6.9), 2.1,  6.5,  53280, 21.9, None,  None, 1.0, 0.08, 0.10]
}

without_cav       = [key for key in all_target.keys() if all_target[key][9] is None]
with_cav          = [key for key in all_target.keys() if all_target[key][9] is not None]
with_companion    = ['SR21','LkCa15','RXJ1615','DoAr44','HD142527','DMTau']
without_companion = [key for key in all_target.keys() if key not in with_companion]


def alpha_1au(psi, Lstar, Mstar, Mdot, Sigma):
    T_disk = (psi*Lstar*Lsun/(8*np.pi*SB*au**2))**(1/4)
    a = Mdot*(Msun/yr)*2.3*mp*np.sqrt(G*Mstar*Msun/(au**3))/(3*np.pi*Sigma*kB*T_disk)
    return a

def vr_1au(Mstar, Mdot, Sigma):
    vr_cms  = Mdot*(Msun/yr)/(2*np.pi*Sigma*au)
    kep_cms = np.sqrt(G*Mstar*Msun/au)
    return vr_cms/kep_cms

def sigma_at_1au(rc, sigmac, delta_gas, rcav, gamma = 1):
    if delta_gas is None and rcav is None:
        delta_gas, rcav = 1, 0
    density = delta_gas * sigmac * ((1/rc)**(-gamma)) * np.exp(-(1/rc)**(2-gamma))
    return density



def make_dict(target = with_cav + without_cav):
    
    cavity_rc        = []    # AU
    cavity_sigmac    = []    # gcm-2
    cavity_rcavg     = []    # AU
    cavity_deltag    = []    # 
    cavity_mdot      = []    # M_sun/yr
    cavity_mstar     = []    # M_sun
    cavity_age       = []    # Myr
    cavity_md        = []    # M_earth
    cavity_lum       = []    # L_sun
    cavity_rcavd     = []    # AU
    cavity_deltad    = []    # 
    cavity_existence = []    #
    cavity_companion = []    #
    cavity_vr        = []    # kep
    cavity_alpha     = []    #
    cavity_gamma     = []    #
    cavity_psi       = []    #
    cavity_disk_life = []    # yr
    cavity_hc        = []
    
    for t in target:
        cavity_rc.append(all_target[t][0])
        cavity_sigmac.append(all_target[t][1])
        cavity_rcavg.append(all_target[t][2])
        cavity_deltag.append(all_target[t][3])
        cavity_mdot.append(all_target[t][4])
        cavity_mstar.append(all_target[t][5])
        cavity_age.append(all_target[t][6])
        cavity_md.append(all_target[t][7])
        cavity_lum.append(all_target[t][8])
        cavity_rcavd.append(all_target[t][9])
        cavity_deltad.append(all_target[t][10])
        cavity_gamma.append(all_target[t][11])
        cavity_psi.append(all_target[t][12])
        cavity_disk_life.append(all_target[t][7]/(all_target[t][4]*333000))
        cavity_hc.append(all_target[t][13])
        
        if all_target[t][9] is not None:
            cavity_existence.append(True)
        else:
            cavity_existence.append(False)
        if t in with_companion:
            cavity_companion.append(True)
        else:
            cavity_companion.append(False)
        cavity_vr.append(
            vr_1au(all_target[t][5],
                   all_target[t][4],
                   sigma_at_1au(
                    all_target[t][0],
                    all_target[t][1],
                    all_target[t][3],
                    all_target[t][2],
                    all_target[t][11]
                               )
                   )
            )
        cavity_alpha.append(
            alpha_1au(all_target[t][12],
                      all_target[t][8],
                      all_target[t][5],
                      all_target[t][4],
                      sigma_at_1au(
                            all_target[t][0],
                            all_target[t][1],
                            all_target[t][3],
                            all_target[t][2],
                            all_target[t][11]
                                  )
                     )
            )
        
        
    cavity = {
        'name'      : target,
        'Mdot'      : cavity_mdot,
        'rc'        : cavity_rc,
        'sigmac'    : cavity_sigmac,
        'rcavg'     : cavity_rcavg,
        'rcavd'     : cavity_rcavd,
        'deltag'    : cavity_deltag,
        'deltad'    : cavity_deltad,
        'Mstar'     : cavity_mstar,
        'age'       : cavity_age,
        'Md'        : cavity_md,
        'Lstar'     : cavity_lum,
        'existence' : cavity_existence,
        'companion' : cavity_companion,
        'vr'        : cavity_vr,
        'alpha'     : cavity_alpha,
        'gamma'     : cavity_gamma,
        'psi'       : cavity_psi,
        'td'        : cavity_disk_life,
        'hc'        : cavity_hc
    }
    return cavity

key_to_text = {
    'Mdot'   : r'$\dot{M}$',
    'rc'     : r'$r_{c}$',
    'sigmac' : r'$\Sigma_{c}$',
    'rcavg'  : r'$r_{cav,g}$',
    'rcavd'  : r'$r_{cav,d}$',
    'deltag' : r'$\delta_{gas}$',
    'deltad' : r'$\delta_{dust}$',
    'Mstar'  : r'$M_{*}$',
    'age'    : r'$Age$',
    'Md'     : r'$M_{d}$',
    'Lstar'  : r'$L_{*}$',
    'vr'     : r'$v_{r}$',
    'alpha'  : r'$\alpha_{@1AU}$',
    'gamma'  : r'$\gamma$',
    'psi'    : r'$\psi$',
    'td'     : r'$t_{d}$'
}

key_to_unit = {
    'Mdot'   : r'$[M_{\odot}/yr]$',
    'rc'     : r'$[AU]$',
    'sigmac' : r'$[gcm^{-2}]$',
    'rcavg'  : r'$[AU]$',
    'rcavd'  : r'$[AU]$',
    'deltag' : '',
    'deltad' : '',
    'Mstar'  : r'$[M_{\odot}]$',
    'age'    : r'$[Myr]$',
    'Md'     : r'$[M_{\odot}]$',
    'Lstar'  : r'$[L_{\odot}$]',
    'vr'     : r'$[c_{s}]$',
    'alpha'  : '',
    'gamma'  : '',
    'psi'    : '',
    'td'     : r'$[Myr]$'
}

