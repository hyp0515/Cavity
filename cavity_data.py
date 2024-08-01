#                  rc,  Σc, rcg,   δg,      Mdot,     M*  Age,   Md,   L*, rcd
target = {
    'SR21'      :[ 15,  400, 25,  1e-2,  10**(-7.9),   2, 1.5,  1900,   11, 56],
    'HD135344B' :[ 25,  120, 20,  2e-4,  10**(-7.4), 1.6, 5.7,  4900,  6.7, 52],
    'LkCa15'    :[ 85,   34, 45,  1e-1,  10**(-8.4), 0.8, 1.5,  3300,  1.3, 76],
    'RXJ1615'   :[115,   50, 20,  1e-4,  10**(-8.5), 0.8, 1.5,  6500,  0.9, 17],
    'J1604'     :[ 60,   12, 70,  1e-5, 10**(-10.5), 1.2,  10,  2100,  0.7, 87],
    'CQTau'     :[ 56,  2.5, 20,  1e-2,    10**(-7), 1.6, 8.9,  4100,   10, 50],
    'DoAr44'    :[ 25,   60, 16,  1e-4,  10**(-8.3),   1, 1.5,  1900,  1.9, 40],
    'IRS48'     :[ 60,  0.5, 25,  1e-3,  10**(-8.4),   2, 1.5,  1300, 17.8, 70],
    'TWHya'     :[ 35,   30,  4,  1e-2,  10**(-8.9), 0.6,   9, 16650, 0.28, 3],
    'HD169142'  :[100,  6.5, 35, 0.025,  10**(-8.7),   2,   9,  2700,    8, 24],
    'Sz91'      :[ 75,    7, 50,  1e-5,  10**(-8.7), 0.5, 2.5,   300,  0.2, 67],
    'J16083070' :[ 50,   39, 60,  1e-4,  10**(-9.1), 1.5, 2.5,  1000,  1.8, 62],
    'Sz111'     :[ 50, 1500, 45,  1e-2,  10**(-9.1), 0.5, 2.5,  1600,  0.2, 55],
    'RYLup'     :[ 25,  200, 50,  1e-1,  10**(-8.2), 1.5, 2.5,  2300,   15, 27],
    'Sz118'     :[ 25,  100, 40,  1e-3,  10**(-9.2),   1, 2.5,   700,  0.7, 64],
    'Sz123A'    :[ 25,   36, 30,  1e-3,  10**(-9.2), 0.6, 2.5,   400,  0.1, 39],
    'HD100546'  :[ 30,  250, 15,  1e-5,    10**(-7), 2.2, 5.5,  4800,   25, 25],
    'HD142527'  :[200,   22, 90,  2e-2,  10**(-7.5), 2.3, 6.6, 30600,  9.9, 185],
    'HD139614'  :[  6, 0.17,  6,  1e-2,    10**(-8), 1.8, 8.8,  3330,   11, 6]
}

cavity_rc       = []    # AU
cavity_sigmac   = []    # gcm-2
cavity_rcav     = []    # AU
cavity_deltagas = []    # 
cavity_mdot     = []    # M_sun/yr
cavity_mstar    = []    # M_sun
cavity_age      = []    # Myr
cavity_md       = []    # M_earth
cavity_lum      = []    # L_sun
cavity_rcavd    = []    # AU

for t in target:
    cavity_rc.append(target[t][0])
    cavity_sigmac.append(target[t][1])
    cavity_rcav.append(target[t][2])
    cavity_deltagas.append(target[t][3])
    cavity_mdot.append(target[t][4])
    cavity_mstar.append(target[t][5])
    cavity_age.append(target[t][6])
    cavity_md.append(target[t][7])
    cavity_lum.append(target[t][8])
    cavity_rcavd.append(target[t][9])


cavity = {
    r'$\dot{M}$'      : cavity_mdot,
    r'$r_{c}$'        : cavity_rc,
    r'$\Sigma_{c}$'   : cavity_sigmac,
    r'$r_{cav}$'      : cavity_rcav,
    r'$\delta_{gas}$' : cavity_deltagas,
    r'$M_{*}$'        : cavity_mstar,
    r'$Myr$'          : cavity_age,
    r'$M_{d}$'        : cavity_md,
    r'$L_{*}$'        : cavity_lum,
    r'$r_{cav,d}$'    : cavity_rcavd,
}