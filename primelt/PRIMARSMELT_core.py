import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


def batch_melt(mgo_s, feo_s, f):
    mgo_l = np.arange(2, 30.1, 0.1)
    D_mgo = ((mgo_s / mgo_l) - f) / (1 - f)
    kd = 0.35
    D_feo = kd * D_mgo
    feo_l = feo_s / (f * (1 - D_feo) + D_feo)

    mgo_r = (mgo_s - f * mgo_l) / (1 - f)
    feo_r = (feo_s - f * feo_l) / (1 - f)

    return mgo_l, feo_l, mgo_r, feo_r


def afm_melt(mgo_s, feo_s, f):
    mgo_l, feo_l, mgo_r, feo_r = batch_melt(mgo_s, feo_s, f)
    D_mgo = mgo_r / mgo_l
    D_feo = feo_r / feo_l

    feo_l_afm = feo_s * (1 - (1 - f) ** (1 / D_feo)) / f
    mgo_l_afm = mgo_s * (1 - (1 - f) ** (1 / D_mgo)) / f

    return mgo_l_afm, feo_l_afm


def olivine_line(mgo_s, feo_s, mgo_l, feo_l, f):
    # This calculates olivine line
    if f != 0:
        kd_ol = 0.35

        fe_mg_l = 0.561 * feo_l / mgo_l
        fe_mg_ol = fe_mg_l * kd_ol
        mg_n_ol = 1 / (1 + fe_mg_ol)
        mg_ol = 2 * mg_n_ol
        fe_ol = 2 - mg_ol

        feo_ol = fe_ol * 71.844
        mgo_ol = mg_ol * 40.3044
        sio2_ol = 60.08
        sum_ol = feo_ol + mgo_ol + sio2_ol

        feo_ol = 100 * feo_ol / sum_ol
        mgo_ol = 100 * mgo_ol / sum_ol
        sio2_ol = 100 * sio2_ol / sum_ol

        _mgo_l = (mgo_s - (1 - f) * mgo_ol) / f
        _feo_l = (feo_s - (1 - f) * feo_ol) / f

        _kd_ol = (feo_ol / _feo_l) / (mgo_ol / _mgo_l)

        er = abs(kd_ol - _kd_ol)
        ind = np.where(er == np.min(er))
        ind = ind[0]
        ind = ind[0]
        return _mgo_l[ind], _feo_l[ind]


def make_figure_batch(ax1):
    mgo_s = 30.2
    feo_s = 17.9   ### XFo = 0.75
    #feo_s = 13.45   ### XFo = 0.80
    #feo_s = 15.18  ### XFo = 0.80

    for f in np.linspace(0, 0.6, 7):
        if f == 0:
            a, b, e, f = batch_melt(mgo_s, feo_s, f)
            ax1.plot(a, b, color='k', lw=1.5, label='Solidus')
        else:
            a, b, e, f = batch_melt(mgo_s, feo_s, f)
            ax1.plot(a, b, color='gray', lw=0.2)

    mg_ol = []
    fe_ol = []
    for f2 in np.linspace(0, 0.99, 100):
        if f2 > 0:
            a, b, e, f = batch_melt(mgo_s, feo_s, f2)
            c, d = olivine_line(mgo_s, feo_s, a, b, f2)
            if c >= 0 and d >= 0:
                mg_ol.append(c)
                fe_ol.append(d)

    ax1.fill_between(mg_ol, fe_ol, color='w', zorder=2)
    ax1.plot(mg_ol, fe_ol, color='k', lw=1, zorder=5, label='Liq + Ol')
    ax1.plot(mgo_s, feo_s, 'o', c='g', zorder=6, markersize=20, label='Peridotite Source')
    ax1.set(ylim=[5, 25], xlim=[2, 31], ylabel='FeO (wt%)', xlabel='MgO (wt%)', title='BATCH MELTING')
    ax1.legend(loc="lower right")


def make_figure_afm(ax2):
    mgo_s = 30.2
    feo_s = 17.9

    for f in np.linspace(0, 0.8, 9):
        if f == 0:
            a, b, c, d = batch_melt(mgo_s, feo_s, 0)
            ax2.plot(a, b, color='k', lw=1.5, label='Solidus')
        else:
            a, b = afm_melt(mgo_s, feo_s, f)
            ax2.plot(a, b, color='gray', lw=0.2)

    mg_ol = []
    fe_ol = []
    for f2 in np.linspace(0, 0.96, 96):
        if f2 > 0:
            a, b = afm_melt(mgo_s, feo_s, f2)
            c, d = olivine_line(mgo_s, feo_s, a, b, f2)
            if c >= 0 and d >= 0:
                mg_ol.append(c)
                fe_ol.append(d)

    ax2.fill_between(mg_ol, fe_ol, color='w', zorder=2)
    ax2.plot(mg_ol, fe_ol, color='k', lw=1, zorder=5, label='L + Ol')
    ax2.plot(mgo_s, feo_s, 'o', c='g', zorder=6, markersize=20, label='Peridotite Source')
    ax2.set(ylim=[5, 25], xlim=[2, 31], ylabel='FeO (wt%)', xlabel='MgO (wt%)', title='AFM MELTING')
    ax2.legend(loc="lower right")



fig = plt.figure(figsize=(18, 6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
make_figure_batch(ax1)
make_figure_afm(ax2)

##### HERE ARE CORE CALCULATIONS #####

# Inside your core module

class ConfigManager:
    def __init__(self):
        self.config = {}

    def update_settings(self, **kwargs):
        self.config.update(kwargs)

    def get_setting(self, key, default=None):
        return self.config.get(key, default)

# Initialize a global config manager instance
config_manager = ConfigManager()

mgo_s = config_manager.get_setting('mgo_s')
feo_s = config_manager.get_setting('feo_s')

#mgo_s = 30.2
#feo_s = 17.9   ### XFo = 0.75
#mgo_s = 31.0
#feo_s = 14.7  ### XFo = 0.79
#eo_s = 17.3  ### XFo = 0.78
#mgo_s = 32.8  ### XFo = 0.77 + 0.1 Fe2O3


def read_data(filename):
    # Read data
    data = pd.read_csv(filename, sep="\t|;|,", header=0, decimal='.', na_filter='-', engine='python')
    data.replace('0', 0)
    data.fillna(0, inplace=True)
    data.insert(13, 'H2O', np.zeros(len(data)))   #### add 0% H2O
    print('read_data OK')
    return data.apply(pd.to_numeric, errors='ignore')


def normalize(data):
    # Normalize data and calculate Fe2O3

    x = data['Fe2Fet']
    x = 1 - x
    data = data.drop('Fe2Fet', axis=1)
    y = data.sum(axis=1)
    data = 100 * data / y.values[0]
    z = x * data['FeO'] / 0.8998
    data.insert(4, 'Fe2O3', z)
    data['FeO'] = data['FeO'] * (1 - x)
    y = data.sum(axis=1)
    data = 100 * data / y.values[0]
    return data


def test_data(data):
    while True:
        try:
            data = data.drop('Name', axis=1)
            data = data.drop('P', axis=1)
            data = normalize(data)
        except:
            return False
        else:
            return True


def mol_cat_oxides(data):
    mw = np.asarray([60.0848, 79.8988, 101.96128, 151.9902, 159.6922, 71.8464, 70.9374, 40.3044, 56.0794, 61.97894, 94.1954, 74.7094, 141.94452, 18.0152])
    data_mol = data / mw
    y = data_mol.sum(axis=1)
    data_mol = 100 * data_mol / y.values[0]

    cp = np.asarray(
        [60.0848, 79.8988, 101.96128 / 2, 151.9902 / 2, 159.6922 / 2, 71.8464, 70.9374, 40.3044, 56.0794, 61.97894 / 2, 94.1954 / 2,
         74.7094, 141.94452 / 2, 18.0152 / 2])
    data_cat = data / cp
    data_cat = data_cat.set_axis(['Si', 'Ti', 'Al', 'Cr', 'Fe3', 'Fe', 'Mn', 'Mg', 'Ca', 'Na', 'K', 'Ni', 'P', 'H'],
                                 axis=1, copy=False)
    z = data_cat.sum(axis=1)
    data_cat = data_cat / z.values[0]

    return data_mol, data_cat


def temp_olivine(cat, p):
    #R = 8.31446261815324
    R = 8.3143
    dh = 113100
    ds = 52.05
    dv = 0.00000411

    a_ol_liq_m = [0, 0, 0, 0, 0.279, 0.259, 1, 0.0056, 0, 0, 3.346, 0]
    b_ol_liq_m = [0.0001, 0.0001, 0.0001, 0.0001, 0.031, -0.049, 0, 0.0135, 0.0001, 0.0001, -3.665, 0.0001]

    suma = np.sum(cat.loc[:, 'Ti':'P'].values * a_ol_liq_m, axis=1)
    sumb = np.sum(cat.loc[:, 'Ti':'P'].values * b_ol_liq_m, axis=1)
    d_mg = (2 / 3 - sumb) / suma
    d_mg = d_mg.astype(float)

    c_nm = np.sum(cat[['Fe', 'Mn', 'Mg', 'Ca', 'Ni']].values, axis=1)
    c_nm = c_nm.astype(float)
    c_si = np.sum(cat[['Si']].values, axis=1)
    c_si = c_si.astype(float)
    c_al = np.sum(cat[['Al']].values, axis=1)
    c_al = c_al.astype(float)
    c_ti = np.sum(cat[['Ti']].values, axis=1)
    c_ti = c_ti.astype(float)

    nf = 3.5 * np.log(1 - c_al) + 7 * np.log(1 - c_ti)
    t_ol = (dh / R + p * dv / R) / (ds / R + 2 * np.log(d_mg) + 2 * np.log(1.5 * c_nm) + 2 * np.log(3 * c_si) - nf)
    t_ol = t_ol - 273.15
    t_ol_a = t_ol + 54 * (p / 10000) - 2 * (p / 10000) ** 2
    t_ol_1b = t_ol + 54 * (1 / 10000) - 2 * (1 / 10000) ** 2
    print('olivine temperature OK!')
    return t_ol_a[0], t_ol_1b[0]

##### Have to add a new function that uses the Putirka (2008; Eq. 15) thermometer in Collinet et al. (2021) paper, instead of Beattie (1993)
def temp_olivine_p(data, p):
    #mg_number = (data.MgO/40.3044) / (data.MgO/40.3044 + data.FeO/71.844)
    #t_ol = 815.3 + 265.5*mg_number + 15.37*data.MgO + 8.61*data.FeO + 6.646*(data.Na2O + data.K2O) + 39.16*p/10000 - 12.83*data.H2O
    #t_ol_1b = 815.3 + 265.5*mg_number + 15.37*data.MgO + 8.61*data.FeO + 6.646*(data.Na2O + data.K2O) + 39.16/10000 - 12.83*data.H2O
    #a, b, c, d = 1.08327401e+03, 1.80750881e+01, 3.26805278e-02, 3.53941757e+01
    a, b, c, d = 1042.406169413418, 27.254916493709445, -0.3752746389512436, 36.92053769428407
    t_ol = a + b * data.MgO + c * data.MgO ** 2 + d * p/10000
    t_ol_1b = a + b * data.MgO + c * data.MgO ** 2 + d/10000
    return t_ol[0], t_ol_1b[0]


def sio2pound(comp, comp_mol):
    comp_anh = comp_mol[0]
    comp_anh = comp_anh.drop(columns=['H2O'])

    y = comp_anh.sum(axis=1)
    comp_anh = 100 * comp_anh / y.values[0]

    if comp_anh.SiO2[0] <= 60:
        psi = (0.46 * (100 / (100 - comp_anh.SiO2[0])) - 0.93) * (comp_anh.Na2O[0] + comp_anh.K2O[0]) + (
                -5.33 * (100 / (100 - comp_anh.SiO2[0])) + 9.69)
        sio2a = comp_anh.SiO2[0] + psi * (comp_anh.Na2O[0] + comp_anh.K2O[0])
        sio2n = sio2a + 0.8 * (comp.H2O[0])
    else:
        psi = (11 - 5.5 * (100 / (100 - comp_anh.SiO2[0]))) * np.exp(-0.13 * (comp_anh.Na2O[0] + comp_anh.K2O[0]))
        sio2a = comp_anh.SiO2[0] + psi * (comp_anh.Na2O[0] + comp_anh.K2O[0])
        sio2n = sio2a + 0.8 * (comp.H2O[0])

    return sio2n


def calc_kd(data, t_ol, sio2n, p):
    #R = 8.31446261815324
    R = 8.3143
    T = t_ol + 273.15
    sio2n = sio2n
    Fe = data.FeO[0] / 71.844
    Mg = data.MgO[0] / 40.3044
    xfo = 0.75
    kd_old = 0.35
    it = 0

    kd = np.exp((-6766 / (R * T) - 7.34 / R) + np.log(0.036 * sio2n - 0.22) + (3000 * (1 - 2 * xfo) / (R * T)) + (
            0.035 * (p - 1) / (R * T)))

    xfo = 1 / (kd * (Fe / Mg) + 1)

    while np.abs(kd - kd_old) > 1e-4:
        it += 1
        kd_old = kd
        kd = np.exp((-6766 / (R * T) - 7.34 / R) + np.log(0.036 * sio2n - 0.22) + (3000 * (1 - 2 * xfo) / (R * T)) + (
                0.035 * (p - 1) / (R * T)))
        xfo = 1 / (kd * (Fe / Mg) + 1)


    return kd, xfo


def eq_olivine(comp_mol, xfo):
    comp_mol = comp_mol[1]
    d_mgo = xfo / comp_mol.Mg[0]

    cat = comp_mol.loc[:, 'Ti':'P'].values
    coef = np.array([[0.0, 0.0035, 0.067, 0.0, 0.263, 0.214, 1.0, 0.0071, 0.0, 0.0, 3.346, 0.0],
                     [0.03, -0.031, 0.183, 0.0001, 0.196, 0.118, 0.0, -0.019, 0.0001, 0.0001, -3.665, 0.0001],
                     [0.0, 0.093, 0.0, 0.0, 0.0, 0.0, 0.0, 0.063, 0.0, 0.0, 0.0, 0.0]])

    b = np.array([d_mgo * (2 / 3), 1, 1 / (d_mgo * (2 / 3))])

    c = coef.T @ b.reshape(-1, 1)
    c = c.reshape(1, -1) * cat * 100
    c = c[0]

    sum_1 = 0.5 * (c[1] - c[2] - c[8] - c[9] - c[11]) + c[2] + c[3] + c[4] + c[5] + c[6] + c[7] + c[10] + 2 * (
            c[8] + c[9] + c[11])
    ox_ol = c * (200 / 3) / sum_1

    sum_2 = 0.5 * (ox_ol[1] - ox_ol[2] - ox_ol[8] - ox_ol[9] - ox_ol[11]) + ox_ol[2] + ox_ol[3] + ox_ol[5] + ox_ol[7] + \
            c[10] + 2 * (ox_ol[8] + ox_ol[9] + ox_ol[11])
    ox_ol[4] = (1 - xfo) * (200 / 3 - sum_2)
    ox_ol[6] = (xfo) * (200 / 3 - sum_2)
    ox_ol[10] = ((3.346 * ox_ol[6] / comp_mol.Mg[0] / 100) - 3.665) * 100 * comp_mol.Ni[0]
    sio2 = (100 / 3) - ox_ol[0] - ox_ol[2] - 0.5 * (ox_ol[1] - ox_ol[2] - ox_ol[8] - ox_ol[9] - ox_ol[11])
    h2o = 0
    ox_ol = np.concatenate([[sio2], ox_ol, [h2o]])
    cp = np.asarray(
        [60.0848, 79.8988, 101.96128 / 2, 151.9902 / 2, 159.6922 / 2, 71.8464, 70.9374, 40.3044, 56.0794, 61.97894 / 2, 94.1954 / 2,
         74.7094, 141.94452 / 2, 18.0152 / 2])

    z = ox_ol * cp
    ox_ol = 100 * ox_ol * cp / z.sum()

    ol_comp = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'NiO', 'P2O5',
                 'H2O'])
    ox_ol = pd.Series(ox_ol, index=ol_comp.columns)

    ol_comp = pd.concat([ol_comp, ox_ol.to_frame().T], ignore_index=True)
    print('equilibrium olivine OK!')
    return ol_comp


def f_values(comp_mol):
    coef = np.array([[0.0, -1.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, -0.5, -0.5, 3, 0.5, 0.0, 0.0],
                     [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     [1.0, 0.0, -0.5, -0.5, 0.0, -0.5, -0.5, -0.5, -1.5, -3.0, -3.0, -0.5, 0.0, 0.0]])

    comp_mol = comp_mol[0].values
    ol_an_q = coef @ comp_mol.T
    suma = ol_an_q.sum()
    ol_an_q = ol_an_q / suma
    ol_an_q = ol_an_q.astype(float)

    #a,b,c,d = -0.16493050706635692, 1.3920490156986904, 15.467495095918055, 0.056720656838951555
    # f1 = a + (np.exp(b * ol_an_q[1])) / (c * ol_an_q[1] + d)

    a,b,c,d = 1.5659541137528663, 0.22734661860160088, 32.90876426562045, -0.2877726247264274
    f1 = (a * np.exp(b*ol_an_q[1])) / (c*ol_an_q[1] + d)

    return f1, ol_an_q

def kd_f_values(comp, ol_an_q, mgo_s, feo_s):
    print(f"kd_f_values called with mgo_s={mgo_s}, feo_s={feo_s}")
    if mgo_s is None or feo_s is None:
        raise ValueError("mgo_s or feo_s is None")

    sio2 = comp.SiO2[0]
    al2o3 = comp.Al2O3[0]

    cr2o3 = comp.Cr2O3[0]
    cao = comp.CaO[0]
    mgo = comp.MgO[0]
    feo = comp.FeO[0]
    mg_number = (mgo / 40.3044) / ((mgo / 40.3044) + (feo / 71.844))
    an = ol_an_q[1]

    #a, b, c, d, e, f, g, h = -1.8595100e+01,  2.9239500e+00, -4.4504400e+00, -1.1953500e+00, 3.6635247e+02, -1.8340790e+02,  1.6530000e-02, 3.6848000e-01  # with no error in Kd
    #x = np.arctan(a * cr2o3 + b * feo + c * mgo + d * cao + e * mg_number + f) * np.array([1])
    # kd_new = g * x + h

    # Fe3 0.05
    a, b, c, d, e, f, g, h = -0.0019333795116506563, -0.0017160591478513774, -0.052740612773712926, 0.018191632074304014, -0.034321224555874795, -0.004130761603503998, 2.5225545104988276, -0.7399807560604785 ### Note coefficient g in manuscript was changed to 0.025  to make Mg# 100*Mg/(Mg+Fe)
    # All Fe2 #a,b,c,d,e,f,g,h = -0.0018534041125449394e+00, -0.0013590861979132563e+00, -0.04997584032631637e+00, 0.01629548012598762e+00, -0.0321809814003769e+00, -0.003802989185443774e+00, 2.376500496449525e+00, -0.7019405333477691e+00
    x = (a * sio2 + b * al2o3 + c * cr2o3 + d * feo + e * mgo + f * cao + g * mg_number + h) * np.array([1])

    # Fe3 0.05
    k, x0, LL, UL = 20.95149644752977, 0.34816490786696835, 0.2400000000382186, 0.4572915562664882
    # All Fe2 #k, x0, LL, UL = 26.993603729028518e+00, 0.32907539361670674e+00, 0.24002180043254692e+00, 0.4178953693020628e+00
    kd_new = LL + ((UL - LL) / (1 + np.exp(-k * ( x - x0 ))))

    fafm = 0
    f_femg = 0

    if feo == 0:
        f_femg = 0

    else:
        fb = (mgo * (feo_s / feo) - kd_new * mgo_s) / (mgo * (feo_s / feo) - mgo * kd_new)
        fa = 0
        fc = 0

        i = 0
        while np.abs(fc - fb) > 1e-4:
            fa = (mgo * (((feo_s - fb * feo) / (1 - fb)) / feo) - kd_new * mgo_s) / (
                    mgo * (((feo_s - fb * feo) / (1 - fb)) / feo) - mgo * kd_new)
            if fc != 0 and (fa - fb) / (fc - fb) > 0.9:
                fa = (fa + fb) / 2
            fc = fb
            fb = fa
            f_femg = fa
            i += 1

    if f_femg == 0:
        f_femg = 1e-5

    if f_femg >= 0:
        fafm = 1 / (0.98 / f_femg - 0.90 * (np.log(f_femg) ** 2) + 0.07 * (np.log(f_femg) ** 4))
    else:
       fafm = f_femg

    print('kd calculated OK!')
    return kd_new, f_femg, fafm


#def kd_f_values(comp, ol_an_q,  feo_s, mgo_s):

#    sio2 = comp.SiO2[0]
#    al2o3 = comp.Al2O3[0]
#
#    cr2o3 = comp.Cr2O3[0]
#    cao = comp.CaO[0]
#    mgo = comp.MgO[0]
#    feo = comp.FeO[0]
#    mg_number = (mgo / 40.3044) / ((mgo / 40.3044) + (feo / 71.844))
#    an = ol_an_q[1]
#
#    #a, b, c, d, e, f, g, h = -1.8595100e+01,  2.9239500e+00, -4.4504400e+00, -1.1953500e+00, 3.6635247e+02, -1.8340790e+02,  1.6530000e-02, 3.6848000e-01  # with no error in Kd
#    #x = np.arctan(a * cr2o3 + b * feo + c * mgo + d * cao + e * mg_number + f) * np.array([1])
#    # kd_new = g * x + h
#
#    # Fe3 0.05
#    a, b, c, d, e, f, g, h = -0.0019333795116506563, -0.0017160591478513774, -0.052740612773712926, 0.018191632074304014, -0.034321224555874795, -0.004130761603503998, 2.5225545104988276, -0.7399807560604785 ### Note coefficient g in manuscript was changed to 0.025  to make Mg# 100*Mg/(Mg+Fe)
#    # All Fe2 #a,b,c,d,e,f,g,h = -0.0018534041125449394e+00, -0.0013590861979132563e+00, -0.04997584032631637e+00, 0.01629548012598762e+00, -0.0321809814003769e+00, -0.003802989185443774e+00, 2.376500496449525e+00, -0.7019405333477691e+00
#    x = (a * sio2 + b * al2o3 + c * cr2o3 + d * feo + e * mgo + f * cao + g * mg_number + h) * np.array([1])
#
#    # Fe3 0.05
#    k, x0, LL, UL = 20.95149644752977, 0.34816490786696835, 0.2400000000382186, 0.4572915562664882
#    # All Fe2 #k, x0, LL, UL = 26.993603729028518e+00, 0.32907539361670674e+00, 0.24002180043254692e+00, 0.4178953693020628e+00
#    kd_new = LL + ((UL - LL) / (1 + np.exp(-k * ( x - x0 )))) * np.array([1])

#    try:
#        if feo == 0:
#            f_femg = 1e-05
#            fafm = 1e-05
#            return kd_new, f_femg, fafm
#        else:
#            fb = (mgo * (feo_s / feo) - kd_new * mgo_s) / (mgo * (feo_s / feo) - mgo * kd_new)
#            fa = 0.0
#            fc = 0.0
#            iter = 0
#            while (abs(fc - fb) > precision) and (iter < 30):
#                iter += 1
#                fa = (mgo * (((feo_s - fb * feo) / (1 - fb)) / feo) - kd_new * mgo_s) / \
#                     (mgo * (((feo_s - fb * feo) / (1 - fb)) / feo) - mgo * kd_new)
#                if fc != 0.0 and (fa - fb) / (fc - fb) > 0.9:
#                    fa = (fa + fb) / 2
#                fc = fb
#                fb = fa
#            result = fa
#            f_femg = max(result, 1e-05)
#            fafm = 1 / (0.98 / f_femg - 0.90 * (np.log(f_femg) ** 2) + 0.07 * (np.log(f_femg) ** 4))
#            return kd_new, f_femg, fafm
#    except:
#        f_femg = 1e-05
#        fafm = 1e-05
#        return kd_new, f_femg, fafm



def calculate(comp, ol_frac, ope, p):
    comp_mol = mol_cat_oxides(comp)
    #t_ol, t_ol_1b = temp_olivine(comp_mol[1], p)  ## -> Beattie 1993
    t_ol, t_ol_1b = temp_olivine_p(comp, p)
    sio2n = sio2pound(comp, comp_mol)
    kd, xfo = calc_kd(comp, t_ol, sio2n, p)
    comp_ol = eq_olivine(comp_mol, xfo)
    f1, ol_an_q = f_values(comp_mol)
    kd_new, f_femg, fafm = kd_f_values(comp, ol_an_q, mgo_s = config_manager.get_setting('mgo_s'), feo_s = config_manager.get_setting('feo_s'))

    if ope == 'add':
        new_comp = (comp + ol_frac * comp_ol) / (1 + ol_frac)
    if ope == 'sub':
        new_comp = (comp - ol_frac * comp_ol) / (1 - ol_frac)

    return comp, comp_ol, f1, f_femg, fafm, t_ol, t_ol_1b, kd, xfo, new_comp


def find_intersection(f_1, f_2, fc_1, fc_2):
    def make_line(p1, p2):
        A = (p1[1] - p2[1])
        B = (p2[0] - p1[0])
        C = (p1[0] * p2[1] - p2[0] * p1[1])
        return A, B, -C

    L1 = make_line((0, f_1), (1, f_2))
    L2 = make_line((0, fc_1), (1, fc_2))

    D = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        if x >= 0 and x <= 1:
            return x
        else:
            return False


def fractionate_olivine(comp, ol_frac, ope, typ, p):
    j = 0
    i = 0

    magma = pd.DataFrame(columns=comp.columns)
    olivine = pd.DataFrame(columns=comp.columns)

    t_olivine = []
    t_olivine_1b = []
    kd_value = []
    xfo_value = []
    ol_add = []

    x = 0
    sol = False

    if ope == 'add':
        while j <= 0.75 + ol_frac:
            melt, oliv, f1, f_femg, fafm, t_ol, t_ol_1b, kd, xfo, new_comp = calculate(comp, ol_frac, ope, p)
            _melt, _oliv, f1_2, f_femg_2, fafm_2, _t_ol, _t_ol_1b, _kd, _xfo, _new_comp = calculate(new_comp, ol_frac, ope, p)

            comp = new_comp
            magma = pd.concat([magma, melt.iloc[[0]]], ignore_index=True)
            olivine = pd.concat([olivine, oliv.iloc[[0]]], ignore_index=True)
            ol_add.append(j)

            if typ == 'afm':
                f_femg = fafm
                f_femg_2 = fafm_2

            if abs(f_femg) >= 1.0 and abs(f_femg_2 / f_femg) > 1:
                break

            check_f1 = find_intersection(f1, f1_2, f_femg, f_femg_2)

            t_olivine.append(t_ol)
            t_olivine_1b.append(t_ol_1b)
            kd_value.append(kd)
            xfo_value.append(xfo)

            if (check_f1 != False) and (f_femg >= 0):
                magma = pd.concat([magma, _melt.iloc[[0]]], ignore_index=True)
                olivine = pd.concat([olivine, _oliv.iloc[[0]]], ignore_index=True)
                ol_add.append(j + ol_frac)

                x = check_f1

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            j = j + ol_frac
            i += 1

    if ope == 'sub':
        while j <= 0.75 + ol_frac:
            melt, oliv, f1, f_femg, fafm, t_ol, t_ol_1b, kd, xfo, new_comp = calculate(comp, ol_frac, ope, p)
            _melt, _oliv, f1_2, f_femg_2, fafm_2, _t_ol, _t_ol_1b, _kd, _xfo, _new_comp = calculate(new_comp, ol_frac, ope, p)

            comp = new_comp
            #magma = magma.append(melt.iloc[0])
            #olivine = olivine.append(oliv.iloc[0])
            magma = pd.concat([magma, melt.iloc[[0]]], ignore_index=True)
            olivine = pd.concat([olivine, oliv.iloc[[0]]], ignore_index=True)
            ol_add.append(j)

            if typ == 'afm':
                f_femg = fafm
                f_femg_2 = fafm_2

            if abs(f_femg) >= 1.0 and abs(f_femg_2 / f_femg) > 1:
                break

            check_f1 = find_intersection(f1, f1_2, f_femg, f_femg_2)

            t_olivine.append(t_ol)
            t_olivine_1b.append(t_ol_1b)
            kd_value.append(kd)
            xfo_value.append(xfo)

            if (check_f1 != False) and (f_femg >= 0):
                #magma = magma.append(_melt.iloc[0])
                #olivine = olivine.append(_oliv.iloc[0])
                magma = pd.concat([magma, _melt.iloc[[0]]], ignore_index=True)
                olivine = pd.concat([olivine, _oliv.iloc[[0]]], ignore_index=True)
                ol_add.append(j + ol_frac)

                x = check_f1

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            j = j + ol_frac
            i += 1

    return magma, olivine, sol, ol_add


def run_all(data):
    ol_frac = 0.01
    sol_batch = []
    sol_afm = []
    ol_add_batch = []
    ol_add_afm = []
    fafm = []
    fbatch = []
    op_batch = []
    op_afm = []
    melt_sol_batch = melt_sol_afm = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'NiO', 'P2O5',
                 'H2O'])

    name = []
    pres = []
    input_error = []

    fig = plt.figure(figsize=(18, 6))
    ax_batch = fig.add_subplot(121)
    ax_afm = fig.add_subplot(122)
    make_figure_batch(ax_batch)
    make_figure_afm(ax_afm)

    for i in range(len(data)):
        sample = data.iloc[i]

        comp = pd.DataFrame(columns=data.columns)
        #comp = comp.append(sample)
        comp = pd.concat([comp, data.iloc[[i]]], ignore_index=True)
        comp = comp.reset_index(drop=True)

        name.append(comp.Name[0])
        p = comp.P[0]
        pres.append(p)

        px_error = ''
        vt_error = ''

        if (np.exp((comp.MgO[0]/40.3044)/(comp.CaO[0]/56.0774 + comp.FeO[0]/71.845)) - 1) / ((comp.MgO[0]/40.3044)/(comp.MgO[0]/40.3044 + comp.FeO[0]/71.845)) < 2.0 and ((comp.MgO[0]/40.3044)/(comp.MgO[0]/40.3044 + comp.FeO[0]/71.845)) < 0.4:
            px_error = ' PX '
        if comp.CaO[0] > 2.318 * comp.SiO2[0] - 93.626:
            vt_error = ' VT '
        if px_error != '' or vt_error != '':
            in_error = px_error + vt_error
        else:
            in_error = ' None '

        input_error.append(in_error)

        comp = comp.drop('Name', axis=1)
        comp = comp.drop('P', axis=1)
        comp = normalize(comp)

        # BATCH
        magma_batch, olivine_batch, sol, ol_add = fractionate_olivine(comp, ol_frac, 'add', 'femg', p)

        if sol == True:
            sol_batch.append('Y')
            op_batch.append(1)
            #melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])
            melt_sol_batch = pd.concat([melt_sol_batch, magma_batch.iloc[[len(magma_batch) - 1]]], ignore_index=True)
            ax_batch.plot(magma_batch.MgO, magma_batch.FeO, c='k', lw=0.5)
            ax_batch.plot(magma_batch['MgO'].iloc[-1], magma_batch['FeO'].iloc[-1], 'o', mfc='w', mec='k')
            ax_batch.plot(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0], 'o', mfc='k', mec='None')
            ax_batch.annotate(text=name[i], xy=(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0]),
                              xycoords='data', xytext=(-8, -4), textcoords='offset points', ha='center', va='top',
                              size='x-small')

        if sol == False:
            magma_batch, olivine, sol, ol_add = fractionate_olivine(comp, ol_frac, 'sub', 'femg', p)
            if sol == True:
                sol_batch.append('Y')
                op_batch.append(-1)
                #melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])
                melt_sol_batch = pd.concat([melt_sol_batch, magma_batch.iloc[[len(magma_batch) - 1]]], ignore_index=True)
                ax_batch.plot(magma_batch.MgO, magma_batch.FeO, c='k', lw=0.5)
                ax_batch.plot(magma_batch['MgO'].iloc[-1], magma_batch['FeO'].iloc[-1], 'o', mfc='w', mec='k')
                ax_batch.plot(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0], 'o', mfc='k', mec='None')
                ax_batch.annotate(text=name[i], xy=(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0]),
                                  xycoords='data', xytext=(8, 4), textcoords='offset points', ha='center', va='bottom',
                                  size='x-small')
            else:
                sol_batch.append('N')
                op_batch.append(0)
                #melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])
                melt_sol_batch = pd.concat([melt_sol_batch, magma_batch.iloc[[len(magma_batch) - 1]]], ignore_index=True)
        
        ol_add_batch.append(100 * ol_add[len(ol_add) - 1])

        # AFM
        magma_afm, olivine_afm, sol, ol_add = fractionate_olivine(comp, ol_frac, 'add', 'afm', p)

        if sol == True:
            sol_afm.append('Y')
            op_afm.append(1)
            #melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])
            melt_sol_afm = pd.concat([melt_sol_afm, magma_afm.iloc[[len(magma_afm) - 1]]], ignore_index=True)
            ax_afm.plot(magma_afm.MgO, magma_afm.FeO, c='k', lw=0.5)
            ax_afm.plot(magma_afm['MgO'].iloc[-1], magma_afm['FeO'].iloc[-1], 'o', mfc='w', mec='k')
            ax_afm.plot(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0], 'o', mfc='k', mec='None')
            ax_afm.annotate(text=name[i], xy=(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0]), xycoords='data',
                            xytext=(-8, -4), textcoords='offset points', ha='center', va='top', size='x-small')

        if sol == False:
            magma_afm, olivine_afm, sol, ol_add = fractionate_olivine(comp, ol_frac, 'sub', 'afm', p)
            if sol == True:
                sol_afm.append('Y')
                op_afm.append(-1)
                #melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])
                melt_sol_afm = pd.concat([melt_sol_afm, magma_afm.iloc[[len(magma_afm) - 1]]], ignore_index=True)
                ax_afm.plot(magma_afm.MgO, magma_afm.FeO, c='k', lw=0.5)
                ax_afm.plot(magma_afm['MgO'].iloc[-1], magma_afm['FeO'].iloc[-1], 'o', mfc='w', mec='k')
                ax_afm.plot(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0], 'o', mfc='k', mec='None')
                ax_afm.annotate(text=name[i], xy=(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0]), xycoords='data',
                                xytext=(8, 4), textcoords='offset points', ha='center', va='bottom', size='x-small')
            else:
                sol_afm.append('N')
                op_afm.append(0)
                #melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])
                melt_sol_afm = pd.concat([melt_sol_afm, magma_afm.iloc[[len(magma_afm) - 1]]], ignore_index=True)

        ol_add_afm.append(100 * ol_add[len(ol_add) - 1])

    melt_sol_batch.index = name
    melt_sol_batch.index.name = 'Name'
    melt_sol_batch.insert(0, '%Ol', [a * b for a, b in zip(ol_add_batch, op_batch)])
    melt_sol_batch.insert(0, 'P', pres)
    melt_sol_batch.insert(0, 'Solution', sol_batch)
    melt_sol_batch.insert(0, 'Warning', input_error)

    melt_sol_afm.index = name
    melt_sol_afm.index.name = 'Name'
    melt_sol_afm.insert(0, '%Ol', [a * b for a, b in zip(ol_add_afm, op_afm)])
    melt_sol_afm.insert(0, 'P', pres)
    melt_sol_afm.insert(0, 'Solution', sol_afm)
    melt_sol_afm.insert(0, 'Warning', input_error)


    return melt_sol_batch, melt_sol_afm, pres, input_error, fig


def pressure(afm_solution):
    Pi_ = []

    for i in range(len(afm_solution)):
        mgo = afm_solution.MgO[i]

        #pi = 0.93730275 + 0.1417003*mgo + 0.00463463*mgo**2
        pi = 0.08188644 + 0.21344162*mgo + 0.0028519*mgo**2

        Pi = pi

        Pi_.append(Pi)

    return Pi_


def final_results(data):
    batch_solution, afm_solution, p, input_error, figure = run_all(data)
    batch_sol = batch_solution.drop(['Solution', '%Ol', 'Warning', 'P'], axis=1)
    afm_sol = afm_solution.drop(['Solution', '%Ol', 'Warning', 'P'], axis=1)

    f_batch = []
    f_afm = []
    x_fo_batch = []
    x_fo_afm = []
    kd_batch = []
    kd_afm = []
    tol_batch = []
    tol_batch_1b = []
    tol_afm = []
    tol_afm_1b = []
    source = []
    f1_batch = []
    f1_afm = []

    for i in range(len(batch_sol)):
        melt, oliv, f1, f_femg, fafm, t_ol, t_ol_1b, kd, xfo, new_comp = calculate(batch_sol.iloc[[i]], 0, 'add',
                                                                          p[i])
        f_batch.append(f_femg[0])
        f1_batch.append(f1)
        x_fo_batch.append(xfo)
        kd_batch.append(kd)
        tol_batch.append(t_ol)
        tol_batch_1b.append(t_ol_1b)

        melt, oliv, f1, f_femg, fafm, t_ol, t_ol_1b, kd, xfo, new_comp = calculate(afm_sol.iloc[[i]], 0, 'add',
                                                                          p[i])
        f_afm.append(fafm[0])
        f1_afm.append(f1)
        x_fo_afm.append(xfo)
        kd_afm.append(kd)
        tol_afm.append(t_ol)
        tol_afm_1b.append(t_ol_1b)

    batch_solution.insert(4, 'Xfo', np.round(x_fo_batch,2))
    afm_solution.insert(4, 'Xfo', np.round(x_fo_afm,2))

    batch_solution.insert(5, 'Mg#', 1/(1 + (((1/np.array(x_fo_batch)) - 1)/np.array(kd_batch))))
    afm_solution.insert(5, 'Mg#', 1/(1 + (((1/np.array(x_fo_afm)) - 1)/np.array(kd_batch))))

    batch_solution.insert(6, 'KD', kd_batch)
    afm_solution.insert(6, 'KD', kd_afm)

    batch_solution.insert(7, 'T Ol/Liq (°C)', tol_batch)
    afm_solution.insert(7, 'T Ol/Liq (°C)', tol_afm)

    batch_solution.insert(8, 'F (Fe/Mg)', f_batch)
    afm_solution.insert(8, 'F AFM', f_afm)

    #batch_solution.insert(8, 'F1', f1_batch)
    #afm_solution.insert(8, 'F1', f1_afm)

    Pi_batch = pressure(batch_solution)
    Pi_afm = pressure(afm_solution)

    #batch_solution.insert(6, 'TP (°C)', [t_ol_batch_1b + 39.16 * pi - 14.266 * pi for t_ol_batch_1b, pi in zip(tol_batch_1b, Pi)])
    #afm_solution.insert(6, 'TP (°C)', [t_ol_afm_1b + 39.16 * pi - 14.266 * pi for t_ol_afm_1b, pi in zip(tol_afm_1b, Pi)])
    batch_solution.insert(9, 'TP (°C)', [(-8.5885 * pi ** 2 + 154.08 * pi + 1038.4) - 14.266 * pi for pi in Pi_batch])
    afm_solution.insert(9, 'TP (°C)',[(-8.5885 * pi ** 2 + 154.08 * pi + 1038.4) - 14.266 * pi for pi in  Pi_afm])

    batch_solution.insert(10, 'Pi', Pi_batch)
    afm_solution.insert(10, 'Pi', Pi_afm)
    # afm_solution.insert(8, 'Pf', Pf)
    # afm_solution.insert(9, 'T Pf (°C)', t_pf)
    # afm_solution.insert(7, 'Residue', source)

    batch_solution = batch_solution.apply(pd.to_numeric, errors='ignore')
    afm_solution = afm_solution.apply(pd.to_numeric, errors='ignore')


    return batch_solution, afm_solution, figure  # , fig_pt


def work(data):
    while True:
        try:
            return final_results(data)
        except:
            return False


