import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

plt.rcParams.update({'font.family': 'Arial'})
plt.rcParams.update({'font.size': 12})


def batch_melt(mgo_s, feo_s, f):
    mgo_l = np.arange(2, 40.1, 0.1)
    D_mgo = ((mgo_s / mgo_l) - f) / (1 - f)
    kd = 0.381 - 0.774 / mgo_l + 0.998 / mgo_l ** 2
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
        kd_ol = 0.381 - 0.790 / mgo_l + 1.039 / mgo_l ** 2

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
    mgo_s = 38.12
    feo_s = 8.02

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
    ax1.plot(mgo_s, feo_s, 'o', c='g', zorder=6, markersize=20, label='KR-4003')
    ax1.set(ylim=[4, 12], xlim=[2, 40], ylabel='FeO (wt%)', xlabel='MgO (wt%)', title='BATCH MELTING')
    ax1.legend(loc="lower right")


def make_figure_afm(ax2):
    mgo_s = 38.12
    feo_s = 8.02

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
    ax2.plot(mgo_s, feo_s, 'o', c='g', zorder=6, markersize=20, label='KR-4003')
    ax2.set(ylim=[4, 12], xlim=[2, 40], ylabel='FeO (wt%)', xlabel='MgO (wt%)', title='AFM MELTING')
    ax2.legend(loc="lower right")
    # ax2.text(5.5, 5, 'SOLIDUS', rotation=65)
    # ax2.text(39, 5, 'Accumulated Fractional Melting', ha='right', size=14, zorder=7)
    # ax2.text(18.6, 11.1, 0.1, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(20.5, 11.5, 0.2, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(22.5, 11.5, 0.3, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(24.2, 11.5, 0.4, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(27.5, 11.5, 0.5, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(31.25, 11.25, 0.6, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(37, 10.9, 0.7, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(38, 9.9, 0.8, bbox=dict(ec='k', fc='w'), size='xx-small')
    # ax2.text(39, 11.7, 'Melt Fraction (F)', size='x-small', ha='right')
    # ax2.text(39, 7, 'Kettle River peridotite\nKR-4003', size='small', ha='right', zorder=8)


def figure_pt(results_afm):
    mgo_s = 38.12
    feo_s = 8.02
    plt.rcParams.update({'font.family': 'Arial'})
    plt.rcParams.update({'font.size': 12})

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(111)

    Pf = results_afm.Pf.values
    T_Pf = results_afm['T Pf (°C)'].values
    names = list(results_afm.index.values)

    p_solidus = np.linspace(0, 10, 101)
    p_solidus[0] = 0.01
    t_solidus = 1086 - 5.7 * p_solidus + 390 * np.log(p_solidus)
    t_solidus[0:27] = 132 * p_solidus[0:27] + 1102
    p_liquidus = np.arange(0, 10.1, 0.1)
    t_liquidus = 1020 + 24.4 * mgo_s - 0.161 * mgo_s ** 2 + 54 * p_liquidus - 2 * p_liquidus ** 2

    normalize = matplotlib.colors.Normalize(vmin=1200, vmax=1600)
    ax1.scatter(Pf, T_Pf, c=T_Pf, cmap='turbo', norm=normalize)
    ax1.plot(p_solidus, t_solidus, c='k', label="Solidus")
    ax1.fill_between(p_solidus, t_solidus, color='gray', alpha=0.1)
    ax1.fill_between(p_liquidus, t_liquidus, 10000, color='orange', alpha=0.1)
    ax1.set(xlim=[0, 7], ylim=[1000, 2000], xlabel='P (GPa)', ylabel='T (°C)')
    ax1.plot(np.arange(0, 10.1, 0.1), t_liquidus, c='k', label="Liquidus")
    # ax1.grid(color='k', lw=0.2, alpha=0.5)
    # ax1.annotate('All solid', xy=(6.7, 1020), xycoords='data', clip_on=True, ha='right', size=14)
    # ax1.annotate('All liquid', xy=(0.5, 1900), xycoords='data', clip_on=True, ha='left', size=14)
    # ax1.annotate('Solidus', xy=(0.5, 1100), xytext=(0, 0), xycoords='data', textcoords='offset points', clip_on=True, ha='center', rotation=27, size=12)
    # ax1.annotate('Liquidus', xy=(0.5, 1750), xytext=(0, 0), xycoords='data', textcoords='offset points', clip_on=True, ha='center', rotation=12, size=12)

    for i, txt in enumerate(names):
        ax1.annotate(txt, xy=(Pf[i], T_Pf[i]), xytext=(0, 8), xycoords='data', textcoords='offset points', clip_on=True,
                     ha='center', va='center', size=8)

    ax1.legend(loc='lower right')
    return fig


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
    z = x * data['FeO'] / 0.9
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
    mw = np.asarray(
        [60.085, 79.899, 101.962, 151.99, 159.962, 71.846, 70.937, 40.311, 56.079, 61.979, 94.203, 74.71, 141.945,
         18.0152])
    data_mol = data / mw
    y = data_mol.sum(axis=1)
    data_mol = 100 * data_mol / y.values[0]

    cp = np.asarray(
        [60.085, 79.899, 101.962 / 2, 151.99 / 2, 159.962 / 2, 71.846, 70.937, 40.311, 56.079, 61.979 / 2, 94.203 / 2,
         74.71, 141.945 / 2, 18.0152 / 2])
    data_cat = data / cp
    data_cat = data_cat.set_axis(['Si', 'Ti', 'Al', 'Cr', 'Fe3', 'Fe', 'Mn', 'Mg', 'Ca', 'Na', 'K', 'Ni', 'P', 'H'],
                                 axis=1, inplace=False)
    z = data_cat.sum(axis=1)
    data_cat = data_cat / z.values[0]

    return data_mol, data_cat


def temp_olivine(cat, p):
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
    t_ol = t_ol + 54 * (p / 10000) - 2 * (p / 10000) ** 2

    print('olivine temperature OK!')
    return t_ol[0]


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
    R = 8.3143
    T = t_ol + 273.15
    sio2n = sio2n
    Fe = data.FeO[0] / 71.846
    Mg = data.MgO[0] / 40.311
    xfo = 0.9
    kd_old = 0.0
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
        [60.085, 79.899, 101.962 / 2, 151.99 / 2, 159.962 / 2, 71.846, 70.937, 40.311, 56.079, 61.979 / 2, 94.203 / 2,
         74.71, 141.945 / 2, 18.0152 / 2])

    z = ox_ol * cp
    ox_ol = 100 * ox_ol * cp / z.sum()

    ol_comp = pd.DataFrame(
        columns=['SiO2', 'TiO2', 'Al2O3', 'Cr2O3', 'Fe2O3', 'FeO', 'MnO', 'MgO', 'CaO', 'Na2O', 'K2O', 'NiO', 'P2O5',
                 'H2O'])
    ox_ol = pd.Series(ox_ol, index=ol_comp.columns)

    ol_comp = ol_comp.append(ox_ol, ignore_index=True)

    print('equilibrium olivine OK!')
    return ol_comp


def kd_f_values(comp, mgo_s, feo_s):
    print(f"kd_f_values called with mgo_s={mgo_s}, feo_s={feo_s}")
    if mgo_s is None or feo_s is None:
        raise ValueError("mgo_s or feo_s is None")
    mgo = comp.MgO[0]
    feo = comp.FeO[0]
    kd_new = 0.3813 - 0.7896 / mgo + 1.0389 / (mgo ** 2)
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


def f_values(comp_mol):
    coef = np.array([[0.0, -1.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, -0.5, -0.5, 3, 0.5, 0.0, 0.0],
                     [0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     [1.0, 0.0, -0.5, -0.5, 0.0, -0.5, -0.5, -0.5, -1.5, -3.0, -3.0, -0.5, 0.0, 0.0]])

    comp_mol = comp_mol[0].values
    ol_an_q = coef @ comp_mol.T
    suma = ol_an_q.sum()
    ol_an_q = ol_an_q / suma
    ol_an_q = ol_an_q.astype(float)

    s_cpx = -0.074 + 0.1713 / ol_an_q[0] - 0.0135 / ol_an_q[0] ** 2
    s_gtlhzhz = 0
    if ol_an_q[1] < 700:
        s_gtlhzhz = 1 / (16.843 + 28.733 * ol_an_q[1] - 14.183 * np.exp(ol_an_q[1]))

    dom = 2
    if ol_an_q[0] > 0.5 and ol_an_q[2] < s_gtlhzhz:
        dom = 3
    elif ol_an_q[2] > s_cpx:
        dom = 1

    f1 = 6.2819 * ol_an_q[1] ** 2 - 14.7789 * ol_an_q[1] ** 3 + 0.00825 * (1 / ol_an_q[1]) ** 2
    f2_temp = 0
    f3_temp = -2.5345 + 5.329 * (ol_an_q[2] + ol_an_q[0] * 0.348) + 0.3012 / (ol_an_q[2] + ol_an_q[0] * 0.348)
    if ol_an_q[2] > 0:
        f2_temp = ((-1.994 + 2.25 * ol_an_q[2] + 0.041 / ol_an_q[2]) + (
                -1.183 - 3.005 * ol_an_q[2] + 13.774 * ol_an_q[2] ** 2 - 12.615 * ol_an_q[2] ** 3)) / 2 + np.exp(
            0.931 + 1.623 * ol_an_q[2]) * ol_an_q[2] ** 0.245 * ol_an_q[0] + np.exp(0.769 - 7.514 * ol_an_q[2]) * \
                  ol_an_q[2] ** 0.577 / ol_an_q[0]

    f2 = f2_temp
    if f2_temp > 0 and ol_an_q[2] < (-0.1773 + 0.1404 / ol_an_q[0] - 0.008434 / ol_an_q[0] ** 2):
        f2 = -f2_temp

    f3 = f3_temp
    if f3_temp > 0 and ol_an_q[2] < (-0.0896 + 0.02002 / ol_an_q[0] + 0.02989 / ol_an_q[0] ** 2):
        f3 = -f3_temp

    return dom, f1, f2, f3, s_cpx, s_gtlhzhz


def calculate(comp, ol_frac, ope, p):
    ol_frac = 0.01
    comp_mol = mol_cat_oxides(comp)
    t_ol = temp_olivine(comp_mol[1], p)
    sio2n = sio2pound(comp, comp_mol)
    kd, xfo = calc_kd(comp, t_ol, sio2n, p)
    comp_ol = eq_olivine(comp_mol, xfo)
    kd_new, f_femg, fafm = kd_f_values(comp, mgo_s = config_manager.get_setting('mgo_s'), feo_s = config_manager.get_setting('feo_s'))
    dom, f1, f2, f3, s_cpx, s_gtlhzhz = f_values(comp_mol)

    if ope == 'add':
        new_comp = (comp + ol_frac * comp_ol) / (1 + ol_frac)
    if ope == 'sub':
        new_comp = (comp - ol_frac * comp_ol) / (1 - ol_frac)

    return comp, comp_ol, dom, f1, f2, f3, f_femg, fafm, t_ol, kd, xfo, new_comp


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

    dom_list = []
    t_olivine = []
    kd_value = []
    xfo_value = []
    ol_add = []

    x = 0
    sol = False

    if ope == 'add':
        while j <= 0.60 + ol_frac:
            melt, oliv, dom, f1, f2, f3, f_femg, fafm, t_ol, kd, xfo, new_comp = calculate(comp, ol_frac, ope, p)
            _melt, _oliv, _dom, f1_2, f2_2, f3_2, f_femg_2, fafm_2, _t_ol, _kd, _xfo, _new_comp = calculate(new_comp,
                                                                                                            ol_frac,
                                                                                                            ope, p)
            comp = new_comp
            magma = magma.append(melt.iloc[0])
            olivine = olivine.append(oliv.iloc[0])
            ol_add.append(j)

            if typ == 'afm':
                f_femg = fafm
                f_femg_2 = fafm_2

            if abs(f_femg) >= 1:
                break

            check_f1 = find_intersection(f1, f1_2, f_femg, f_femg_2)
            check_f2 = find_intersection(f2, f2_2, f_femg, f_femg_2)
            check_f3 = find_intersection(f3, f3_2, f_femg, f_femg_2)

            dom_list.append(dom)
            t_olivine.append(t_ol)
            kd_value.append(kd)
            xfo_value.append(xfo)

            if (1 in dom_list or dom == 1) and (check_f1 != False) and (f_femg >= 0):
                # if not (dom == 2 and check_f2 != False) or (dom == 3 and check_f3 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f1

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            if (2 in dom_list or dom == 2) and (check_f2 != False) and f_femg >= 0:
                # if not (dom == 1 and check_f1 != False) or (dom == 3 and check_f3 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f2

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            if (3 in dom_list or dom == 3) and (check_f3 != False) and f_femg >= 0:
                # if not (dom == 1 and check_f1 != False) or (dom == 2 and check_f2 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f3

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            j = j + ol_frac
            i += 1

    if ope == 'sub':
        while j <= 0.40 + ol_frac:
            melt, oliv, dom, f1, f2, f3, f_femg, fafm, t_ol, kd, xfo, new_comp = calculate(comp, ol_frac, ope, p)
            _melt, _oliv, _dom, f1_2, f2_2, f3_2, f_femg_2, fafm_2, _t_ol, _kd, _xfo, _new_comp = calculate(new_comp,
                                                                                                            ol_frac,
                                                                                                            ope, p)
            comp = new_comp
            magma = magma.append(melt.iloc[0])
            olivine = olivine.append(oliv.iloc[0])
            ol_add.append(j)

            if typ == 'afm':
                f_femg = fafm
                f_femg_2 = fafm_2

            if abs(f_femg) >= 1:
                break

            check_f1 = find_intersection(f1, f1_2, f_femg, f_femg_2)
            check_f2 = find_intersection(f2, f2_2, f_femg, f_femg_2)
            check_f3 = find_intersection(f3, f3_2, f_femg, f_femg_2)

            dom_list.append(dom)
            t_olivine.append(t_ol)
            kd_value.append(kd)
            xfo_value.append(xfo)

            if (1 in dom_list or dom == 1) and (check_f1 != False) and (f_femg >= 0):
                # if not (dom == 2 and check_f2 != False) or (dom == 3 and check_f3 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f1

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            if (2 in dom_list or dom == 2) and (check_f2 != False) and f_femg >= 0:
                # if not (dom == 1 and check_f1 != False) or (dom == 3 and check_f3 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f2

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            if (3 in dom_list or dom == 3) and (check_f3 != False) and f_femg >= 0:
                # if not (dom == 1 and check_f1 != False) or (dom == 2 and check_f2 != False):
                magma = magma.append(_melt.iloc[0])
                olivine = olivine.append(_oliv.iloc[0])
                ol_add.append(j + ol_frac)

                x = check_f3

                solution = magma.iloc[len(magma) - 1] * (x) + magma.iloc[len(magma) - 2] * (1 - x)
                magma.iloc[len(magma) - 1] = solution

                ol = ol_add[len(magma) - 1] * (x) + ol_add[len(magma) - 2] * (1 - x)
                ol_add[len(magma) - 1] = ol[0]

                sol = True
                break

            j = j + ol_frac
            i += 1

    return magma, olivine, sol, ol_add  # , check_f1, check_f2, check_f3, dom


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
        comp = comp.append(sample)
        comp = comp.reset_index(drop=True)

        name.append(comp.Name[0])

        p = comp.P[0]
        pres.append(p)

        px_error = ''
        vt_error = ''

        if comp.CaO[0] < 13.81 - 0.274 * comp.MgO[0]:
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
            melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])
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
                melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])
                ax_batch.plot(magma_batch.MgO, magma_batch.FeO, c='k', lw=0.5)
                ax_batch.plot(magma_batch['MgO'].iloc[-1], magma_batch['FeO'].iloc[-1], 'o', mfc='w', mec='k')
                ax_batch.plot(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0], 'o', mfc='k', mec='None')
                ax_batch.annotate(text=name[i], xy=(magma_batch['MgO'].iloc[0], magma_batch['FeO'].iloc[0]),
                                  xycoords='data', xytext=(8, 4), textcoords='offset points', ha='center', va='bottom',
                                  size='x-small')
            else:
                sol_batch.append('N')
                op_batch.append(0)
                melt_sol_batch = melt_sol_batch.append(magma_batch.iloc[len(magma_batch) - 1])

        ol_add_batch.append(100 * ol_add[len(ol_add) - 1])

        # AFM
        magma_afm, olivine_afm, sol, ol_add = fractionate_olivine(comp, ol_frac, 'add', 'afm', p)

        if sol == True:
            sol_afm.append('Y')
            op_afm.append(1)
            melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])
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
                melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])
                ax_afm.plot(magma_afm.MgO, magma_afm.FeO, c='k', lw=0.5)
                ax_afm.plot(magma_afm['MgO'].iloc[-1], magma_afm['FeO'].iloc[-1], 'o', mfc='w', mec='k')
                ax_afm.plot(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0], 'o', mfc='k', mec='None')
                ax_afm.annotate(text=name[i], xy=(magma_afm['MgO'].iloc[0], magma_afm['FeO'].iloc[0]), xycoords='data',
                                xytext=(8, 4), textcoords='offset points', ha='center', va='bottom', size='x-small')
            else:
                sol_afm.append('N')
                op_afm.append(0)
                melt_sol_afm = melt_sol_afm.append(magma_afm.iloc[len(magma_afm) - 1])

        ol_add_afm.append(100 * ol_add[len(ol_add) - 1])

    melt_sol_batch.index = name
    melt_sol_batch.index.name = 'Name'
    melt_sol_batch.insert(0, '%Ol', [a * b for a, b in zip(ol_add_batch, op_batch)])
    melt_sol_batch.insert(0, 'Sol', sol_batch)
    melt_sol_batch.insert(0, 'P', pres)
    melt_sol_batch.insert(0, 'Input Error', input_error)

    melt_sol_afm.index = name
    melt_sol_afm.index.name = 'Name'
    melt_sol_afm.insert(0, '%Ol', [a * b for a, b in zip(ol_add_afm, op_afm)])
    melt_sol_afm.insert(0, 'Sol', sol_afm)
    melt_sol_afm.insert(0, 'P', pres)
    melt_sol_afm.insert(0, 'Input Error', input_error)

    return melt_sol_batch, melt_sol_afm, pres, input_error, fig


def pressure(afm_solution):
    Pi_ = []
    Pf_ = []
    t_pf = []
    dp = []

    a1, a2, a3, b1, b2, b3 = -0.13, 0.0101, 120, 1.74, -46.2, 555
    x1, x2, x3, y1, y2, y3, z1, z2, z3 = -21.37, -72.7, 1168, 9.079, 6097.2, -18506.4, -0.0144, -5956, 52903

    for i in range(len(afm_solution)):
        mgo = afm_solution.MgO[i]
        fafm = afm_solution['F AFM'].iloc[i]

        pi = 2.5
        if mgo >= 12.3:
            pi = 11.248 * mgo - 13700 / (mgo ** 3) - 8.13 * (np.log(mgo) ** 3)

        A = a1 + a2 * pi ** 3 + a3 / pi ** 3
        B = b1 * pi + b2 / pi + b3 / pi ** 3

        X = x1 + x2 / pi + x3 / pi ** 3
        Y = y1 * pi + y2 / pi ** 2 + y3 / pi ** 3
        Z = z1 * pi ** 4 + z2 / pi ** 2 + z3 / pi ** 4

        dPPD = A * fafm + B * fafm ** 2
        dPHZ = X + Y * fafm + Z * fafm ** 2
        dPPDHZ = (0.0831 * pi ** 2) + (20.15 / pi ** 2) - (1.19 * np.log(pi))

        dP = dPPD
        if dPPD > dPPDHZ:
            dP = dPHZ

        Pi = -0.481 * mgo + 0.000454 * mgo ** 3 + 1.166 * np.log(mgo) ** 2
        Pf = Pi - dP

        t_p = 1020 + 24.4 * mgo - 0.161 * mgo ** 2 + 54 * Pf - 2 * Pf ** 2

        if afm_solution.Sol[i] == 'N':
            Pi = Pf = t_p = np.nan

        Pi_.append(Pi)
        Pf_.append(Pf)
        t_pf.append(t_p)
        dp.append(dP)

    return Pi_, Pf_, t_pf


def final_results(data):
    batch_solution, afm_solution, p, input_error, figure = run_all(data)
    batch_sol = batch_solution.drop(['Sol', '%Ol', 'Input Error', 'P'], axis=1)
    afm_sol = afm_solution.drop(['Sol', '%Ol', 'Input Error', 'P'], axis=1)

    f_batch = []
    f_afm = []
    x_fo_batch = []
    x_fo_afm = []
    kd_batch = []
    kd_afm = []
    tol_batch = []
    tol_afm = []
    source = []

    for i in range(len(batch_sol)):
        melt, oliv, dom, f1, f2, f3, f_femg, fafm, t_ol, kd, xfo, new_comp = calculate(batch_sol.iloc[[i]], 0, 'add',
                                                                                       p[i])
        f_batch.append(f_femg)
        x_fo_batch.append(xfo)
        kd_batch.append(kd)
        tol_batch.append(t_ol)

        melt, oliv, dom, f1, f2, f3, f_femg, fafm, t_ol, kd, xfo, new_comp = calculate(afm_sol.iloc[[i]], 0, 'add',
                                                                                       p[i])
        f_afm.append(fafm)
        x_fo_afm.append(xfo)
        kd_afm.append(kd)
        tol_afm.append(t_ol)

        # DETECTS SOURCE
        if dom == 1:
            src = 'Harz'
        if dom == 2:
            src = 'Sp-Per'
        if dom == 3:
            src = 'Grt-Per'
        source.append(src)

    batch_solution.insert(2, 'Xfo', x_fo_batch)
    afm_solution.insert(2, 'Xfo', x_fo_afm)

    batch_solution.insert(3, 'KD', kd_batch)
    afm_solution.insert(3, 'KD', kd_afm)

    batch_solution.insert(4, 'T Ol/Liq (°C)', tol_batch)
    afm_solution.insert(4, 'T Ol/Liq (°C)', tol_afm)

    batch_solution.insert(5, 'F (Fe/Mg)', f_batch)
    afm_solution.insert(5, 'F AFM', f_afm)

    batch_solution.insert(6, 'TP (°C)', 1025 + 28.6 * batch_solution.MgO - 0.084 * batch_solution.MgO ** 2)
    afm_solution.insert(6, 'TP (°C)', 1025 + 28.6 * afm_solution.MgO - 0.084 * afm_solution.MgO ** 2)

    Pi, Pf, t_pf = pressure(afm_solution)
    afm_solution.insert(7, 'Pi', Pi)
    afm_solution.insert(8, 'Pf', Pf)
    afm_solution.insert(9, 'T Pf (°C)', t_pf)
    afm_solution.insert(10, 'Residue', source)

    sol_afm = afm_solution.pop('Sol')
    afm_solution.insert(0, 'Solution', sol_afm)

    sol_batch = batch_solution.pop('Sol')
    batch_solution.insert(0, 'Solution', sol_batch)
    batch_solution = batch_solution.drop('H2O', axis=1)

    batch_solution = batch_solution.apply(pd.to_numeric, errors='ignore')
    afm_solution = afm_solution.apply(pd.to_numeric, errors='ignore')
    afm_solution = afm_solution.drop('H2O', axis=1)

    fig_pt = figure_pt(afm_solution)

    return batch_solution, afm_solution, figure, fig_pt

def work(data):
    while True:
        try:
            return final_results(data)
        except:
            return False

