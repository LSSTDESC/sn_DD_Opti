import pandas as pd
from __init__ import plt
import numpy as np
from optparse import OptionParser
from wrapper import DD_Budget
from scipy.interpolate import interp1d


def plot(df, xvar='year', yvar='sigma_w', yvartwin='', legx='Year', legy='$\sigma_w$', tag_budget=[], tag_marker=[]):

    vva = ['deep_rolling_early', 'deep_rolling_ten_years', 'universal']
    vvp = ['Early Deep Rolling (EDR$^3$)',
           'Deep Rolling 10 years (DR$^{10Y}$)', 'Deep Universal (DU$^{0.65}$)']
    colors = ['red', 'blue', 'magenta']
    corresp = dict(zip(vva, vvp))
    ccolors = dict(zip(vva, colors))
    fig, ax = plt.subplots(figsize=(15, 8))
    lls = ['solid', 'dashed', 'dotted', 'dotted', 'solid', 'dashed', ]
    keys = list(df.keys())
    ls = dict(zip(keys, lls[:len(keys)]))
    for key, vals in df.items():
        print('elelelele', vals.columns)
        ax.plot(vals[xvar], vals[yvar], marker='o',
                label=corresp[key], ls=ls[key], ms=5, color=ccolors[key])
        if len(tag_budget) > 0:
            interp_bud = interp1d(
                vals['budget'], vals[xvar], bounds_error=False, fill_value=0.)
            interp_var = interp1d(
                vals[xvar], vals[yvar], bounds_error=False, fill_value=0.)
        print('Survey', key)
        for i, val in enumerate(tag_budget):
            ttime = interp_bud(val)
            valres = interp_var(ttime)
            ax.plot(ttime, valres,
                    marker=tag_marker[i], color='k', ms=15., mfc='None', markeredgewidth=2)
            print('resu budget', val, xvar, ttime, yvar, valres)

    if yvartwin != '':
        axb = ax.twinx()
        for key, vals in df.items():
            axb.plot(vals[xvar], vals[yvartwin],
                     marker='o', label=key, ls=ls[key])

    if xvar == 'year':
        ax.set_xlim([0.9, 10.1])
    if yvar == 'sigma_w':
        ax.set_ylim([0.005, None])
    ax.grid()
    ax.legend()
    ax.set_xlabel(legx)
    ax.set_ylabel(legy)


def plot_syste(df, df_syste, xvar='year', yvar='sigma_w', yvartwin='', legx='Year', legy='$\sigma_w$', tag_budget=[], tag_marker=[], norm=False):

    vva = ['deep_rolling_early', 'deep_rolling_ten_years', 'universal']
    vvp = ['Early Deep Rolling',
           'Deep Rolling 10 years', 'Deep Universal']
    vvp = ['Early Deep Rolling (EDR$^3$)',
           'Deep Rolling 10 years (DR$^{10Y}$)', 'Deep Universal (DU$^{0.65}$)']
    colors = ['red', 'blue', 'magenta']
    corresp = dict(zip(vva, vvp))
    ccolors = dict(zip(vva, colors))
    fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(15, 8))
    lls = ['solid', 'dashed', 'dotted', 'dotted', 'solid', 'dashed', ]
    keys = list(df.keys())
    ls = dict(zip(keys, lls[:len(keys)]))
    for key, vals in df.items():
        # ax.plot(vals[xvar], vals[yvar], marker='o',
        #        label=corresp[key], ls=ls[key], ms=5, color=ccolors[key])
        # get systematic (if exist)
        key_syste = []
        for keyb in df_syste.keys():
            if key in keyb:
                key_syste.append(keyb)

        for io, keysyst in enumerate(key_syste):

            dfsyst = df_syste[keysyst]
            dfsyst['year'] = dfsyst['year'].astype(int)
            vals['year'] = vals['year'].astype(int)

            dfb = vals.merge(dfsyst, left_on=[xvar], right_on=[xvar])
            varx = '{}_x'.format(yvar)
            vary = '{}_y'.format(yvar)
            #print('hello', dfb[[varx, vary]])
            dfb['diff'] = dfb[varx]-dfb[vary]
            if norm:
                dfb['diff'] /= dfb['{}_x'.format(yvar)]
            ax[io].plot(dfb[xvar], dfb['diff'], marker='o',
                        ls=ls[key], ms=5, color=ccolors[key], label=corresp[key])

            if len(tag_budget) > 0:
                interp_bud = interp1d(
                    vals['budget'], vals[xvar], bounds_error=False, fill_value=0.)
                interp_varm = interp1d(
                    vals[xvar], vals[yvar], bounds_error=False, fill_value=0.)
                interp_var = interp1d(
                    dfb[xvar], dfb['diff'], bounds_error=False, fill_value=0.)

            print('survey', key, keysyst)
            for i, val in enumerate(tag_budget):
                ttime = interp_bud(val)
                valres = interp_var(ttime)
                valresm = interp_varm(ttime)
                ax[io].plot(ttime, valres,
                            marker=tag_marker[i], color='k', ms=15., mfc='None', markeredgewidth=2)
                print('resu budget', val, xvar,
                      ttime, yvar, valres, valresm)
            if xvar == 'year':
                ax[io].set_xlim([0.9, 10.1])
                # if yvar == 'sigma_w':
                #    ax.set_ylim([0.005, None])
            ax[io].grid()
            if io == 0:
                ax[io].legend(fontsize=20)
                # ax[io].legend(loc='upper left', bbox_to_anchor=(
                #    -0.1, 1.3), ncol=3, frameon=False, fontsize=20)
            ax[io].set_xlabel(legx)
            ax[io].set_ylabel(legy)


def plot_multiple(df, xvar='year', yvar='sigma_w', yvartwin=''):

    nplots = len(df.keys())
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(
        15, 8), sharex=True, sharey=True)
    ccorr = dict(zip(range(nplots), [(0, 0), (0, 1), (1, 0), (1, 1)]))
    ppos = [(0, 0), (0, 1), (1, 0), (1, 1)]
    nn = ['deep_rolling_early', 'universal_0.65',
          'deep_rolling_ten_years', 'universal_0.60']
    ccorr = dict(zip(nn, [(0, 0), (0, 1), (1, 0), (1, 1)]))
    keys = list(df.keys())
    for key, vals in df.items():
        idx = vals['budget'] <= 0.12
        vals = vals[idx]
        ipos = ccorr[key][0]
        jpos = ccorr[key][1]
        ax[ipos, jpos].plot(vals[xvar], vals[yvar],
                            marker='o', color='r')
        ax[ipos, jpos].text(2, 0.05, key)
        if yvartwin != '':
            axb = ax[ipos, jpos].twinx()
            axb.plot(vals[xvar], vals[yvartwin],
                     marker='o', color='k', ls='dashed')

        ax[ipos, jpos].grid()
        # ax[ipos, jpos].legend()
        xlims = ax[ipos, jpos].get_xlim()
        axb.plot(xlims, [0.05]*len(xlims), ls='dotted', color='b')
        axb.plot(xlims, [0.1]*len(xlims), ls='dotted', color='b')
        # axb.plot(xlims, [0.15]*len(xlims), ls='dotted', color='b')
        ax[ipos, jpos].set_xlim([0.95, 10])
        if 'deep' in key:
            axb.set_yticklabels('')
            ax[ipos, jpos].set_ylabel('$\sigma_w$')
        else:
            axb.set_ylabel('DD budget')
        if ipos == 1:
            ax[ipos, jpos].set_xlabel('Year')

    fig.tight_layout()


def plot_diff(df, ref=['deep_rolling_early', 'deep_rolling_ten_years', 'universal'], xvar='year', yvar='sigma_w', yvartwin='', legx='Year', legy='$\sigma_w$'):

    fig, ax = plt.subplots(figsize=(15, 8))

    for refval in ref:
        df_ref = df[refval]
        for key, vals in df.items():
            if refval in key and refval != key:
                print('boo', refval, key)
                ax.plot(df_ref[xvar], df_ref[yvar]-vals[yvar], label=key)

    ax.legend()
    ax.grid()


def load(fi):
    """
    Function to load file in pandas df

    Parameters
    ---------------
    fi: str
      name of the file (full path) to be loaded

    Returns
    ----------
    pandas df corresponding to the file
    """
    df = pd.read_hdf(fi)

    df['year'] = df['conf'].str.split('_').str.get(1).astype(int)
    df = df.sort_values(by=['year'])

    print('ici', df)
    return df


def load_b(fi, zcomp=0.65):
    """
    Function to load file in pandas df

    Parameters
    ---------------
    fi: str
      name of the file (full path) to be loaded
    zcomp: float, opt
      redshift completeness (default: 0.65)

    Returns
    ----------
    pandas df corresponding to the file
    """
    df = pd.read_hdf(fi)

    df['zcomp_dd_unique'] = df.apply(lambda x: np.mean(x['zcomp_dd']), axis=1)
    idx = np.abs(df['zcomp_dd_unique']-zcomp) < 1.e-5

    sel = df[idx]
    sel = sel.sort_values(by=['sigma_w'], ascending=False)
    sel['year'] = sel.reset_index().index+1

    print(sel)

    return sel


parser = OptionParser()
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option("--cosmoDir", type=str, default='cosmo_files',
                  help="directory where cosmo files are located[%default]")
parser.add_option('--cadence', type=float, default=1.,
                  help='cadence of observation[%default]')
parser.add_option('--config', type='str', default='config_survey_time.csv',
                  help='configuration file [%default]')
parser.add_option('--prefix_Nvisits', type='str', default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help='prefix for Nvisits file[%default]')

opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
visitsDir = opts.visitsDir
cadence = opts.cadence
prefix_Nvisits = opts.prefix_Nvisits

conf = pd.read_csv(opts.config, comment='#')

df = {}
df_syste = {}

"""
for i in range(len(cosmoFiles)):
    fi = '{}/{}.hdf5'.format(cosmoDir, cosmoFiles[i])
    nName = nickNames[i]
    if 'universal' in nName:
        # df['universal_0.60'] = load_b(fi, zcomp=0.60)
        df['universal_0.65'] = load_b(fi, zcomp=0.65)
    else:
        df[nName] = load(fi)
"""
for index, row in conf.iterrows():
    fi = '{}/{}.hdf5'.format(cosmoDir, row['cosmoFile'])
    nName = row['name']
    ttype = row['type']
    if 'universal' in nName:
        # df['universal_0.60'] = load_b(fi, zcomp=0.60)
        fload = load_b(fi, zcomp=0.65)
    else:
        fload = load(fi)
    fload['sigma_w_rel'] = np.abs(fload['sigma_w']/fload['w'])
    if ttype == 'nominal':
        df[nName] = fload
    else:
        df_syste[nName] = fload

# estimate Budget for all of these
# budget class instance
for key, vals in df.items():
    Ny = '{}'.format(int(np.median(vals['Ny'])))
    fName = '{}_{}.npy'.format(prefix_Nvisits, Ny)
    budget = DD_Budget(visitsDir, fName)
    # add budget colmuns
    cosmo_bud = budget(vals)
    df[key] = cosmo_bud
    df[key]['frac_high_z'] = df[key]['nsn_z_09']/df[key]['nsn_DD']
    df[key]['nsn_DD_ultra'] = df[key]['nsn_DD_COSMOS']+df[key]['nsn_DD_XMM-LSS']

# plot_multiple(df, yvartwin='budget')
# plot(df, yvar='w', legy='w')
# plot(df, yvar='Om', legy='Om')

"""
plot(df, yvar='sigma_w', legy='$\sigma_w$', tag_budget=[
     0.05, 0.0788, 0.10], tag_marker=['*', '^', 's'])
plot(df, xvar='nsn_DD', legx='$N_{SN}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])

plot(df, xvar='nsn_DD_ultra', legx='$N_{SN}^{ultra}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])
"""
plot_syste(df, df_syste, yvar='sigma_w', tag_budget=[
    0.05, 0.0788, 0.10], legy='$\Delta\sigma_w$', tag_marker=['*', '^', 's'], norm=True)


"""
plot_syste(df, df_syste, yvar='nsn_DD', tag_budget=[
    0.05, 0.08, 0.10], legy='$\Delta N_{SN}$', tag_marker=['*', '^', 's'], norm=True)
"""
"""
plot_syste(df, df_syste, yvar='w',
           legy='$w$', tag_marker=['*', '^', 's'])

plot(df, xvar='nsn_z_09', legx='$N_{SN}^{z\geq 0.9}$', tag_budget=[
     0.05, 0.08, 0.10], tag_marker=['*', '^', 's'])
plot(df, xvar='nsn_DD', legx='$N_{SN}$', tag_budget=[
     0.05, 0.08, 0.10], tag_marker=['*', '^', 's'])
plot(df, yvar='nsn_ultra', legy='$N_{SN}^{ultra}$', tag_budget=[
     0.05, 0.08, 0.10], tag_marker=['*', '^', 's'])
"""

# plot(df, yvar='sigma_Om', legy='$\sigma_{\Omega_m}$', tag_budget=[
#     0.05, 0.07, 0.09], tag_marker=['*', '^', 'v'])
# plot(df, yvar='budget', legy='budget')
# plot_diff(df)
# plot_diff(df, yvar='w')
# plot(df, yvar='budget', legy='budget')
# plot(df, xvar='nsn_z_09')
# plot(df, yvar='budget')
plt.subplots_adjust(hspace=0, wspace=0)
plt.show()