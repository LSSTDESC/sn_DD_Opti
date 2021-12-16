import pandas as pd
from __init__ import plt
import numpy as np
from optparse import OptionParser
from wrapper import DD_Budget


def plot(df, xvar='year', yvar='sigma_w', yvartwin='', legx='Year', legy='$\sigma_w$'):

    fig, ax = plt.subplots(figsize=(15, 8))
    lls = ['solid', 'dashed', 'dotted', 'dotted', 'solid', 'dashed', ]
    keys = list(df.keys())
    ls = dict(zip(keys, lls[:len(keys)]))
    for key, vals in df.items():
        ax.plot(vals[xvar], vals[yvar], marker='o', label=key, ls=ls[key])

    if yvartwin != '':
        axb = ax.twinx()
        for key, vals in df.items():
            axb.plot(vals[xvar], vals[yvartwin],
                     marker='o', label=key, ls=ls[key])

    if xvar == 'year':
        ax.set_xlim([0.9, 10.1])
    ax.grid()
    ax.legend()
    ax.set_xlabel(legx)
    ax.set_ylabel(legy)


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
parser.add_option('--cosmoFiles', type='str', default='cosmoSN_deep_rolling_2_2_mini_yearly_Ny_40,cosmoSN_deep_rolling_2_2_mini_0.65_yearly_Ny_40,cosmoSN_deep_rolling_2_2_mini_0.60_yearly_Ny_40,cosmoSN_deep_rolling_0.80_0.80_yearly_Ny_40,cosmoSN_universal_yearly_Ny_40',
                  help='cosmo files to process[%default]')
parser.add_option('--nickNames', type='str', default='deep_rolling_early,deep_rolling_early_0.65,deep_rolling_early_0.60,deep_rolling_ten_years,universal',
                  help='nicknames corresponding to cosmofiles to process[%default]')
parser.add_option('--prefix_Nvisits', type='str', default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help='prefix for Nvisits file[%default]')

opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
visitsDir = opts.visitsDir
cadence = opts.cadence
cosmoFiles = opts.cosmoFiles.split(',')
nickNames = opts.nickNames.split(',')
prefix_Nvisits = opts.prefix_Nvisits

df = {}

for i in range(len(cosmoFiles)):
    fi = '{}/{}.hdf5'.format(cosmoDir, cosmoFiles[i])
    nName = nickNames[i]
    if 'universal' in nName:
        # df['universal_0.60'] = load_b(fi, zcomp=0.60)
        df['universal_0.65'] = load_b(fi, zcomp=0.65)
    else:
        df[nName] = load(fi)

# estimate Budget for all of these
# budget class instance
for key, vals in df.items():
    Ny = '{}'.format(int(np.median(vals['Ny'])))
    fName = '{}_{}.npy'.format(prefix_Nvisits, Ny)
    budget = DD_Budget(visitsDir, fName)
    # add budget colmuns
    cosmo_bud = budget(vals)
    df[key] = cosmo_bud

# plot_multiple(df, yvartwin='budget')
plot(df)
plot(df, yvar='budget')
# plot(df, xvar='nsn_z_09')
# plot(df, yvar='budget')
plt.subplots_adjust(hspace=0, wspace=0)
plt.show()
