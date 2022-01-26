import pandas as pd
from __init__ import plt
import numpy as np
from optparse import OptionParser
from wrapper import DD_Budget
from scipy.interpolate import interp1d


class Budget_Time:
    """
    class to estimate for each configuration some infos corresponding to a set of budgets

    Parameters
    ---------------
    config: str
      list of files to process (csv file)
    budget: list(float), opt
      list of budget to consider (default: [0.05,0.08,0.1])


    """

    def __init__(self, config, cosmoDir='cosmo_files_yearly',
                 prefix_Nvisits='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                 visitsDir='visits_files', budget=[0.05, 0.08, 0.1]):

        conf = pd.read_csv(config, comment='#')
        self.cosmoDir = cosmoDir
        self.prefix_Nvisits = prefix_Nvisits
        self.visitsDir = visitsDir
        self.budget = budget

        self.data = self.get_data(conf)

        print(self.data)

    def get_data(self, conf):
        """
        Method to load data in a pandas df

        Parameters
        --------------
        conf: str
          configuration file to process (csv file)

        Returns
        ----------
        pandas df of data

        """
        df = pd.DataFrame()
        for index, row in conf.iterrows():
            fi = '{}/{}.hdf5'.format(self.cosmoDir, row['cosmoFile'])
            nName = row['name']
            ttype = row['type']
            fload = self.load(fi)
            fload['confName'] = nName
            df = pd.concat((df, fload))

        df = self.add_infos(df)
        df = df.fillna(0)
        return df

    def load(self, fi):
        """
        Method  to load file in pandas df

        Parameters
        ---------------
        fi: str
          name of the file (full path) to be loaded

        Returns
        ----------
        pandas df corresponding to the file
        """
        df = pd.read_hdf(fi)

        # df['year'] = df['conf'].str.split('_').str.get(1).astype(int)
        df = df.sort_values(by=['year'])

        # print('ici', df)
        return df

    def add_infos(self, df):
        """
        Method to add infos to df

        Parameters
        ---------------
        df: pandas df
          original df to append

        Returns
        ----------
        pandas df with added values

        """
        Ny = '{}'.format(int(np.median(df['Ny'])))

        fName = '{}_{}.npy'.format(self.prefix_Nvisits, Ny)
        budget = DD_Budget(self.visitsDir, fName)
        # add budget colmuns
        cosmo_bud = budget(df)
        df = cosmo_bud
        df['frac_high_z'] = df['nsn_z_09']/df['nsn_DD']
        df['nsn_DD_ultra'] = df['nsn_DD_COSMOS'] + \
            df['nsn_DD_XMM-LSS']
        if 'nsn_ultra' in df.columns:
            df['nsn_ultra_zless_08'] = df['nsn_ultra'] - \
                df['nsn_ultra_z_08']
            df['nsn_dd_zless_05'] = df['nsn_dd']-df['nsn_dd_z_05']
            df['nsn_ultra_spectro'] = 200*df['year']
            df['nsn_dd_spectro'] = 200*df['year']
            df['nsn_spectroid'] = df[['nsn_ultra_zless_08', 'nsn_ultra_spectro']].min(axis=1) +\
                df[['nsn_dd_zless_05', 'nsn_dd_spectro']].min(axis=1)
            df['nsn_photid'] = df['nsn_DD']-df['nsn_spectroid']
            df['frac_photid'] = df['nsn_photid']/df['nsn_DD']
        df['sigma_w_rel'] = df['sigma_w']/df['w']
        return df

    def get_infos(self, xvar='year', yvar=['sigma_w', 'nsn_DD', 'nsn_DD_ultra'], budget=[0.05, 0.08, 0.10]):

        rtot = []

        print(self.data.columns)
        self.data = self.data.round(
            {'zcomp_dd_unique': 2, 'zcomp_ultra_unique': 2})
        self.data['zcomp_dd_str'] = self.data['zcomp_dd_unique'].astype(
            str).apply(
            lambda x: x.ljust(4, '0'))

        self.data['zcomp_ultra_str'] = self.data['zcomp_ultra_unique'].astype(
            str).apply(lambda x: x.ljust(4, '0'))
        self.data['confName'] = self.data['confName']+'_' + \
            self.data['zcomp_ultra_str']+'_'+self.data['zcomp_dd_str']
        print(self.data[['confName', 'zcomp_dd_str', 'zcomp_ultra_str']])

        self.data = self.modif_early()
        # confref = 'universal_0.00_0.65'
        confrefs = ['deep_rolling_early_0.80_0.70']
        # self.plot(confrefs)
        rr = self.data.groupby(['confName', 'zcomp_dd_str', 'zcomp_ultra_str']).apply(
            lambda x: self.ana_grp(x, yvar=yvar)).reset_index()

        return rr
        """
        print(rr)

        idx = rr['confName'] == confrefs[0]
        sel = rr[idx]
        print(sel[['time_0.05', 'sigma_w_0.05', 'nsn_DD_0.05']])
        """

    def plot(self, confNames):

        sel = pd.DataFrame()
        for confName in confNames:
            idx = self.data['confName'] == confName
            sel = pd.concat((sel, self.data[idx]))
        sel = sel.sort_values(by=['year'])
        fig, ax = plt.subplots()
        ax.plot(sel['year'], sel['budget'])

        plt.show()

    def modif_early(self):

        idx = self.data['confName'].str.contains("deep_rolling_early")

        sel = self.data[idx]
        df_new = pd.DataFrame(self.data[~idx])
        print('alors man', sel['confName'])
        for confName in sel['confName'].unique():
            ido = sel['confName'] == confName
            selb = sel[ido]
            z_ultra = selb['zcomp_ultra_str'].to_list()[0]
            z_dd = selb['zcomp_dd_str'].to_list()[0]
            if z_dd != '0.00':
                cconfName = '_'.join(confName.split('_')[:3])
                str_search = '{}_{}_0.00'.format(cconfName, z_ultra)
                print('searching string', str_search)
                idata = self.data['confName'] == str_search
                seldata = self.data[idata]
                if len(seldata) > 0:
                    print('yes man', confName, seldata)
                    dfcp = pd.DataFrame(seldata)
                    dfcp['confName'] = confName
                    dfcp['zcomp_ultra_str'] = z_ultra
                    dfcp['zcomp_dd_str'] = z_dd
                    dfcon = pd.concat((dfcp, selb))
                    df_new = pd.concat((df_new, dfcon))
        return df_new

    def ana_grp(self, grp, xvar='year', yvar=['sigma_w', 'nsn_DD']):

        resdict = {}

        print('hello man', grp.name)
        # budget interpolator
        interp_bud = interp1d(
            grp['budget'], grp[xvar], bounds_error=False, fill_value=0.)
        for i, val in enumerate(self.budget):
            ttime = interp_bud(val)
            tName = 'time_{}'.format(val)
            resdict[tName] = [ttime.tolist()]
            for yv in yvar:
                vName = '{}_{}'.format(yv, val)
                interp_var = interp1d(
                    grp[xvar], grp[yv], bounds_error=False, fill_value=0.)
                valres = interp_var(ttime)
                resdict[vName] = [valres.tolist()]

        print(grp.name, resdict)
        return pd.DataFrame.from_dict(resdict)


def plot_new(res, varx='zcomp', vary='sigma_w', budget=[0.05, 0.08, 0.10]):

    fig, ax = plt.subplots()

    io = res['zcomp'] > 0.
    sel = res[io]
    if len(sel) > 0:
        for bud in budget:
            iob = np.abs(sel['budget']-bud) < 1.e-5
            selb = sel[iob]
            ax.plot(selb[varx], selb[vary])

    # plt.show()


def get_infos(df, xvar='year', yvar=['sigma_w', 'nsn_DD'], budget=[0.05, 0.08, 0.10]):

    rtot = []
    for key, vals in df.items():
        r = []
        print('Survey', key, vals[['conf', 'year', 'sigma_w']])
        r.append(key)
        zcomp = -1
        if 'universal' in key:
            zcomp = key.split('_')[-1]
            zcomp = float(zcomp)
        r.append(zcomp)
        print('iii', vals['year'])
        interp_bud = interp1d(
            vals['budget'], vals[xvar], bounds_error=False, fill_value=0.)
        for i, val in enumerate(budget):
            ttime = interp_bud(val)
            print('eee', ttime)
            rp = list(r)
            rp += [val, ttime]
            for yv in yvar:
                interp_var = interp1d(
                    vals[xvar], vals[yv], bounds_error=False, fill_value=0.)
                valres = interp_var(ttime)
                print('resu budget', val, xvar, ttime, yv, valres)
                rp.append(valres)
                print('boo', rp)
            rtot.append(rp)

    return np.rec.fromrecords(rtot, names=['survey', 'zcomp', 'budget', 'time']+yvar)


def add_infos(df):

    for key, vals in df.items():

        Ny = '{}'.format(int(np.median(vals['Ny'])))

        fName = '{}_{}.npy'.format(prefix_Nvisits, Ny)
        budget = DD_Budget(visitsDir, fName)
        # add budget colmuns
        cosmo_bud = budget(vals)
        df[key] = cosmo_bud
        df[key]['frac_high_z'] = df[key]['nsn_z_09']/df[key]['nsn_DD']
        df[key]['nsn_DD_ultra'] = df[key]['nsn_DD_COSMOS'] + \
            df[key]['nsn_DD_XMM-LSS']
        if 'nsn_ultra' in df[key].columns:
            df[key]['nsn_ultra_zless_08'] = df[key]['nsn_ultra'] - \
                df[key]['nsn_ultra_z_08']
            df[key]['nsn_dd_zless_05'] = df[key]['nsn_dd']-df[key]['nsn_dd_z_05']
            df[key]['nsn_ultra_spectro'] = 200*df[key]['year']
            df[key]['nsn_dd_spectro'] = 200*df[key]['year']
            df[key]['nsn_spectroid'] = df[key][['nsn_ultra_zless_08', 'nsn_ultra_spectro']].min(axis=1) +\
                df[key][['nsn_dd_zless_05', 'nsn_dd_spectro']].min(axis=1)

            df[key]['nsn_photid'] = df[key]['nsn_DD']-df[key]['nsn_spectroid']
            df[key]['frac_photid'] = df[key]['nsn_photid']/df[key]['nsn_DD']
    return df


def plot(df, xvar='year', yvar='sigma_w', yvartwin='', legx='Year', legy='$\sigma_w$', tag_budget=[], tag_marker=[], zcomp=0.7):

    vva = ['deep_rolling_early', 'deep_rolling_ten_years',
           'universal_{}'.format(zcomp)]
    vvp = ['Early Deep Rolling (EDR$^3$)',
           'Deep Rolling 10 years (DR$^{10Y}$)', 'Deep Universal (DU$^{0.65}$)']
    colors = ['red', 'blue', 'magenta']
    corresp = dict(zip(vva, vvp))
    ccolors = dict(zip(vva, colors))
    fig, ax = plt.subplots(figsize=(15, 8))
    lls = ['solid', 'dashed', 'dotted', 'dotted', 'solid', 'dashed', ]
    keys = list(df.keys())
    ls = dict(zip(keys, lls[:len(keys)]))
    # for key, vals in df.items():
    for key in vva:
        vals = df[key]
        go_plot = True
        label = corresp[key]
        if 'universal' in key:
            idx = vals['year'] == 1
            zcomph = np.round(np.mean(vals[idx]['zcomp_dd'].to_list()), 2)
            label = label.replace('0.65', str(zcomp))
            if np.abs(zcomp-zcomph) >= 1.e-4:
                go_plot = False

        vals = vals.to_records(index=False)
        if go_plot:
            ax.plot(vals[xvar], vals[yvar], marker='o',
                    label=label, ls=ls[key], ms=5, color=ccolors[key])
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
            # print('hello', dfb[[varx, vary]])
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
                # print('boo', refval, key)
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

    # print('ici', df)
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

    # print(sel)

    return sel


def gime_data(conf, zcomps=[0.65, 0.70]):

    df = {}
    df_syste = {}
    for index, row in conf.iterrows():
        fi = '{}/{}.hdf5'.format(cosmoDir, row['cosmoFile'])
        nName = row['name']
        ttype = row['type']
        if 'universal' in nName:
            # df['universal_0.60'] = load_b(fi, zcomp=0.60)
            for zcomp in zcomps:
                nNameb = '{}_{}'.format(nName, zcomp)
                fload = load_b(fi, zcomp=zcomp)
                if ttype == 'nominal':
                    df[nNameb] = fload
                else:
                    df_syste[nNameb] = fload
        else:
            fload = load(fi)
            if ttype == 'nominal':
                df[nName] = fload
            else:
                df_syste[nName] = fload

    for key, vals in df.items():
        vals['sigma_w_rel'] = vals['sigma_w']/vals['w']

    return df, df_syste


def ana_res(infos, xvar='sigma_w', budget=0.05):

    print('hello', infos)
    vartot = '{}_{}'.format(xvar, budget)

    infos = infos.sort_values(by=[vartot])
    vvar = ['confName']+[vartot]
    print(infos[vvar])


def plot_info(df, xvar='zcomp_dd', yvar='sigma_w', xleg='$z_{complete}^{DD}$', yleg='$\sigma_w$', budget=0.05):

    yplot = '{}_{}'.format(yvar, budget)

    fig, ax = plt.subplots(figsize=(15, 8))
    fig.suptitle('Budget: {}'.format(budget))
    confPlot = ['deep_rolling_early_0.70',
                'deep_rolling_early_0.75',
                'deep_rolling_early_0.80',
                'universal_0.00', 'deep_rolling_ten_years_0.75']
    # namePlot = ['EDR$^3_{0.70}$',
    #            'EDR$^3_{0.75}$', 'EDR$^3_{0.80}$', 'DU', 'DR$^{10Y}_{0.75}$']
    namePlot = ['EDR$_{0.70}$',
                'EDR$_{0.75}$', 'EDR$_{0.80}$', 'DU', 'DR$_{0.75}$']
    corresp = dict(zip(confPlot, namePlot))
    markers = ['o', 's', '*', 'p', '>']
    mm = dict(zip(confPlot, markers))

    for confName in confPlot:
        idx = df['confName'].str.contains(confName)
        sel = df[idx].to_records(index=False)
        sel = sel[sel[yplot] > 1.e-5]
        print('aooo', sel[yplot], sel[xvar])
        ax.plot(sel[xvar], sel[yplot], marker=mm[confName],
                ls='None', ms=15, label=corresp[confName])

    ax.grid()
    ax.set_xlabel(xleg)
    ax.set_ylabel(yleg)
    ax.legend(loc='upper left')


parser = OptionParser()
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option("--cosmoDir", type=str, default='cosmo_files',
                  help="directory where cosmo files are located[%default]")
parser.add_option('--cadence', type=float, default=1.,
                  help='cadence of observation[%default]')
parser.add_option('--zcomp', type=float, default=0.70,
                  help='redshift completeness for universal survey [%default]')
parser.add_option('--config', type='str', default='config_survey_time.csv',
                  help='configuration file [%default]')
parser.add_option('--prefix_Nvisits', type='str', default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help='prefix for Nvisits file[%default]')

opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
visitsDir = opts.visitsDir
cadence = opts.cadence
prefix_Nvisits = opts.prefix_Nvisits
zcomp = opts.zcomp

bb = Budget_Time(opts.config, cosmoDir=cosmoDir,
                 prefix_Nvisits=prefix_Nvisits,
                 visitsDir=visitsDir, budget=[0.05, 0.06, 0.08, 0.10])
infos = bb.get_infos()

ana_res(infos)
print(infos.columns)
infos['zcomp_dd'] = infos['zcomp_dd_str'].astype(float)
infos['zcomp_ultra'] = infos['zcomp_ultra_str'].astype(float)
budget = 0.05
plot_info(infos, xvar='zcomp_dd', yvar='sigma_w',
          xleg='$z_{complete}^{DD}$', yleg='$\sigma_w$', budget=budget)
plot_info(infos, xvar='zcomp_dd', yvar='nsn_DD',
          xleg='$z_{complete}^{DD}$', yleg='$N_{SN}$', budget=budget)
plot_info(infos, xvar='zcomp_dd', yvar='time',
          xleg='$z_{complete}^{DD}$', yleg='Time budget [y]', budget=budget)
plot_info(infos, xvar='zcomp_dd', yvar='nsn_DD_ultra',
          xleg='$z_{complete}^{DD}$', yleg='$N_{SN}^{ultra}$', budget=budget)
plt.show()
print(test)

conf = pd.read_csv(opts.config, comment='#')
zcomps = [0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]
zcomps = [0.65, 0.7]
df, df_syste = gime_data(
    conf, zcomps=zcomps)
print('rrrr', df)
df = add_infos(df)
res = get_infos(df)

"""
plot_new(res, varx='zcomp', vary='time')
plot_new(res, varx='zcomp', vary='sigma_w')
plt.show()
print(res)
print(test)
"""
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


# estimate Budget for all of these
# budget class instance


# plot_multiple(df, yvartwin='budget')
# plot(df, yvar='w', legy='w')
# plot(df, yvar='Om', legy='Om')


plot(df, yvar='w', legy='$\sigma_w$', tag_budget=[
     0.05, 0.0788, 0.10], tag_marker=['*', '^', 's'], zcomp=zcomp)

plot(df, xvar='nsn_DD', legx='$N_{SN}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'], zcomp=zcomp)
"""
plot(df, xvar='frac_photid', legx='$frac_{photid}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])

plot(df, xvar='nsn_photid', legx='$N_{SN}^{photid}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])

plot(df, xvar='nsn_spectroid', legx='$N_{SN}^{spectroid}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])
"""
"""
plot(df, xvar='nsn_DD_ultra', legx='$N_{SN}^{ultra}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])
"""
"""
plot_syste(df, df_syste, yvar='sigma_w', tag_budget=[
    0.05, 0.0788, 0.10], legy='$\Delta\sigma_w$', tag_marker=['*', '^', 's'], norm=True)
"""

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
