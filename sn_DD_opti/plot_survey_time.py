import pandas as pd
from __init__ import plt, linestyles
import numpy as np
from optparse import OptionParser
from wrapper import DD_Budget
from scipy.interpolate import interp1d, make_interp_spline, UnivariateSpline
from scipy.ndimage.filters import gaussian_filter


class CosmoData:
    def __init__(self, config, cosmoDir='cosmo_files_yearly',
                 prefix_Nvisits='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                 visitsDir='visits_files'):

        conf = pd.read_csv(config, comment='#')
        self.cosmoDir = cosmoDir
        self.prefix_Nvisits = prefix_Nvisits
        self.visitsDir = visitsDir

        self.data = self.get_data(conf)

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
        df = df.round(
            {'zcomp_dd_unique': 2, 'zcomp_ultra_unique': 2})
        df = df.fillna(0.0)
        print('hello', df[['zcomp_dd_unique', 'zcomp_ultra_unique']])
        df['zcomp_dd_str'] = df['zcomp_dd_unique'].astype(
            str).apply(
            lambda x: x.ljust(4, '0'))

        df['zcomp_ultra_str'] = df['zcomp_ultra_unique'].astype(
            str).apply(lambda x: x.ljust(4, '0'))
        df['confName'] = df['confName']+'_' + \
            df['zcomp_ultra_str']+'_'+df['zcomp_dd_str']
        df = self.modif_early(df)
        return df

    def modif_early(self, data):

        idx = data['confName'].str.contains("deep_rolling_early")

        sel = data[idx]
        df_new = pd.DataFrame(data[~idx])
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
                idata = data['confName'] == str_search
                seldata = data[idata]
                if len(seldata) > 0:
                    print('yes man', confName, seldata)
                    dfcp = pd.DataFrame(seldata)
                    dfcp['confName'] = confName
                    dfcp['zcomp_ultra_str'] = z_ultra
                    dfcp['zcomp_dd_str'] = z_dd
                    dfcon = pd.concat((dfcp, selb))
                    df_new = pd.concat((df_new, dfcon))
        return df_new


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


"""
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
"""


def plot(df, xvar='year', yvar='sigma_w', legx='Year', legy='$\sigma_w$',
         surveyName=['deep_rolling_early_0.80_0.60', 'deep_rolling_early_0.75_0.60', 'deep_rolling_early_0.70_0.60',
                     'deep_rolling_ten_years_0.75_0.65', 'universal_0.00_0.65', ],
         plotName=['EDR$_{0.80}^{0.60}$', 'EDR$_{0.75}^{0.60}$',
                   'EDR$_{0.70}^{0.60}$',  'DR$_{0.75}^{0.65}$', 'DU$^{0.65}$'],
         lls=['solid', linestyles['densely dashdotdotted'],
              'dashdot', 'dotted', 'dashed'],
         colors=['red', 'red', 'red', 'magenta', 'blue'], tag_budget=[], tag_marker=[], smooth_it=True, figtitle=''):

    corresp = dict(zip(surveyName, plotName))
    print('all', corresp)
    ccolors = dict(zip(surveyName, colors))
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle(figtitle)
    ls = dict(zip(surveyName, lls))
    # for key, vals in df.items():
    for sName in surveyName:
        idx = df['confName'] == sName

        # get budget interpolator
        interp_bud, interp_var = interp_budget(
            df, xvar, yvar, None, None, smooth_it=False)
        # get time budget=11%
        ttime = interp_bud(0.11)
        # remove points after ttime
        idx &= df[xvar] <= ttime
        vals = df[idx]
        if not smooth_it:
            ax.plot(vals[xvar], vals[yvar], marker='None',
                    label=corresp[sName], ls=ls[sName], ms=5, color=ccolors[sName], lw=2)

        xnew = None
        spl_smooth = None
        if smooth_it:
            xmin, xmax = np.min(vals[xvar]), np.max(vals[xvar])
            xnew = np.linspace(xmin, xmax, 100)
            spl = make_interp_spline(
                vals[xvar], vals[yvar], k=1)  # type: BSpline
            spl = UnivariateSpline(vals[xvar], vals[yvar], k=5)
            spl.set_smoothing_factor(0.5)
            spl_smooth = spl(xnew)
            ax.plot(xnew, spl_smooth, marker='None',
                    label=corresp[sName], ls=ls[sName], ms=5, color=ccolors[sName], lw=2)

        if len(tag_budget) > 0:
            interp_bud, interp_var = interp_budget(
                vals, xvar, yvar, xnew, spl_smooth, smooth_it)
            for i, val in enumerate(tag_budget):
                ttime = interp_bud(val)
                valres = interp_var(ttime)
                ax.plot(ttime, valres,
                        marker=tag_marker[i], color='k', ms=7, markeredgewidth=2)
                print('resu budget', val, xvar, ttime, yvar, valres)

    if xvar == 'year':
        ax.set_xlim([1., 6.])
    if yvar == 'sigma_w':
        ax.set_ylim([0.010, 0.05])
    if yvar == 'detfom':
        ax.set_ylim([40, None])
    ax.grid()
    ax.legend(ncol=3, frameon=False)
    ax.set_xlabel(legx)
    ax.set_ylabel(legy)


def interp_budget(vals, xvar, yvar, xnew, ynew, smooth_it):

    interp_bud = interp1d(
        vals['budget'], vals[xvar], bounds_error=False, fill_value=0.)
    if not smooth_it:
        interp_var = interp1d(
            vals[xvar], vals[yvar], bounds_error=False, fill_value=0.)
    if smooth_it:
        interp_var = interp1d(
            xnew, ynew, bounds_error=False, fill_value=0.)

    return interp_bud, interp_var


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
                            marker=tag_marker[i], color='k', ms=10, markeredgewidth=2)
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

    ax.legend(frameon=False)
    ax.grid()


parser = OptionParser()
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option("--cosmoDir", type=str, default='cosmo_files_yearly_zspectro_x1c_1sigma_uy',
                  help="directory where cosmo files are located[%default]")
parser.add_option('--cadence', type=float, default=1.,
                  help='cadence of observation[%default]')
parser.add_option('--zcomp', type=float, default=0.70,
                  help='redshift completeness for universal survey [%default]')
parser.add_option('--config', type='str', default='config_survey_time.csv',
                  help='configuration file [%default]')
parser.add_option('--prefix_Nvisits', type='str', default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help='prefix for Nvisits file[%default]')
parser.add_option('--nspectro', type='str', default='200',
                  help='nspectro tag [%default]')
parser.add_option('--figtitle', type='str', default='$N_{host}^{spectro}$=200/year',
                  help='plot title [%default]')


opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
visitsDir = opts.visitsDir
cadence = opts.cadence
prefix_Nvisits = opts.prefix_Nvisits
zcomp = opts.zcomp
config = opts.config
nspectro = opts.nspectro
config = opts.config
figtitle = opts.figtitle

# get the data
cDir = '{}_{}'.format(cosmoDir, nspectro)
df = CosmoData(config, cDir, prefix_Nvisits, visitsDir).data

print('hello', df['confName'])
print(df.columns)
"""
df['detfomb'] = 1./(df['sigma_wp']*df['sigma_wa'])
df['detfomb'] /= 6.17
df['detfom'] /= 6.17
plot(df, tag_budget=[0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
"""
CL = 6.17
df['detfom'] /= CL
plot(df, yvar='detfom', legy='DET FoM [95$\%$ C.L.]', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
"""
plot(df, yvar='nsn_ultra', legy='$N_{SN}^{ultra}$', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
plot(df, yvar='nsn_dd', legy='$N_{SN}^{deep}$', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
"""
"""
plot(df, yvar='sigma_w', legy='$\sigma_{w}$', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
plot(df, yvar='fom', legy='FoM', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
plot(df, yvar='detfom', legy='DetFoM', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
plot(df, yvar='detfomb', legy='DetFoMb', tag_budget=[
     0.05, 0.0788], tag_marker=['o', 's'], figtitle=figtitle)
"""
plt.show()

"""
plot(df, yvar='w', legy='$\sigma_w$', tag_budget=[
     0.05, 0.0788, 0.10], tag_marker=['*', '^', 's'], zcomp=zcomp)

plot(df, xvar='nsn_DD', legx='$N_{SN}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'], zcomp=zcomp)

plot(df, xvar='nsn_DD_ultra', legx='$N_{SN}^{ultra}$', tag_budget=[
     0.05, 0.0785,  0.10], tag_marker=['*', '^', 's'])
"""
