import pandas as pd
from __init__ import plt, linestyles
import numpy as np
from optparse import OptionParser
import os
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

    def get_infos(self, xvar='year', yvar=['sigma_w', 'nsn_DD', 'nsn_DD_ultra', 'nsn_spectro_ultra_yearly'], budget=[0.05, 0.08, 0.10]):

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
        rt = []
        interp_bud = interp1d(
            grp['budget'], grp[xvar], bounds_error=False, fill_value=0.)
        for i, val in enumerate(self.budget):
            ttime = interp_bud(val)
            tName = 'time_{}'.format(val)
            resdict[tName] = [ttime.tolist()]
            ro = [val, ttime.tolist()]
            for yv in yvar:
                vName = '{}_{}'.format(yv, val)
                interp_var = interp1d(
                    grp[xvar], grp[yv], bounds_error=False, fill_value=0.)
                valres = interp_var(ttime)
                resdict[vName] = [valres.tolist()]
                ro.append(valres.tolist())
            rt.append(ro)
        print(grp.name, resdict)
        resdf = pd.DataFrame(rt, columns=['budget', 'time']+yvar)
        # return pd.DataFrame.from_dict(resdict)
        return resdf


def load_data(fi):

    df = pd.read_hdf(fi)

    return df


def plot(df, xvar='nsn_spectro_ultra_yearly', yvar='sigma_w', xleg='$N_{host}^{spectro}$/year', yleg='$\sigma_w$', budget=0.05):

    fig, ax = plt.subplots(figsize=(12, 8))
    # fig.subplots_adjust(hspace=0, wspace=0)
    fig.subplots_adjust(bottom=0.15)
    fig.suptitle('Budget: {}'.format(budget))
    confPlot = ['deep_rolling_early_0.70_0.60',
                'deep_rolling_early_0.75_0.60',
                'deep_rolling_early_0.80_0.60',
                'universal_0.00_0.65', 'deep_rolling_ten_years_0.75_0.65']
    llsty = ['solid', 'solid', 'solid', 'dashed', 'dotted']
    lsty = dict(zip(confPlot, llsty))
    # namePlot = ['EDR$^3_{0.70}$',
    #            'EDR$^3_{0.75}$', 'EDR$^3_{0.80}$', 'DU', 'DR$^{10Y}_{0.75}$']
    namePlot = ['EDR$_{0.70}$',
                'EDR$_{0.75}$', 'EDR$_{0.80}$', 'DU', 'DR$_{0.75}$']
    corresp = dict(zip(confPlot, namePlot))
    markers = ['o', 's', '*', 'p', '>']
    mm = dict(zip(confPlot, markers))

    idm = np.abs(df['budget']-budget) < 1.e-5
    dfb = df[idm]
    for confName in confPlot:
        idx = dfb['confName'].str.contains(confName)
        sel = dfb[idx].to_records(index=False)
        sel = sel[sel[yvar] > 1.e-5]
        print('aooo', sel[yvar], sel[xvar])
        ax.plot(sel[xvar], sel[yvar], marker=mm[confName],
                ls=lsty[confName], ms=15, label=corresp[confName], mfc='None')

    ax.grid()
    ax.set_xlabel(xleg)
    ax.set_ylabel(yleg)
    ax.legend(ncol=3, frameon=False)
    # fig.tight_layout()


def plot_info(df, xvar='zcomp_dd', yvar='sigma_w', xleg='$z_{complete}^{DD}$', yleg='$\sigma_w$', budget=0.05):

    fig, ax = plt.subplots(figsize=(12, 8))
    fig.suptitle('Budget: {}'.format(budget))
    fig.subplots_adjust(bottom=0.15)

    confPlot = ['deep_rolling_early_0.70',
                'deep_rolling_early_0.75',
                'deep_rolling_early_0.80',
                'universal_0.00', 'deep_rolling_ten_years_0.75']
    llsty = ['solid', 'solid', 'solid', 'dashed', 'dotted']
    lsty = dict(zip(confPlot, llsty))
    # namePlot = ['EDR$^3_{0.70}$',
    #            'EDR$^3_{0.75}$', 'EDR$^3_{0.80}$', 'DU', 'DR$^{10Y}_{0.75}$']
    namePlot = ['EDR$_{0.70}$',
                'EDR$_{0.75}$', 'EDR$_{0.80}$', 'DU', 'DR$_{0.75}$']
    corresp = dict(zip(confPlot, namePlot))
    markers = ['o', 's', '*', 'p', '>']
    mm = dict(zip(confPlot, markers))

    idm = np.abs(df['budget']-budget) < 1.e-5
    dfb = df[idm]
    for confName in confPlot:
        idx = dfb['confName'].str.contains(confName)
        sel = dfb[idx].to_records(index=False)
        sel = sel[sel[yvar] > 1.e-5]
        print('aooo', sel[yvar], sel[xvar])
        ax.plot(sel[xvar], sel[yvar], marker=mm[confName],
                ls=lsty[confName], ms=15, label=corresp[confName], mfc='None')

    ax.grid()
    ax.set_xlabel(xleg)
    ax.set_ylabel(yleg)
    ax.legend(ncol=3, frameon=False)


def plot_vs_budget(df, xvar='budget', yvar='sigma_w', xleg='budget', yleg='$\sigma_w$',
                   surveyName=['deep_rolling_early_0.80_0.60', 'deep_rolling_early_0.75_0.60', 'deep_rolling_early_0.70_0.60',
                               'deep_rolling_ten_years_0.75_0.65', 'universal_0.00_0.65', ],
                   plotName=['EDR$_{0.80}^{0.60}$', 'EDR$_{0.75}^{0.60}$',
                             'EDR$_{0.70}^{0.60}$',  'DR$_{0.75}^{0.65}$', 'DU$^{0.65}$'],
                   lls=['solid', linestyles['densely dashdotdotted'],
                        'dashdot', 'dotted', 'dashed'],
                   colors=['red', 'red', 'red', 'magenta', 'blue']):

    corresp = dict(zip(surveyName, plotName))
    ccolors = dict(zip(surveyName, colors))
    ls = dict(zip(surveyName, lls))
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.subplots_adjust(bottom=0.12, top=0.8)
    """
    confPlot = ['deep_rolling_early_0.70_0.60',
                'deep_rolling_early_0.75_0.60',
                'deep_rolling_early_0.80_0.60',
                'universal_0.00_0.65', 'deep_rolling_ten_years_0.75_0.65']
    llsty = ['solid', 'solid', 'solid', 'dashed', 'dotted']
    lsty = dict(zip(confPlot, llsty))
    # namePlot = ['EDR$^3_{0.70}$',
    #            'EDR$^3_{0.75}$', 'EDR$^3_{0.80}$', 'DU', 'DR$^{10Y}_{0.75}$']
    namePlot = ['EDR$_{0.70}^{0.60}$',
                'EDR$_{0.75}^{0.60}$', 'EDR$_{0.80}^{0.60}$', 'DU$^{0.65}$', 'DR$_{0.75}^{0.65}$']

    

    
    corresp = dict(zip(confPlot, namePlot))
    markers = ['o', 's', '*', 'p', '>']
    mm = dict(zip(confPlot, markers))
    """
    for confName in surveyName:
        idx = df['confName'].str.contains(confName)
        sel = df[idx].to_records(index=False)
        sel = sel[sel[yvar] > 1.e-5]
        sel = sel[sel[xvar] > 1.e-5]
        if confName == 'deep_rolling_ten_years_0.75_0.65':
            print('aooo', sel[yvar], sel[xvar])
        ax.plot(sel[xvar], sel[yvar], marker='None',
                ls=ls[confName], label=corresp[confName], color=ccolors[confName], lw=2)

    ax.grid()
    ax.set_xlabel(xleg)
    ax.set_ylabel(yleg)
    #ax.legend(ncol=3, frameon=False)
    ax.legend(loc='upper left', bbox_to_anchor=(
        0., 1.3), ncol=3, frameon=False)
    # fig.tight_layout()


def plot_vs_nspectro(config, visitsDir, prefix_Nvisits, nspectro, cosmoDir):

    totdf = pd.DataFrame()

    for nsp in nspectro:
        df = load_Summary(config, cosmoDir, prefix_Nvisits, visitsDir, nsp)
        totdf = pd.concat((totdf, df))

        print(totdf.columns)

        plot(totdf, budget=budget)


def load_Summary(config, cosmoDir, prefix_Nvisits, visitsDir, nsp):

    dirFile = '{}_{}'.format(cosmoDir, nsp)
    fName = '{}/budget_summary.hdf5'.format(dirFile)
    print('loading', fName)
    if not os.path.isfile(fName):
        process_Summary(config, dirFile, prefix_Nvisits, visitsDir, fName)
    df = load_data(fName)
    df['zcomp_dd'] = df['zcomp_dd_str'].astype(float)
    df['zcomp_ultra'] = df['zcomp_ultra_str'].astype(float)
    return df


def process_Summary(config, cosmoDir, prefix_Nvisits, visitsDir, outName):

    bb = Budget_Time(config, cosmoDir=cosmoDir,
                     prefix_Nvisits=prefix_Nvisits,
                     visitsDir=visitsDir, budget=np.arange(0.03, 0.11, 0.01))
    infos = bb.get_infos()

    # outName = '{}/budget_summary.hdf5'.format(cosmoDir)
    infos.to_hdf(outName, key='summary')


parser = OptionParser()
parser.add_option("--cosmoDir", type=str, default='cosmo_files_yearly_zspectro_x1c_1sigma_uy',
                  help="directory where cosmo files are located[%default]")
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option('--config', type='str', default='config_survey_time.csv',
                  help='configuration file [%default]')
parser.add_option('--nspectro', type='str', default='200,300,400,500,600,700,800',
                  help='nspectro host per year [%default]')
parser.add_option("--budget", type=float, default=0.05,
                  help="budget [%default]")
parser.add_option('--prefix_Nvisits', type='str', default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help='prefix for Nvisits file[%default]')

opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
config = opts.config
nspectro = opts.nspectro.split(',')
budget = opts.budget
visitsDir = opts.visitsDir
prefix_Nvisits = opts.prefix_Nvisits

# reading the configuration file
conf = pd.read_csv(config, comment='#')

# make a big df from all the files
if len(nspectro) > 1:
    plot_vs_nspectro(config, visitsDir, prefix_Nvisits, nspectro, cosmoDir)
else:
    df = load_Summary(config, cosmoDir, prefix_Nvisits, visitsDir, nspectro[0])

    # plot_vs_budget(df, xvar='budget', yvar='nsn_DD',
    #               xleg='budget', yleg='$N_{SN}$')

    plot_vs_budget(df, yvar='budget', xvar='nsn_DD',
                   yleg='budget', xleg='$N_{SN}$')

    # plot_vs_budget(df, yvar='time', xvar='budget',
    #                yleg='Time budget [year]', xleg='budget')

    plot_vs_budget(df, xvar='time', yvar='budget',
                   xleg='Time budget [year]', yleg='budget')

    """
    plot_vs_budget(df, xvar='budget', yvar='sigma_w',
                   xleg='budget', yleg='$\sigma_w$')
    
    print(df.columns)
    plot_info(df, xvar='zcomp_dd', yvar='sigma_w',
              xleg='$z_{complete}^{deep}$', yleg='$\sigma_w$', budget=budget)
    """
    """
    plot_info(df, xvar='zcomp_dd', yvar='nsn_DD',
              xleg='$z_{complete}^{DD}$', yleg='$N_{SN}$', budget=budget)
    plot_info(df, xvar='zcomp_dd', yvar='time',
              xleg='$z_{complete}^{DD}$', yleg='Time budget [y]', budget=budget)
    plot_info(df, xvar='zcomp_dd', yvar='nsn_DD_ultra',
              xleg='$z_{complete}^{DD}$', yleg='$N_{SN}^{ultra}$', budget=budget)
    plot_info(df, xvar='zcomp_dd', yvar='nsn_spectro_ultra_yearly',
              xleg='$z_{complete}^{DD}$', yleg='$N_{SN}^{ultra}$', budget=budget)
    """

plt.show()
