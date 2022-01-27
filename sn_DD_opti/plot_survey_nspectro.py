import pandas as pd
from __init__ import plt
import numpy as np
from optparse import OptionParser


def load_data(fi):

    df = pd.read_hdf(fi)

    return df


def plot(df, xvar='nsn_spectro_ultra_yearly', yvar='sigma_w', xleg='$N_{host}^{spectro}$/year', yleg='$\sigma_w$', budget=0.05):

    fig, ax = plt.subplots(figsize=(12, 8))
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


parser = OptionParser()
parser.add_option("--cosmoDir", type=str, default='cosmo_files_yearly_zspectro_x1c_1sigma_uy',
                  help="directory where cosmo files are located[%default]")
parser.add_option('--config', type='str', default='config_survey_time.csv',
                  help='configuration file [%default]')
parser.add_option('--nspectro', type='str', default='200,300,400,500,600,700,800',
                  help='nspectro host per year [%default]')
parser.add_option("--budget", type=float, default=0.05,
                  help="budget [%default]")

opts, args = parser.parse_args()

cosmoDir = opts.cosmoDir
config = opts.config
nspectro = opts.nspectro.split(',')
budget = opts.budget


# reading the configuration file
conf = pd.read_csv(config, comment='#')

# make a big df from all the files
totdf = pd.DataFrame()

for nsp in nspectro:
    dirFile = '{}_{}'.format(cosmoDir, nsp)
    fName = '{}/budget_summary.hdf5'.format(dirFile)
    print('loading', fName)
    df = load_data(fName)
    totdf = pd.concat((totdf, df))

print(totdf.columns)

plot(totdf, budget=budget)

plt.show()
