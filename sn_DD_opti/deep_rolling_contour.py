from scipy.interpolate import interp1d
import numpy.lib.recfunctions as rf
import copy
from sn_rate import SN_Rate
import matplotlib.pyplot as plt
from optparse import OptionParser
import numpy as np
from scipy.ndimage import zoom
from scipy.ndimage.filters import gaussian_filter
from matplotlib.ticker import MaxNLocator

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['legend.fontsize'] = 15
#plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.family'] = 'Arial'


class DDF_Scenario:
    """
    class to estimate the DD budget for a set of scenarios

    Parameters
    ---------------
      zfields: dict
       dict of fields with fixed zlim,nseasons
    fDir: str
      location dir of (nvisits, zlim) file
    fName: str
      (nvisits, zlim) filename
    cadence: int, opt
      cadence of observation (default: 1 day)
    H0: float, opt
      H0 parameter (default: 70)
    Om0: float, opt
      Omega_0 parameter (default: 0.3)
    min_rf_phase: float, opt
      min_rf_phase (default: -15.)
    max_rf_phase: float, opt
      max_rf_phase (default: 30.)
    slDir: str, opt
      location dir of (seasonlength vs nvisits) file (default: input)
    slName: str, opt
      (seasonlength vs nvisits) (default: seasonlength_nvisits.py)
    """

    def __init__(self, zfields, fDir, fName,
                 cadence=1,
                 H0=72., Om0=0.30,
                 min_rf_phase=-15, max_rf_phase=30.,
                 slDir='input', slName='seasonlength_nvisits.npy',
                 plot_refs=False):

        self.zlim_visit, self.visit_zlim = self.load(fDir, fName)

        self.nvisits_seasonlength = self.load_seasonlength_nvisits(
            slDir, slName)

        if plot_refs:
            self.plot_zlim_visits(self.visit_zlim)
            self.plot_zlim_visits(
                self.zlim_visit, varx='nvisits', legx='Nvisits', legy='$z_{complete}$')
            self.plot_seasonlength_nvisits(self.nvisits_seasonlength)
            plt.show()

        self.rateSN = SN_Rate(H0=H0, Om0=Om0,
                              min_rf_phase=min_rf_phase,
                              max_rf_phase=max_rf_phase)
        print(vars(self.rateSN))

        self.cadence = cadence

        # get the number of visits and season length for fields with fixed zlim

        for vv in zfields.keys():
            zlim = zfields[vv]['zlimit']
            nvisits = self.visit_zlim[cadence](zlim)
            season_length = self.nvisits_seasonlength[vv](nvisits)

            zfields[vv]['Nvisits'] = nvisits.item()
            zfields[vv]['season_length'] = season_length.item(
            )/30.  # season length has to be in months!

        self.ultraFields = zfields

    def load(self, fDir, fName):
        """
        Method to load a file fDir/fName

        Parameters
        --------------
        fDir: str
          location dir of the file
        fName: str
          file name

        Returns
        -----------
        dict of interp1d
          key: cadence
          values: interp1d(nvisits, z)

        """

        from wrapper import Mod_z
        # tab = np.load(fName, allow_pickle=True)
        tab = Mod_z('{}/{}'.format(fDir, fName)).nvisits
        res = {}
        resb = {}
        print(tab)
        for cad in np.unique(tab['cadence']):
            idx = tab['cadence'] == cad
            sel = tab[idx]
            res[cad] = interp1d(sel['Nvisits'], sel['z'],
                                bounds_error=False, fill_value=0.)
            resb[cad] = interp1d(sel['z'], sel['Nvisits'],
                                 bounds_error=False, fill_value=0.)
        return res, resb

    def load_seasonlength_nvisits(self, slDir, slName):
        """
        Method to load (seasonlength, nvisits) file and make interp1d out of it

        Parameters
        ---------------
        slDir: str
          location dir of (seasonlength vs nvisits) file
        slName: str
          (seasonlength vs nvisits)

        Returns
        -----------
        dict of interp1d(nvisits, seasonlength); key: field name
        """

        tab = np.load('{}/{}'.format(slDir, slName), allow_pickle=True)

        res = {}
        for fieldName in np.unique(tab['name']):
            idx = tab['name'] == fieldName
            sel = tab[idx]
            res[fieldName] = interp1d(sel['nvisits'], sel['season_length'],
                                      bounds_error=False, fill_value=0.)
        return res

    def nights(self, season_length, nseasons, cadence):
        """
        Function to estimate the number of observing nights

        Parameters
        ---------------
        season_length: float
           season length (in months)
        nseasons: int
           number of seasons
        cadence: float
          cadence of observation

        Returns
        ----------
        number of observing nights

        """
        return season_length*nseasons*30./cadence

    def nvisits_field(self, nvisits, season_length, nseasons, cadence):
        """
        Method to estimate the total number of visits

        Parameters
        ---------------
        nvisits: int
          number of visits per observing night
        season_length: float
        season length (in months)
        nseasons: int
           number of seasons
        cadence: float
          cadence of observation
        Returns
        ----------
        total number of visits

        """
        return nvisits*self.nights(season_length, nseasons, cadence)

    def zlim(self, budget, Nfields, Nseasons=2, season_length=6.):

        Nvisits_nonDD = 2122176
        Nvisits_ultraDD = 0.
        for key, vals in self.ultraFields.items():
            print('io', key, vals)
            Nvisits_ultraDD += self.nvisits_field(
                vals['Nvisits'], vals['season_length'], vals['nseasons'], self.cadence)

        print('booo', Nvisits_ultraDD)
        Nvisits_DD = (budget*Nvisits_nonDD+(budget-1.)
                      * Nvisits_ultraDD)/(1.-budget)

        nvisits = Nvisits_DD/(Nfields*Nseasons*season_length*30.)/self.cadence
        print('aoooo', nvisits)
        zlimit = self.zlim_visit[self.cadence](nvisits)

        return np.round(zlimit, 4)

    def budget(self, Nfields, zlim, nseasons=2, season_length=6.):

        Nvisits_nonDD = 2122176

        Nvisits = self.visit_zlim[self.cadence](zlim)
        print('hhhhh', Nvisits)
        Nvisits_DD = Nfields*self.nvisits_field(
            Nvisits, season_length, nseasons, self.cadence)

        for key, vals in self.ultraFields.items():
            print('io', key, vals)
            Nvisits_DD += self.nvisits_field(
                vals['Nvisits'], vals['season_length'], vals['nseasons'], self.cadence)

        return Nvisits_DD/(Nvisits_DD+Nvisits_nonDD)

    def nvisits_tot(self, config):
        """
        Method to estimate the total number of visits for a config file

        Parameters
        ---------------
        config: dict
          configuration dict

        Returns
        -----------
        the total number of visits

        """
        bud = 0.

        for field in config['Fields']:
            if field != 'AGN':
                nnvisits = self.nvisits_field(config[field])
                bud += config[field]['npointings'] * nnvisits
                print('bud', field, nnvisits, bud)

        return bud

    def nsn(self, config, nADFS=1):
        """
        Method to estimate the number of supernovae corresponding to a config file

        Parameters
        ---------------
        config: dict
          configuration of the survey
        nADFS: int, opt
          number of pointings in the ADFS region (default: 1)

        Returns
        -----------
        The total number of supernovae

        """
        zmin = 0.01
        dz = 0.01
        survey_area = 9.6
        nsn_tot = 0
        for field in config['Fields']:
            zmax = config[field]['zlim']
            seasonlength = config[field]['seasonlength']*30.
            nseasons = config[field]['nseasons']

            nsn_field = self.nsn_indiv(
                zmin, zmax, dz, seasonlength, survey_area)
            print('yer', field, zmax, seasonlength, zmin, dz,
                  nseasons, nsn_field)
            if field != 'ADFS':
                nsn_tot += nsn_field*nseasons
            else:
                nsn_tot += nsn_field*nseasons*nADFS

        return nsn_tot

    def nsn_new(self, zlim, Nfields, nseasons=2, season_length=6.):

        zmin = 0.01
        dz = 0.01
        survey_area = 9.6
        nsn_tot = 0

        vec_nsn = np.vectorize(self.nsn_indiv)
        nsn_field = Nfields*nseasons*vec_nsn(
            zmin, zlim, dz, season_length*30., survey_area)

        # add ultradeep fields
        for key, vals in self.ultraFields.items():
            nsn_field += vals['nseasons']*self.nsn_indiv(
                zmin, np.round(vals['zlimit'], 2), dz, np.round(vals['season_length'], 1)*30., survey_area)

        return nsn_field

    def nsn_indiv(self, zmin, zmax, dz, seasonlength, survey_area):
        """
        Method to estimate the number of supernovae

        Parameters
        ---------------
        zmin: float
          min redshift
        zmax: float
         max redshift
        dz: float
          z binning
        seasonlength: float
         season length (in month)
        survey_area: float
          survey area (deg2)

        Returns
        ----------
        nsn: float, the total number of supernovae

        """

        zz, rate, err_rate, nsn_r, err_nsn = self.rateSN(zmin=zmin,
                                                         zmax=zmax,
                                                         dz=dz,
                                                         duration=seasonlength,
                                                         survey_area=survey_area,
                                                         account_for_edges=True)

        return np.cumsum(nsn_r)[-1]

    def plot_zlim_visits(self, tab, varx='zlim', legx='$z_{complete}$', legy='Nvisits'):
        """
        Method to plot zlim/nvisits vs nvisits/zlim

        Parameters
        ---------------
        tab: dict
          dict of 1d interpolator; key: cadence
        varx: str, opt
          variable to plot (default: zlim)
        legx: str, opt
          x-axis legend (default: z_complete)
        legy: str, opt
          y-axis legend (default: Nvisits)


        """

        fig, ax = plt.subplots(figsize=(14, 8))

        if varx == 'zlim':
            xvals = np.arange(0.5, 0.951, 0.01)
        if varx == 'nvisits':
            xvals = np.arange(10, 300, 10)

        for key, vals in tab.items():
            yvals = vals(xvals)
            ax.plot(xvals, yvals, label='cadence={}'.format(key))

        fontsize = 15
        ax.set_xlabel(legx, fontsize=fontsize)
        ax.set_ylabel(legy, fontsize=fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
        ax.grid()
        ax.legend()

    def plot_seasonlength_nvisits(self, tab):

        fig, ax = plt.subplots(figsize=(12, 8))
        nvisits = np.arange(10, 500, 10)

        for key, vals in tab.items():
            sl = vals(nvisits)
            ax.plot(nvisits, sl, label='{}'.format(key))

        fontsize = 15

        ax.set_xlabel('Nvisits', fontsize=fontsize)
        ax.set_ylabel('Season length [days]', fontsize=fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
        ax.tick_params(axis='y', labelsize=fontsize)
        ax.grid()
        ax.legend()


def plotContourBudget(zfields, fDir, slDir='input', slName='seasonlength_nvisits.npy'):

    #fig, ax = plt.subplots(nrows=2, figsize=(6, 8))
    leg = ''
    for kk in zfields.keys():
        leg += kk+'$^{'+str(zfields[kk]['zlimit']) + \
            '}_{'+str(zfields[kk]['nseasons'])+'}$'
        if kk == 'COSMOS':
            leg += ' + '
    """
    fig, ax = plt.subplots(figsize=(8, 6))
 
    fig.suptitle('N$_{\mathrm{visits}}^{y}$ = 80 - '+leg, weight='bold')
    fName = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_80.npy'
    plotContour(ax, zfields, fDir, fName,
                slDir=slDir, slName=slName, color='b')
    ax.set_xlabel('Number of deep fields$_2$')
    ax.set_ylabel('\mathrm{$z_{complete}$}')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    """
    figb, axb = plt.subplots(figsize=(8, 6))
    figb.suptitle('N$_{\mathrm{visits}}^{y}$ = 20 - '+leg)
    fName = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_20.npy'
    plotContour(axb, zfields, fDir, fName,
                slDir=slDir, slName=slName, color='k')

    axb.set_xlabel('Number of deep fields$_{\mathrm{2}}$')
    axb.set_ylabel('$\mathrm{z_{complete}}$')
    axb.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()


def plotContour(ax, zfields, fDir, fName, slDir='input', slName='seasonlength_nvisits.npy', color='k'):
    """
    Method to get results from the DDF_Scenario class

    Parameters
    ---------------
    zfields: dict
      dict of fields with fixed zlim
    fDir: str
      location dir of (nvisits, zlim) file
    fName: str
      (nvisits, zlim) filename
    slDir: str, opt
      location dir of (seasonlength vs nvisits) file (default: input)
    slName: str, opt
      (seasonlength vs nvisits) (default: seasonlength_nvisits.py)

    Returns
    ----------
    numpy arrays with results

    """

    dd_scen = DDF_Scenario(zfields, fDir, fName,
                           slDir=slDir, slName=slName)

    budmin = 0.06
    budmax = 0.15
    nfmin = 1
    nfmax = 4
    zmin = 0.5
    zmax = 0.9

    budget = np.linspace(budmin, budmax, 1000)
    nfields = np.linspace(nfmin, nfmax, 60)
    zlim = np.linspace(zmin, zmax, 100)
    NF, ZLIMIT = np.meshgrid(nfields, zlim)
    BUD = dd_scen.budget(NF, ZLIMIT)
    NSN = dd_scen.nsn_new(ZLIMIT, NF)
    print(NSN)
    ax.imshow(BUD, extent=(nfmin, nfmax, zmin, zmax),
              aspect='auto', alpha=0.25, cmap='hsv')
    # ax[1].imshow(NSN, extent=(nfmin, nfmax, zmin, zmax),
    #             aspect='auto', alpha=0.25, cmap='hsv')

    zzv = [0.5, 0.6, 0.7, 0.8, 0.9]
    zzv = [0.07, 0.08, 0.10, 0.12, 0.15]
    CS = ax.contour(NF, ZLIMIT, BUD, zzv, colors='k')

    fmt = {}
    strs = ['$%3.2f$' % zz for zz in zzv]
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CS.levels, strs):
        fmt[l] = s
    ax.clabel(CS, inline=True, fontsize=10,
              colors='k', fmt=fmt)

    #axb = ax.twinx()
    zzv = [800, 1000, 1200, 1500, 2000, 2500, 3000]
    CSb = ax.contour(NF, ZLIMIT, gaussian_filter(
        NSN, sigma=1.5), zzv, colors='r')

    fmt = {}
    strs = ['$%1i$' % zz for zz in zzv]
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CSb.levels, strs):
        fmt[l] = s
    ax.clabel(CSb, inline=True, fontsize=10,
              colors='r', fmt=fmt)


parser = OptionParser()

parser.add_option('--ultraDeep', type=str, default='COSMOS,XMM-LSS',
                  help='list of ultra deep fields[%default]')
parser.add_option('--z_ultra', type=str, default='0.9,0.9',
                  help='z complete of ultra deep fields[%default]')
parser.add_option('--nseason_ultra', type=str, default='2,2',
                  help='number of visits for ultra deep fields[%default]')
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")

opts, args = parser.parse_args()

ultraDeepFields = opts.ultraDeep.split(',')
z_ultra = list(map(float, opts.z_ultra.split(',')))
nseasons_ultra = list(map(int, opts.nseason_ultra.split(',')))

zfields = {}

for i, vv in enumerate(ultraDeepFields):
    zfields[vv] = {}
    zfields[vv]['zlimit'] = z_ultra[i]
    zfields[vv]['nseasons'] = nseasons_ultra[i]

print(zfields)

plotContourBudget(zfields, opts.visitsDir)
