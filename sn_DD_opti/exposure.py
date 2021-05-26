import yaml
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
import numpy.lib.recfunctions as rf
import copy
from sn_rate import SN_Rate

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.family'] = 'serif'


class DDF_Scenario:
    """
    class to estimate the DD budget for a set of scenarios

    Parameters
    ---------------
    config: dict
       initial configuration of the survey
      zfields: dict
       dict of fields with fixed zlim
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

    def __init__(self, config, zfields, fDir, fName,
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

        # get the number of visits for the "z fixed" fields
        self.update_config_zlim(config, zfields)
        nvisits_DDF = self.nvisits_tot(config)
        print('estimated budget', nvisits_DDF /
              (config['nvisits_WFD']+nvisits_DDF))

        configb = copy.deepcopy(config)
        self.res = self.budget_multiple_config(config, nADFS=1)
        self.resb = self.budget_multiple_config(config, nADFS=2)

    def update_config_zlim(self, config, zfields):
        """
        Method to update config file

        Parameters
        ---------------
        config: dict
          survey configuration dict
        fields: dict
           dict of fields with zlim values
        """

        for field, zlim in zfields.items():
            nvisits = int(self.visit_zlim[self.cadence](zlim))
            config[field]['nvisits'] = nvisits
            seasonlength = self.nvisits_seasonlength[field](nvisits)
            seasonlength = np.min([seasonlength, 180.])
            config[field]['seasonlength'] = np.round(
                seasonlength.item()/30., 1)
            config[field]['cadence'] = self.cadence

    def update_config_nvisits(self, config, fields=['COSMOS', 'XMM-LSS', 'ELAIS', 'CDFS', 'ADFS']):
        """
        Method to update config file

        Parameters
        ---------------
        config: dict
          survey configuration dict
        fields: list(str), opt
           list of fields to consider (default: COSMOS,XMM-LSS,ELAIS,CDFS,ADFS)
        """

        for field in fields:
            nvisits = config[field]['nvisits']
            config[field]['zlim'] = self.zlim_visit[self.cadence](
                nvisits).item()
            seasonlength = self.nvisits_seasonlength[field](nvisits)
            seasonlength = np.min([seasonlength, 180.])
            config[field]['seasonlength'] = np.round(
                seasonlength.item()/30., 1)
            config[field]['cadence'] = self.cadence

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

    def nights(self, field):
        """
        Function to estimate the number of observing nights

        Parameters
        ---------------
        field: dict
        data for the field

        Returns
        ----------
        number of observing nights

        """
        return field['seasonlength']*field['nseasons']*30./field['cadence']

    def nvisits_field(self, field):
        """
        Method to estimate the total number of visits

        Parameters
        ---------------
        field: dict
        data for the field

        Returns
        ----------
        total number of visits

        """
        return field['nvisits']*self.nights(field)

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

        return

    def budget_config(self, config, fields=['CDFS', 'ELAIS', 'ADFS'], nADFS=1):
        """
        Method to estimate the budget depending on the number of visits per obs night

        Parameters
        ---------------
        config: dict
          configuration dict
        fields: list(str)
        list of fields with varying number of visits per ight
        nADSF: int, opt
          number of pointing for the ADFS field

        Returns
        -----------
        numpy array with three cols: nvisits, zlim, budget
        """
        r = []
        print(config)
        for nv in range(8, 250, 10):
            for ff in fields:
                config[ff]['nvisits'] = nv

            config['ADFS']['npointings'] = nADFS
            self.update_config_nvisits(config, fields)
            zlim_deep = config[fields[0]]['zlim']
            nvisits_DDF = self.nvisits_tot(config)
            budget = nvisits_DDF/(config['nvisits_WFD']+nvisits_DDF)
            print('there man', config, budget)
            r.append((nv, zlim_deep, self.nsn(config, nADFS), budget))

        res = np.rec.fromrecords(r, names=['nvisits', 'zlim', 'nsn', 'budget'])
        return res

    def budget_multiple_config(self, config, nADFS=1):
        """
        Function to estimate the budget depending on multiple configurations

        Parameters
        ---------------
        config: dict
          configuration dict
        nADSF: int, opt
          number of pointing for the ADFS field

        Returns
        -----------
        numpy array

        """

        fields_base = ['COSMOS', 'XMM-LSS']
        fields = [['ADFS'], ['ADFS', 'CDFS'], ['CDFS', 'ELAIS', 'ADFS']]
        rtot = None
        for ff in fields:
            config['Fields'] = ff+fields_base
            rb = self.budget_config(config, ff, nADFS)
            print(rb)
            rb = rf.append_fields(rb, 'Nf', [len(ff)]*len(rb))
            if rtot is None:
                rtot = rb
            else:
                rtot = np.concatenate((rtot, rb))

        return rtot

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


def getRes(config, zfields, fDir, fName, slDir='input', slName='seasonlength_nvisits.npy'):
    """
    Method to get results from the DDF_Scenario class

    Parameters
    ---------------
    config: dict
       initial configuration of the survey
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

    dd_scen = DDF_Scenario(config, zfields, fDir, fName,
                           slDir=slDir, slName=slName)

    res = dd_scen.res
    resb = dd_scen.resb

    return res, resb


def sel(res):

    idx = res['zlim'] <= 0.91
    sel = res[idx]
    sel.sort(order='zlim')
    return sel


def plotFinal(config, zfields, res, resb):

    fig, ax = plt.subplots(nrows=2, ncols=1, sharex='col', figsize=(14, 10))
    fig.subplots_adjust(hspace=0)
    color = dict(zip(range(1, 4), ['r', 'orange', 'b']))
    color = dict(zip(range(1, 5), ['b', 'r', 'orange', 'g']))
    nseas_adfs = config['ADFS']['nseasons']
    nseas_field = config['CDFS']['nseasons']

    legm = []
    for field, val in zfields.items():
        nseasons = str(config[field]['nseasons'])
        zlim = '$^{'+str(val)+'}_{'+nseasons+'}$'
        legm.append('{}{}'.format(field, zlim))

    legm = '+'.join(legm)
    ADFS = 'ADFS$_'+str(nseas_adfs)+'$'
    ADFS = 'Euclid/Roman$_'+str(nseas_adfs)+'$'

    field = 'field$_'+str(nseas_field)+'$'
    fields = 'fields$_'+str(nseas_field)+'$'

    xvar = 'zlim'
    yvara = 'budget'
    yvarb = 'nsn'

    print(res.dtype)
    legg = dict(
        zip(range(1, 4), [ADFS+'(1p)', ADFS+'(1p)+1 '+field, ADFS+'(1p)+2 '+fields]))
    leggb = dict(
        zip(range(1, 4), [ADFS+'(2p)', ADFS+'(2p)+1 '+field, ADFS+'(2p)+2 '+fields]))

    biplot(ax, res, xvar, yvara, yvarb, color, legm, legg, hatch='')

    biplot(ax, resb, xvar, yvara, yvarb, color, legm, leggb, hatch='|')

    ax[0].grid()
    ax[1].grid()

    min_budget = np.min(res['budget'])-0.001
    min_nsn = 0.99*np.min(res['nsn'])
    ax[0].set_xlim([0.5, 0.9])
    ax[0].set_ylim([min_budget, 0.15])
    ax[1].set_ylim([min_nsn, None])
    ax[1].set_xlabel('$z_{complete}$', fontsize=15)
    ax[0].set_ylabel('DD budget', fontsize=15, fontweight='bold')
    ax[1].set_ylabel('$N_{SN} (z\leq z_{complete})$', fontsize=15)
    ax[1].tick_params(axis='x', labelsize=15)
    ax[0].tick_params(axis='y', labelsize=15)
    ax[1].tick_params(axis='y', labelsize=15)

    ax[0].legend(fontsize=12, frameon=True)


def biplot(ax, res, xvar, yvara, yvarb, color, legm, legg, hatch='None'):
    """
    Function to make a biplot

    """
    for nf in np.unique(res['Nf']):
        idx = res['Nf'] == nf
        sel = res[idx]
        idxb = resb['Nf'] == nf
        selb = resb[idxb]
        legt = '{}+{}'.format(legm, legg[nf])
        nfc = nf
        if hatch == '' and nf == 3:
            nfc = 4
        if hatch != '' and nf == 1:
            nfc = 4
        ax[0].plot(sel[xvar], sel[yvara], color=color[nfc], linewidth=0.5)
        ax[0].fill_between(sel[xvar], sel[yvara], color=color[nfc],
                           label=legt, alpha=0.4, hatch=hatch)
        ax[1].plot(sel[xvar], sel[yvarb], color=color[nfc], linewidth=0.5)
        ax[1].fill_between(sel[xvar], sel[yvarb], color=color[nfc],
                           label=legt, alpha=0.4, hatch=hatch)


parser = OptionParser()

parser.add_option('--config', type=str, default='input/config.yaml',
                  help='input config file [%default]')
parser.add_option('--ultraDeep', type=str, default='COSMOS,XMM-LSS',
                  help='list of ultra deep fields[%default]')
parser.add_option('--z_ultra', type=str, default='0.9,0.9',
                  help='z complete of ultra deep fields[%default]')
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")

opts, args = parser.parse_args()

config = yaml.load(open(opts.config))

ultraDeepFields = opts.ultraDeep.split(',')
z_ultra = list(map(float, opts.z_ultra.split(',')))

zfields = dict(zip(ultraDeepFields, z_ultra))

print(config)

resa_1, resa_2 = getRes(config, zfields, opts.visitsDir,
                        'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_80.npy')

resb_1, resb_2 = getRes(config, zfields, opts.visitsDir,
                        'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_20.npy')

resa_1 = sel(resa_1)
resa_2 = sel(resa_2)
resb_1 = sel(resb_1)[::-1]
resb_2 = sel(resb_2)[::-1]

res = np.concatenate((resa_1, resb_1))
resb = np.concatenate((resa_2, resb_2))


plotFinal(config, zfields, res, resb)


plt.show()
