from scipy.interpolate import griddata
# import matplotlib.pyplot as plt
from __init__ import plt
from optparse import OptionParser
import numpy as np
# from scipy.ndimage import zoom
from scipy.ndimage.filters import gaussian_filter
from matplotlib.ticker import MaxNLocator
import pandas as pd
from wrapper import DD_Budget
"""
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 25
plt.rcParams['figure.titlesize'] = 25
plt.rcParams['legend.fontsize'] = 25
plt.rcParams['lines.linewidth'] = 3
# plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.family'] = 'Arial'
"""


class DDF_Scenario_deprecated:
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
                 H0=70., Om0=0.30,
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

        self.cadence_ref = cadence

        # get the number of visits and season length for fields with fixed zlim

        for vv in zfields.keys():
            zlim = zfields[vv]['zlimit']
            nvisits = self.visit_zlim[cadence](zlim)
            season_length = self.nvisits_seasonlength[vv](nvisits)

            zfields[vv]['Nvisits'] = nvisits.item()
            zfields[vv]['season_length'] = season_length.item(
            )/30.  # season length has to be in months!

        self.ultraFields = zfields

        if zfields:
            self.DD_list = list(set(self.DD_list)-set(list(zfields.keys())))

        self.DD_list.append('ADFS')
        print(self.DD_list, list(zfields.keys()))

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
        self.DD_list = np.unique(tab['name']).tolist()

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
                vals['Nvisits'], vals['season_length'], vals['nseasons'], self.cadence_ref)

        print('booo', Nvisits_ultraDD)
        Nvisits_DD = (budget*Nvisits_nonDD+(budget-1.)
                      * Nvisits_ultraDD)/(1.-budget)

        nvisits_night = Nvisits_DD / \
            (Nfields*Nseasons*season_length*30.)/self.cadence_ref

        print('aoooo', nvisits_night)
        zlimit = self.zlim_visit[self.cadence_ref](nvisits_night)

        return np.round(zlimit, 4)

    def budget(self, Nfields, zlim, nseasons=2, season_length=6.):

        Nvisits_nonDD = 2122176

        Nvisits = self.visit_zlim[self.cadence_ref](zlim)

        # this is not correcting the season length vs exp time
        Nvisits_DD = Nfields*self.nvisits_field(
            Nvisits, season_length, nseasons, self.cadence_ref)
        # this is correcting the season length vs exp time
        """
        res = self.get_season_length(Nfields, Nvisits, season_length)
        Nvisits_DD = Nfields*self.nvisits_field(
            Nvisits, res, nseasons, self.cadence_ref)
        """
        Nvisits_DD_ultra = 0
        for key, vals in self.ultraFields.items():
            print('io', key, vals)
            Nvisits_DD_ultra += self.nvisits_field(
                vals['Nvisits'], vals['season_length'], vals['nseasons'], self.cadence_ref)

        Nvisits_DD_tot = Nvisits_DD+Nvisits_DD_ultra
        bud = Nvisits_DD_tot/(Nvisits_DD_tot+Nvisits_nonDD)
        # print('budget', Nfields, Nvisits_DD, Nvisits_DD_ultra, bud)
        return bud

    def get_season_length(self, Nfields, Nvisits, season_length):

        tt = Nfields[0].astype(int)
        res = None
        for vv in tt:
            seas_length_all = None
            for jo in range(vv):
                dd_name = self.DD_list[jo]
                seas_length = np.array(
                    self.nvisits_seasonlength[dd_name](Nvisits[:, 0])/30.)
                seas_length = np.where(
                    seas_length > season_length, season_length, seas_length)
                if seas_length_all is None:
                    seas_length_all = np.copy(seas_length)
                else:
                    seas_length_all += np.copy(seas_length)
            seas_length_all /= vv
            seas_length_all.reshape((len(seas_length_all), 1))
            if res is None:
                res = np.copy(seas_length_all)
            else:
                res = np.column_stack((res, seas_length_all))
        return res

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
        dz = 0.0001
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
        dz = 0.0001
        survey_area = 9.6
        nsn_tot = 0

        vec_nsn = np.vectorize(self.nsn_indiv)
        nsn_field = Nfields*nseasons*vec_nsn(
            zmin, zlim, dz, season_length*30., survey_area)

        # this is to correct for the season length here
        """
        Nvisits = self.visit_zlim[self.cadence_ref](zlim)
        res = self.get_season_length(Nfields, Nvisits, season_length)

        nsn_field = Nfields*nseasons*vec_nsn(
            zmin, zlim, dz, res*30., survey_area)
        """
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

        fig, ax = plt.subplots(figsize=(15, 12))

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


def plotContourBudget_deprecated(zfields, fDir,
                                 slDir='input', slName='seasonlength_nvisits.npy',
                                 nseasons=2, season_length=6., cadence=1.,
                                 nDD_max=4, Ny=20, cosmo_res=pd.DataFrame(),
                                 toplot='sigma_w', xaxis='nddf'):

    # fig, ax = plt.subplots(nrows=2, figsize=(6, 8))
    tt = 'Deep Rolling survey'
    runtype = 'deep_rolling'
    if not zfields:
        tt = 'Deep Universal survey'
        runtype = 'universal'
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
    figb, axb = plt.subplots(figsize=(10, 8))
    figb.subplots_adjust(top=0.85)
    tot_tit = tt
    nvisitsy = 'N$_{\mathrm{visits}}^{y} \leq $'+str(Ny)
    tot_tit += '\n  {}'.format(nvisitsy)
    if zfields:
        tot_tit += ' - {}'.format(leg)
    figb.suptitle(tot_tit)
    fName = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
        Ny)
    plotContour(axb, zfields, fDir, fName,
                slDir=slDir, slName=slName, color='k',
                nseasons=nseasons, season_length=season_length,
                cadence=cadence, nDD_max=nDD_max, cosmo_res=cosmo_res,
                var=toplot, runtype=runtype, xaxis=xaxis)

    if xaxis == 'nddf':
        axb.set_xlabel('Number of deep fields$_{\mathrm{'+str(nseasons)+'}}$')
    axb.set_ylabel('$z_{\mathrm{complete}}$')
    axb.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()


def plotContour_deprecated(ax, zfields, fDir, fName,
                           slDir='input', slName='seasonlength_nvisits.npy', color='k',
                           nseasons=2, season_length=6., cadence=1., nDD_max=4,
                           cosmo_res=pd.DataFrame(), var='nsn_DD', runtype='deep_rolling', xaxis='nddf'):
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
    nseasons: int, opt
      number of seasons of observation (per field) (default: 2)
    season_length: float, opt
      season length of observation (in months) (default: 6)
    cadence: int, opt
      cadence of observations (default: 1)
    nDD_max: int, opt
      max DD number (default: 4)
    xaxis: str, opt
      xaxis variable (default: nddf)
    Returns
    ----------
    numpy arrays with results

    """

    dd_scen = DDF_Scenario(zfields, fDir, fName, cadence=cadence,
                           slDir=slDir, slName=slName)

    budmin = 0.03
    budmax = 0.15
    nfmin = 1
    nfmax = nDD_max
    zmin = 0.55
    zmax = 0.90
    if runtype == 'universal':
        zmax = 0.85

    budget = np.linspace(budmin, budmax, 1000)
    nfields = np.linspace(nfmin, nfmax, 60)
    zlim = np.linspace(zmin, zmax, 100)

    NF, ZLIMIT = np.meshgrid(nfields, zlim)
    BUD = dd_scen.budget(NF, ZLIMIT, nseasons=nseasons,
                         season_length=season_length)
    # this is the number of SN up to zlim
    # NSN = dd_scen.nsn_new(ZLIMIT, NF, nseasons=nseasons,
    #                      season_length=season_length)

    ZLIMITB, NDDF, ZVAR = getVals(
        cosmo_res, 'zcomp', xaxis, var, nbins=100, method='linear')
    ax.imshow(BUD, extent=(nfmin, nfmax, zmin, zmax),
              aspect='auto', alpha=0.25, cmap='hsv')
    # ax[1].imshow(NSN, extent=(nfmin, nfmax, zmin, zmax),
    #             aspect='auto', alpha=0.25, cmap='hsv')

    zzv = [0.5, 0.6, 0.7, 0.8, 0.9]
    if runtype == 'deep_rolling':
        zzv = [0.03, 0.05, 0.06, 0.07,  0.08, 0.10, 0.12, 0.15]
    if runtype == 'universal':
        zzv = [0.03, 0.05, 0.08, 0.12, 0.15]

    BUD = gaussian_filter(BUD, sigma=3)
    CS = ax.contour(NF, ZLIMIT, BUD, zzv, colors='k')

    fmt = {}
    strs = ['$%3.2f$' % zz for zz in zzv]
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CS.levels, strs):
        fmt[l] = s

    fontsize = 20
    ax.clabel(CS, inline=True, fontsize=fontsize,
              colors='k', fmt=fmt)

    # axb = ax.twinx()
    # zzv = [1000, 1200, 1500, 1700, 2000, 2500, 3000, 4000, 5000, 6000]
    if runtype == 'deep_rolling':
        # zzv = [1000, 1500, 2000, 2500, 3000, 3500, 4200, 4500, 5000, 6000]
        zzv = [2000, 3000, 4000, 5000, 6000]
    if runtype == 'universal':
        zzv = [4000, 5000, 6000, 8000, 10000, 14000]
    fmt = {}
    strs = ['$%1i$' % zz for zz in zzv]
    if var == 'sigma_w':
        if runtype == 'deep_rolling':
            zzv = [0.010, 0.011, 0.012, 0.013, 0.014, 0.015,
                   0.016, 0.017, 0.018, 0.020, 0.05, 0.06, 0.07, 0.08]
        if runtype == 'universal':
            zzv = [0.010, 0.013, 0.015, 0.018]
        strs = ['$%3.3f$' % zz for zz in zzv]
    # CSb = ax.contour(NDDF, ZLIMITB, gaussian_filter(
    #    ZVAR, sigma=3), zzv, colors='r')
    ZVAR = gaussian_filter(ZVAR, sigma=3)
    CSb = ax.contour(NDDF, ZLIMITB, ZVAR, zzv, colors='r', linestyles='dashed')
    print(ZVAR)
    print(ZLIMITB)
    print(NDDF)
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CSb.levels, strs):
        fmt[l] = s
    ax.clabel(CSb, inline=True, fontsize=fontsize,
              colors='r', fmt=fmt)
    ax.grid(alpha=0.3)
    ax.set_ylim([zmin, zmax])


def plotContourBudget(cosmo_res, toplot='nsn_DD', xaxis='nseasons_ultra_unique',
                      yaxis='zcomp_ultra_unique', Ny=20):

    fig, ax = plt.subplots(figsize=(12, 12))
    fig.subplots_adjust(top=0.85, left=0.15)

    tt = 'Deep Rolling survey'
    # check the type of run in cosmo_res file
    nddf_ultra = np.mean(cosmo_res['nddf_ultra'])
    leg = ''
    legb = ''
    print(cosmo_res.columns)
    # at least one dd field requested
    idx = cosmo_res['nddf_dd'] >= 1
    cosmo_res = cosmo_res[idx]

    print('hello here', nddf_ultra)
    if nddf_ultra < 1:
        tt = 'Deep Universal survey'
        runtype = 'universal'
    else:
        # at least two DD ultra fields required
        runtype = 'deep_rolling'
        idx = cosmo_res['nddf_ultra'] >= 2
        cosmo_res = cosmo_res[idx]
        ddf_ultra = np.unique(cosmo_res['ddf_ultra'])[0]
        seasons_ultra = np.unique(cosmo_res['nseasons_ultra'])[0]
        zcomp_ultra = np.unique(cosmo_res['zcomp_ultra'])[0]
        ddf_dd = np.unique(cosmo_res['ddf_dd'])[0]
        seasons_dd = np.unique(cosmo_res['nseasons_dd'])[0]
        zcomp_dd = np.unique(cosmo_res['zcomp_dd'])[0]
        print('aoooooooo', ddf_dd, 'jj', np.unique(cosmo_res['ddf_dd']))
        ddf_dd_new = {}
        for ii, vv in enumerate(ddf_dd):
            if vv not in ddf_dd_new.keys():
                ddf_dd_new[vv] = {}
                ddf_dd_new[vv]['nseasons'] = []
                ddf_dd_new[vv]['zcomp'] = []

            ddf_dd_new[vv]['nseasons'].append(seasons_dd[ii])
            ddf_dd_new[vv]['zcomp'].append(zcomp_dd[ii])

        if len(ddf_ultra) == len(seasons_ultra):
            for i in range(len(ddf_ultra)):
                leg += ddf_ultra[i] + \
                    '$_{'+str(seasons_ultra[i])+'}^{'+str(zcomp_ultra[i])+'}$'
        for key, vals in ddf_dd_new.items():
            legb += key + \
                '$_{'+str(np.sum(vals['nseasons'])) + \
                '}^{'+str(np.median(vals['zcomp']))+'}$'

    nseasons_dd = np.unique(cosmo_res['nseasons_dd_unique'].astype(int))[0]
    tot_tit = tt
    nvisitsy = 'N$_{\mathrm{visits}}^{y} \leq $'+str(Ny)
    tot_tit += '\n  {}'.format(nvisitsy)
    if xaxis == 'nddf_dd' and 'Universal' not in tt:
        tot_tit += ' - {}'.format(leg)
    else:
        tot_tit += ' - {}'.format(legb)
    fig.suptitle(tot_tit)
    plotContour(ax, cosmo_res, var=toplot,
                runtype=runtype, xaxis=xaxis, yaxis=yaxis)

    if xaxis == 'nddf_dd':
        ax.set_xlabel(
            'Number of deep fields$_{\mathrm{'+str(nseasons_dd)+'}}$')
    if xaxis == 'nseasons_ultra_unique':
        ax.set_xlabel('N$_{\mathrm{seasons}}$ per ultra-deep field')

    if yaxis == 'zcomp_dd_unique':
        ax.set_ylabel('$z_{\mathrm{complete}}^{\mathrm{deep}}$')

    if yaxis == 'zcomp_ultra_unique':
        ax.set_ylabel('$z_{\mathrm{complete}}^{\mathrm{ultra-deep}}$')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()
    return None
    # fig, ax = plt.subplots(nrows=2, figsize=(6, 8))

    runtype = 'deep_rolling'
    if not zfields:
        tt = 'Deep Universal survey'
        runtype = 'universal'
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
    figb, axb = plt.subplots(figsize=(10, 8))
    figb.subplots_adjust(top=0.85)
    tot_tit = tt
    nvisitsy = 'N$_{\mathrm{visits}}^{y} \leq $'+str(Ny)
    tot_tit += '\n  {}'.format(nvisitsy)
    if zfields:
        tot_tit += ' - {}'.format(leg)
    figb.suptitle(tot_tit)
    fName = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
        Ny)
    plotContour(axb, zfields, fDir, fName,
                slDir=slDir, slName=slName, color='k',
                nseasons=nseasons, season_length=season_length,
                cadence=cadence, nDD_max=nDD_max, cosmo_res=cosmo_res,
                var=toplot, runtype=runtype, xaxis=xaxis)

    if xaxis == 'nddf':
        axb.set_xlabel('Number of deep fields$_{\mathrm{'+str(nseasons)+'}}$')
    axb.set_ylabel('$z_{\mathrm{complete}}$')
    axb.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.show()


def plotContour(ax, cosmo_res, var='sigma_w', runtype='deep_rolling', xaxis='nddf_dd', yaxis='zcomp_ultra_unique'):
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
    nseasons: int, opt
      number of seasons of observation (per field) (default: 2)
    season_length: float, opt
      season length of observation (in months) (default: 6)
    cadence: int, opt
      cadence of observations (default: 1)
    nDD_max: int, opt
      max DD number (default: 4)
    xaxis: str, opt
      xaxis variable (default: nddf)
    Returns
    ----------
    numpy arrays with results

    """
    print(cosmo_res.columns)

    print(
        'hhhhh', cosmo_res[['nddf', 'zcomp_dd_unique', 'zcomp_ultra_unique', 'budget', 'nddf_dd']])
    budmin = 0.03
    budmax = 0.15
    nfmin = np.min(cosmo_res[xaxis])
    nfmax = np.max(cosmo_res[xaxis])
    zmin = np.min(cosmo_res[yaxis])
    zmax = np.max(cosmo_res[yaxis])
    if runtype == 'universal':
        zmax = 0.85

    budget = np.linspace(budmin, budmax, 1000)
    nfields = np.linspace(nfmin, nfmax, 60)
    zlim = np.linspace(zmin, zmax, 100)

    NF, ZLIMIT = np.meshgrid(nfields, zlim)
    # BUD = dd_scen.budget(NF, ZLIMIT, nseasons=nseasons,
    #                     season_length=season_length)
    # this is the number of SN up to zlim
    # NSN = dd_scen.nsn_new(ZLIMIT, NF, nseasons=nseasons,
    #                      season_length=season_length)

    ZLIMITB, NDDF, BUD = getVals(
        cosmo_res, yaxis, xaxis, 'budget', nbins=100, method='linear')
    print('bud', BUD)

    ax.imshow(BUD, extent=(nfmin, nfmax, zmin, zmax),
              aspect='auto', alpha=0.25, cmap='hsv')

    ZLIMITB, NDDF, ZVAR = getVals(
        cosmo_res, yaxis, xaxis, var, nbins=100, method='linear')

    # ax[1].imshow(NSN, extent=(nfmin, nfmax, zmin, zmax),
    #             aspect='auto', alpha=0.25, cmap='hsv')

    zzv = [0.5, 0.6, 0.7, 0.8, 0.9]
    if runtype == 'deep_rolling':
        zzv = [0.03, 0.05, 0.06, 0.07,  0.08, 0.10, 0.12, 0.15]
    if runtype == 'universal':
        zzv = [0.03, 0.05, 0.08, 0.12, 0.15]
    fontsize = 20

    BUD = gaussian_filter(BUD, sigma=3)
    CS = ax.contour(NDDF, ZLIMITB, BUD, zzv, colors='k')

    fmt = {}
    strs = ['$%3.2f$' % zz for zz in zzv]
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CS.levels, strs):
        fmt[l] = s

    ax.clabel(CS, inline=True, fontsize=fontsize,
              colors='k', fmt=fmt)

    # axb = ax.twinx()
    # zzv = [1000, 1200, 1500, 1700, 2000, 2500, 3000, 4000, 5000, 6000]
    if runtype == 'deep_rolling':
        # zzv = [1000, 1500, 2000, 2500, 3000, 3500, 4200, 4500, 5000, 6000]
        zzv = [2000, 2500., 3000, 3500, 4000, 4500., 5000, 6000]

    if runtype == 'universal':
        zzv = [4000, 5000, 6000, 8000, 10000, 14000]
    fmt = {}
    strs = ['$%1i$' % zz for zz in zzv]
    if var == 'sigma_w':
        if runtype == 'deep_rolling':
            zzv = [0.010, 0.011, 0.012, 0.013, 0.014, 0.015,
                   0.016, 0.017, 0.018, 0.020, 0.05, 0.06, 0.07, 0.08]
        if runtype == 'universal':
            zzv = [0.010, 0.013, 0.015, 0.018]
        strs = ['$%3.3f$' % zz for zz in zzv]
    # CSb = ax.contour(NDDF, ZLIMITB, gaussian_filter(
    #    ZVAR, sigma=3), zzv, colors='r')
    ZVAR = gaussian_filter(ZVAR, sigma=3)
    CSb = ax.contour(NDDF, ZLIMITB, ZVAR, zzv, colors='r', linestyles='dashed')
    print(ZVAR)
    print(ZLIMITB)
    print(NDDF)
    # strs = ['{}%'.format(np.int(zz)) for zz in zzvc]
    for l, s in zip(CSb.levels, strs):
        fmt[l] = s
    ax.clabel(CSb, inline=True, fontsize=fontsize,
              colors='r', fmt=fmt)
    ax.grid(alpha=0.3)
    ax.set_ylim([zmin, zmax])


def getVals(res, varx='zcomp', vary='sigma_w', varz='nddf', nbins=800, method='linear'):

    xmin, xmax = res[varx].min(), res[varx].max()
    xlim = np.linspace(xmin, xmax, nbins)
    ymin, ymax = res[vary].min(), res[vary].max()
    ylim = np.linspace(ymin, ymax, nbins)
    X, Y = np.meshgrid(xlim, ylim)
    X_grid = np.c_[np.ravel(X), np.ravel(Y)]
    # Z = griddata((res[varx],res[vary]),res[varz],(X,Y),method=method)
    print('booo', len(res[varx]), len(res[vary]), len(res[varz]))
    Z = griddata((res[varx], res[vary]), res[varz], X_grid, method=method)
    Z = Z.reshape(X.shape)
    return X, Y, Z


parser = OptionParser()

"""
parser.add_option('--ultraDeep', type=str, default='COSMOS,XMM-LSS',
                  help='list of ultra deep fields[%default]')
parser.add_option('--z_ultra', type=str, default='0.9,0.9',
                  help='z complete of ultra deep fields[%default]')
parser.add_option('--nseason_ultra', type=str, default='2,2',
                  help='number of visits for ultra deep fields[%default]')
"""
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option("--cosmoDir", type=str, default='cosmo_files',
                  help="directory where cosmo files are located[%default]")
"""
parser.add_option('--nseasons', type=int, default=2,
                  help='number of seasons (not ultra_deep fields)[%default]')
parser.add_option('--season_length', type=float, default=6.,
                  help='season length (not ultra_deep fields)[%default]')
"""
parser.add_option('--cadence', type=float, default=1.,
                  help='cadence of observation[%default]')
"""
parser.add_option('--nDD_max', type=int, default=4,
                  help='max number of DD fields[%default]')
"""
parser.add_option('--Ny', type=int, default=20,
                  help='max number of y-band visits [%default]')
parser.add_option('--var_to_plot', type='str', default='sigma_w',
                  help='var to plot in addition to the budget (sigma_w,nsn_DD) [%default]')
parser.add_option('--xaxis', type='str', default='nddf_dd',
                  help='xaxis var (nddf,nseasons_ultra) [%default]')
parser.add_option('--yaxis', type='str', default='zcomp_dd_unique',
                  help='xaxis var (nddf,nseasons_ultra) [%default]')
parser.add_option('--cosmoFile', type='str', default='cosmoSN_deep_rolling_0.80_0.80_2_2',
                  help='cosmo file to process[%default]')

opts, args = parser.parse_args()

"""
ultraDeepFields = opts.ultraDeep.split(',')
z_ultra = list(map(float, opts.z_ultra.split(',')))
nseasons_ultra = list(map(int, opts.nseason_ultra.split(',')))
"""
xaxis = opts.xaxis
visitsDir = opts.visitsDir
Ny = opts.Ny
xaxis = opts.xaxis
yaxis = opts.yaxis
var_to_plot = opts.var_to_plot

# load cosmo results
cosmoName = '{}/{}.hdf5'.format(opts.cosmoDir, opts.cosmoFile)
cosmo_res = pd.read_hdf(cosmoName)

# get the max number of y-visits@z=0.9
Ny = '{}'.format(int(np.median(cosmo_res['Ny'])))

fName = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
    Ny)

# budget class instance
budget = DD_Budget(visitsDir, fName)


# add budget colmuns
cosmo_res = budget(cosmo_res)

plotContourBudget(cosmo_res, toplot=var_to_plot,
                  xaxis=xaxis, yaxis=yaxis, Ny=opts.Ny)

"""
zfields = {}

if ultraDeepFields != ['None']:
    for i, vv in enumerate(ultraDeepFields):
        zfields[vv] = {}
        zfields[vv]['zlimit'] = z_ultra[i]
        zfields[vv]['nseasons'] = nseasons_ultra[i]

print(zfields)

# count the number of ultra_deep and get zlim
n_ultra = len(zfields.keys())
zlims = []
for kk in zfields.keys():
    print(zfields[kk])
    zlims.append(zfields[kk]['zlimit'])

# loading cosmology file
prefix = 'cosmoSN'
suffix = 'universal_10'
zli = '_'.join(['{}'.format(zz) for zz in zlims])
if n_ultra >= 2:
    zlia = '{}'.format(z_ultra[0]).ljust(4, '0')
    zlib = '{}'.format(z_ultra[1]).ljust(4, '0')
    suffix = 'deep_rolling_{}_{}_{}_{}'.format(
        zlia, zlib, nseasons_ultra[0], nseasons_ultra[1])

fName = '{}/{}_{}.hdf5'.format(opts.cosmoDir, prefix, suffix)
print('loading', fName)
cosmo_res = pd.read_hdf(fName)

print('hello zfields', zfields)
plotContourBudget(zfields, opts.visitsDir,
                  nseasons=opts.nseasons,
                  season_length=opts.season_length,
                  cadence=opts.cadence, nDD_max=opts.nDD_max,
                  Ny=opts.Ny, cosmo_res=cosmo_res, toplot=opts.var_to_plot, xaxis=xaxis)
"""
