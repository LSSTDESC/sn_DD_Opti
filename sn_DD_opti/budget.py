import yaml
import copy
from scipy import interpolate
import numpy as np
# import matplotlib.pyplot as plt
import pandas as pd
from . import filtercolors
from . import plt
import tkinter as tk
# from tkinter import tkFont
from tkinter import font as tkFont
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from .wrapper import Mod_z
from matplotlib import transforms
import os
from sn_rate import SN_Rate


class DD_Budget:
    def __init__(self, file_visits_ref, file_visits,
                 runtype='Nvisits_single', dir_config='input/sn_studies',
                 configName='DD_init'):
        """
        class estimating the DD budget vs redshift limit

        Parameters
        ---------------
        file_visits_ref: npy file
          with the following columns
          fieldname: name of the DD field
          season: season number
          Nvisits: total number of visits
          cadence: candence
          Nvisits_g,r,i,z,y: number of visits in g,r,i,z,y bands
         file_visits: npy file
          with same infos as df_visits_ref
         runtype: str, opt
           type of SNR run to consider (default: Nvisits_single)
        dir_config: str,opt
          directory where config files are located
        configName: str
         configuration (yaml) file of the DDF
        """
        self.runtype = runtype
        self.bands = 'grizy'
        self.zminp = 0.5
        self.zmaxp = 0.97
        self.colors = dict(zip(self.bands, ['c', 'g', 'y', 'r', 'm']))

        self.df_visits_ref = Mod_z(
            '{}/{}'.format(dir_config, file_visits_ref)).nvisits
        self.df_visits = Mod_z('{}/{}'.format(dir_config, file_visits)).nvisits
        # self.dir_config = dir_config
        # loading the configuration file for this scenario
        name = '{}/{}.yaml'.format(dir_config, configName)
        # config_orig = yaml.load(open(name), Loader=yaml.FullLoader)
        config_orig = yaml.load(open(name))

        # little modif here: convert season param to lists

        self.Fields = np.unique(config_orig['Fields'])
        print('hello', name, self.Fields)

        self.conf_orig = config_orig
        config = copy.deepcopy(config_orig)

        for field in config['Fields']:
            config[field]['seasons'] = self.getSeasons(
                config[field]['seasons'])

         # this is to estimate the season length as a function of the number of visits per obs night
        self.seasonLength_nvisits = self.load_SL_Nv(
            '{}/seasonlength_nvisits.npy'.format(dir_config))

        self.process(config)

        self.sn_rate = SN_Rate()
        # self.gui()

    def load_SL_Nv(self, fi):
        """
        Method to load and interpolate season length (max) vs nvisits

        Parameters
        ---------------
        fi: str
           input file

        Returns
        -----------
        dict of interpolator (key: field name)

        """

        dictout = {}
        if not os.path.exists(fi):
            print(fi, ': this file does not exist')
        else:
            tab = np.load(fi)
            for field in self.Fields:
                idx = tab['name'] == field
                sel = tab[idx]
                if len(sel) == 0:
                    print('problem here: could not find data for', field)
                else:
                    dictout[field] = interpolate.interp1d(
                        sel['nvisits'], sel['season_length'], fill_value=0.0, bounds_error=False)

        return dictout

    def process(self, config):
        """
        Method to process data to get ready for plot
        configuration loaded from config file

        Parameters
        ---------------
        configName: str
          config file name composed of the total number of visits,
        fields and associated infos (cadence, season length, seasons)

        """
        self.conf = config
        self.configName = 'from entries'

        # loading the number of visits for the case one m5 per band per season and per field

        self.nvisits_ref, self.z_ref, self.nvisits_band_ref = self.interp_visits(
            self.df_visits_ref, runtype='Nvisits_single')

        # print('test2', self.nvisits_ref['COSMOS'][1]([0.8, 0.85]))
        # loading the number of visits for the case a single  m5 per band

        self.nvisits, self.z, self.nvisits_band = self.interp_visits(
            self.df_visits, runtype='')

        # estimate the budget

        self.budget = self.budget_calc(self.runtype)

        # print('hello self budget', self.budget)
        # if self.runtype == 'Nvisits_single':
        self.summary_Nvisits_single()

    def interp_visits(self, df_tot, runtype):
        """
        Method to interpolate the number of visits vs z

        Parameters
        ---------
        df_tot: pandas df
         data used to make interpolations
        run_type: str
         type of run: Nvisits_single or Nvisits_adjusted

        Returns
        ------
        nvisits: dict of interp1d
          keys: fieldName, season; parameter: z
        z: dict of interp1d
          keys: fieldName, season; parameter: nvisits
        nvisits_band: dict of interp1d
          keys: fieldName, season, band; parameter: z

        """

        fieldnames = self.Fields
        seasons = range(1, 11)
        nvisits = {}
        nvisits_band = {}
        z = {}
        # print('hello', df_tot)
        for fieldname in fieldnames:
            nvisits[fieldname] = {}
            nvisits_band[fieldname] = {}
            z[fieldname] = {}
            for season in seasons:
                fieldname_ref = fieldname
                season_ref = season
                if runtype == 'Nvisits_single':
                    fieldname_ref = 'all'
                    season_ref = 0
                idx = df_tot['fieldname'] == fieldname_ref
                idx &= df_tot['season'] == season_ref
                idx &= np.abs(df_tot['cadence'] -
                              self.conf[fieldname]['cadence']) < 1.e-5
                sel = df_tot[idx]

                nvisits[fieldname][season] = interpolate.interp1d(
                    sel['z'], sel['Nvisits'], bounds_error=False, fill_value=0.0)
                z[fieldname][season] = interpolate.interp1d(
                    sel['Nvisits'], sel['z'], bounds_error=False, fill_value=0.0)
                nvisits_band[fieldname][season] = {}

                for b in 'grizy':
                    nvisits_band[fieldname][season][b] = interpolate.interp1d(
                        sel['z'], sel['Nvisits_{}'.format(b)], bounds_error=False, fill_value=0.0)
        # if runtype == 'Nvisits_single':
        #    print('test', nvisits['COSMOS'][1]([0.80, 0.85]))
        return nvisits, z, nvisits_band

    def interp_ref(self, df_ref, what='Nvisits'):

        idx = df_ref['fieldname'] == 'all'
        idx &= df_ref['season'] == 0
        sel = df_ref[idx]

        nvisits_ref = interpolate.interp1d(
            sel['z'], sel[what], bounds_error=False, fill_value=0.0)

        return nvisits_ref

    def budget_calc(self, runtype):
        """
        Method to estimate, vs z, the DD budget

        Parameters
        ----------
        run_type: str
         type of run: Nvisits_single or Nvisits_adjusted

        Returns
        -------
        pandas df with the following cols:
         Nvisits_fieldName_season: number of visits
                                   for all field/season considered in the scenario (conf file)
         z_fieldName_season: zlimit for all field/season considered in the scenario (conf file)
         z_ref: redshift limit corresponding to the case same number of visits per field/season/night
         Nvisits: total number of visits
         Nvisits_night: total number of visits per night
         DD_budget: DD budget

        """
        zstep = 0.01
        zr = np.arange(0.5, self.zmaxp+zstep, zstep)
        df_tot = pd.DataFrame()

        """
        if runtype == 'Nvisits_single':
            Nvisits_ref = self.nvisits_ref(zr)
        """
        self.cols = []
        self.cols_night = []
        self.zcols = []
        for fieldname in self.conf['Fields']:
            theconf = self.conf[fieldname]
            for season in theconf['seasons']:
                if runtype == 'Nvisits_single':
                    nvisits_night = self.nvisits_ref[fieldname][season](zr)
                else:
                    nvisits_night = self.nvisits[fieldname][season](zr)

                # estimate the season length here depending on the number of visits
                season_length = np.array(self.seasonLength_nvisits[fieldname](
                    nvisits_night)/30.)
                season_length[season_length >=
                              self.conf_orig[fieldname]['season_length']] = 6.

                nvisits_season = nvisits_night*30 * \
                    season_length / \
                    self.conf[fieldname]['cadence']

                zvals = self.z[fieldname][season](nvisits_night)
                # print('ici', fieldname, season,
                #      nvisits_season, zr, zvals, nvisits_night)
                zvals = zr
                zname = '{}_{}_{}'.format('z', fieldname, season)
                vname = '{}_{}_{}'.format('Nvisits', fieldname, season)
                vname_night = '{}_{}_{}'.format(
                    'Nvisits_night', fieldname, season)
                self.cols.append(vname)
                self.cols_night.append(vname_night)
                self.zcols.append(zname)
                dfa = pd.DataFrame(nvisits_season, columns=[vname])
                dfa.loc[:, zname] = zvals
                dfa.loc[:, vname_night] = nvisits_night
                dfa.loc[:, 'zref'] = zr
                if df_tot.empty:
                    df_tot = dfa.copy()
                else:
                    df_tot = df_tot.merge(
                        dfa, left_on=['zref'], right_on=['zref'])

        df_tot['Nvisits'] = df_tot[self.cols].sum(axis=1)
        df_tot['Nvisits_night'] = df_tot[self.cols_night].median(axis=1)
        df_tot['DD_budget'] = df_tot['Nvisits']/self.conf['Nvisits']

        # print('hello', df_tot[['Nvisits', 'DD_budget']])
        return df_tot

    def summary_Nvisits_single(self):
        """
        Method to estimate a summary of the budget results regarding
        the number of visits and the redshift limits.

        """
        zstep = 0.01
        z = np.arange(0.1, self.zmaxp+zstep, zstep)
        idx = (self.budget['DD_budget'] >= 0.)
        toplot = self.budget[idx]
        """
        print('here',toplot.columns,self.zcols)
        plt.plot(toplot['zref'],toplot['DD_budget'],'ko')
        plt.show()
        """

        # print('hhhh', toplot, self.zcols)
        medz = toplot[self.zcols].median(axis=1)  # median redshift limit
        medbud = toplot['DD_budget']
        medvisits = toplot['Nvisits_night']

        # get zlim_min and zlim_mac cols
        r = []

        # buds = np.arange(0.000001,0.10,0.0001)

        df_bud_z = pd.DataFrame()

        for io, col in enumerate(self.zcols):
            idx = toplot[col] >= 0.25
            idx &= toplot[col] <= self.zmaxp+0.01
            # plt.plot(toplot[idx][col],toplot[idx]['DD_budget'],color='k')
            """
            interp_ddbudget = interpolate.interp1d(toplot[idx]['DD_budget'],
                                                        toplot[idx][col],
                                                        bounds_error=False,fill_value=0)
            interp_night = interpolate.interp1d(toplot[idx]['DD_budget'],
                                                        toplot[idx]['Nvisits_night'],
                                                        bounds_error=False,fill_value=0)
            """
            # print('there we go', col, toplot[idx][col])

            """
            import matplotlib.pyplot as plt
            plt.plot(toplot[idx][col], toplot[idx][col])
            print('hhh', toplot[idx][[col, 'DD_budget']])
            plt.show()
            """
            interp_ddbudget = interpolate.interp1d(toplot[idx][col],
                                                   toplot[idx]['DD_budget'],
                                                   bounds_error=False, fill_value=0)

            interp_night = interpolate.interp1d(toplot[idx]['DD_budget'],
                                                toplot[idx]['Nvisits_night'],
                                                bounds_error=False, fill_value=0)

            df = pd.DataFrame()
            df.loc[:, 'budget'] = interp_ddbudget(z)
            df.loc[:, 'z'] = z
            df.loc[:, 'name'] = col
            df.loc[:, 'Nvisits_night'] = interp_night(df['budget'])

            # print('hhh', df[['z', 'budget']])
            # df_bud_z = pd.concat([df_bud_z, df], sort=False)
            df_bud_z = pd.concat([df_bud_z, df])

            # r.append((col,interp_ddbudget(0.04),interp_ddbudget(0.08)))

        df_buz_z = df_bud_z.sort_values(by=['budget'])

        idx = df_bud_z['budget'] >= 0.

        """
        print(df_bud_z[idx][['budget','z']])
        plt.plot(df_bud_z[idx]['z'],df_bud_z[idx]['budget'])
        plt.show()
        """
        df_min = df_bud_z[idx].groupby(['z']).min().reset_index()
        df_max = df_bud_z[idx].groupby(['z']).max().reset_index()
        df_med = df_bud_z[idx].groupby(['z']).median().reset_index()

        """
        plt.plot(df_min['z'], df_min['budget'])
        plt.plot(df_max['z'], df_max['budget'])
        plt.plot(df_med['z'], df_med['budget'])

        plt.show()
        """
        """
        colors = dict(zip(['COSMOS', 'CDFS', 'XMM-LSS', 'ELAIS',
                           'ADFS1', 'ADFS2'], ['k', 'b', 'r', 'g', 'm', 'orange']))
        for io, col in enumerate(self.zcols):
            idx = toplot[col] >= 0.3
            idx &= toplot[col] <= 1.0
            fieldname = col.split('_')[1]
            plt.plot(toplot[idx][col], toplot[idx]
                     ['DD_budget'], color=colors[fieldname])

        plt.plot(df_min['z'], df_min['budget'])
        plt.plot(df_max['z'], df_max['budget'])
        plt.show()
        """
        self.interpmin = interpolate.interp1d(
            df_min['z'], df_min['budget'], bounds_error=False, fill_value=0.00)
        self.interpmax = interpolate.interp1d(
            df_max['z'], df_max['budget'], bounds_error=False, fill_value=0.0)

        self.interpmin_ddbudget = interpolate.interp1d(
            df_max['budget'], df_max['z'], bounds_error=False, fill_value=0.10)
        self.interpmax_ddbudget = interpolate.interp1d(
            df_min['budget'], df_min['z'], bounds_error=False, fill_value=0.0)

        self.interp_ddbudget = interpolate.interp1d(
            df_med['budget'], df_med['z'], bounds_error=False, fill_value=0.0)
        self.interp_z_ddbudget = interpolate.interp1d(
            df_med['z'], df_med['budget'], bounds_error=False, fill_value=0.0)
        self.interp_z = interpolate.interp1d(
            df_med['z'], df_med['Nvisits_night'], bounds_error=False, fill_value=0.0)
        self.interp_nvisits_z = interpolate.interp1d(
            df_med['Nvisits_night'], df_med['z'], bounds_error=False, fill_value=0.0)

        """
        self.interpmin = interpolate.interp1d(
            toplot[colmin],medbud,bounds_error=False,fill_value=0.10)
        self.interpmax = interpolate.interp1d(
            toplot[colmax],medbud,bounds_error=False,fill_value=0.0)

        self.interpmin_ddbudget = interpolate.interp1d(
            medbud,toplot[colmin],bounds_error=False,fill_value=0.10)
        self.interpmax_ddbudget = interpolate.interp1d(
            medbud,toplot[colmax],bounds_error=False,fill_value=0.0)

        self.interp_ddbudget = interpolate.interp1d(
            medbud,medz,bounds_error=False,fill_value=0.0)
        self.interp_z =  interpolate.interp1d(
            medz,medvisits,bounds_error=False,fill_value=0.0)
        """
        self.medbud = df_med['budget'].values
        self.medz = df_med['z'].values
        self.medvisits = df_med['Nvisits_night'].values
        self.zmax = df_max['z'].max()
        self.zmax = df_bud_z[idx]['z'].max()

    def zlim_Nvisits_single(self, dd_value):
        """
        Method to estimate some results corresponding to a given DD budget

        Parameters
        ----------
        dd_value: float
         DD budget

        Returns:
        -------
        zlim_median: float
          median redshift limit
        zlim_min: float
          min redshift limit
        zlim_max: float
          max redshift limit
        nvisits_choice: float
          number of visits (per night) corresponding to zlim_median
        Nvisits_band: dict
          number of visits per night and per band (key)
        """

        if self.runtype == 'Nvisits_single':
            zlim_median = self.interp_ddbudget(dd_value)
            zlim_min = self.interpmin_ddbudget(dd_value)
            zlim_max = self.interpmax_ddbudget(dd_value)
        else:
            nVisits = self.nVisits_Fields(dd_value)
            for key, vals in nVisits.items():
                for keyb, valb in vals.items():
                    zlim_median = np.median(valb['zref'])
                    zlim_min = np.min(valb['zref'])
                    zlim_max = np.max(valb['zref'])
                    # self.zmax = zlim_max

        """
        nvisits_choice = self.interp_z(zlim_median)
        Nvisits_band = {}
        for b in 'rizy':
            myinterp = self.interp_ref(
                self.df_visits_ref, 'Nvisits_{}'.format(b))
            Nvisits_band[b] = myinterp(zlim_median)

        return zlim_median, zlim_min, zlim_max, nvisits_choice, Nvisits_band
        """
        return zlim_median, zlim_min, zlim_max

    def nVisits_Fields(self, dd_value):
        """
        Method to estimate the number of visits per fields
        depending on the dd value

        Parameters
        ---------------
        dd_value: float
          DD budget

        Returns
        ----------
        nVisits: dict
          dict with the number of visits (per obs night)

        """
        nVisits = {}
        # print(self.budget.columns)
        for fieldName in self.conf['Fields']:
            nVisits[fieldName] = {}
            theconf = self.conf[fieldName]
            for seas in theconf['seasons']:
                nVisits[fieldName][seas] = {}
                zName = 'z_{}_{}'.format(fieldName, seas)
                nvisitsName = 'Nvisits_{}_{}'.format(fieldName, seas)
                myinterp = interpolate.interp1d(
                    self.budget['DD_budget'].values, self.budget['Nvisits_night_{}_{}'.format(fieldName, seas)].values)
                nVisits[fieldName][seas]['all'] = np.asscalar(
                    myinterp(dd_value))
                myinterpz = interpolate.interp1d(
                    self.budget['DD_budget'].values, self.budget[zName].values)
                nVisits[fieldName][seas]['zlim'] = np.asscalar(
                    myinterpz(dd_value))

                myinterpz = interpolate.interp1d(
                    self.budget['DD_budget'].values, self.budget['zref'].values)
                nVisits[fieldName][seas]['zref'] = np.asscalar(
                    myinterpz(dd_value))

                for b in 'grizy':
                    if self.runtype == 'Nvisits_single':
                        nVisits[fieldName][seas][b] = np.asscalar(self.nvisits_band_ref[fieldName][seas][b](
                            myinterpz(dd_value)))
                    else:
                        nVisits[fieldName][seas][b] = np.asscalar(self.nvisits_band[fieldName][seas][b](
                            myinterpz(dd_value)))

        return nVisits

    def plotBudget_zlim(self):
        """
        Plot to display DD budget results as a function of the redshift limit

        """
        zminval = 0.3
        z = np.arange(zminval, self.zmax, 0.05)

        self.ax1.set_ylim(ymax=0.50)

        self.ax1.set_xlim([zminval+0.01, self.zmax])

        # self.ax1.set_ylim([self.interp_z_ddbudget(zminval), np.min(
        #    [100., self.interp_z_ddbudget(self.zmax)])])
        zb = np.arange(zminval, self.zmax, 0.01)
        self.ax1.fill_between(zb, self.interpmin(
            zb), self.interpmax(zb), color='yellow', alpha=0.5)
        self.ax1.plot(self.medz, self.medbud, color='k')
        self.ax1.set_ylabel(r'DD budget')
        self.ax1.set_xlabel(r'z$_{lim}$')
        self.ax1.grid()
        ymax = self.interpmax(self.zmax)
        ymax = np.max(self.interpmax(zb))
        print('ooo', ymax, self.zmax)
        self.ax1.set_ylim(ymax=1.1*ymax)
        self.ax1.set_xlim(xmax=self.zmax)

    def plotBudget_zlim_budget(self, dd_budget):
        """
        Method to plot zlimits (min, median, max) depending on the dd_budget

        Parameters
        ---------------
        dd_budget: float
          dd budget

        """
        fontsize = 12
        zlim_median, zlim_min, zlim_max = self.zlim_Nvisits_single(
            dd_budget)
        self.ax1.plot(self.ax1.get_xlim(), [
                      dd_budget]*2, color='r', ls='--')
        alltext = ''
        n_SN = int(self.nSN(self.conf, zlim_median))
        # for ip, val in enumerate([('min', zlim_min), ('median', zlim_median), ('max', zlim_max)]):
        for ip, val in enumerate([('$z_{complete}^{median}$', zlim_median), ('$N_{SN}$', n_SN)]):
            self.ax1.arrow(val[1], dd_budget, 0., -dd_budget,
                           length_includes_head=True, color='b',
                           head_length=0.005, head_width=0.01)
            # self.ax1.text(0.35, 0.05-0.01*ip,
            #              '$z_{lim}^{'+val[0]+'}$='+str(np.round(val[1], 2)), fontsize=fontsize)
            tt = 'lim'
            tt = 'complete'
            #alltext += '$z_{'+tt+'}^{'+val[0]+'}$='+str(np.round(val[1], 2))
            alltext += '{}={}'.format(val[0], str(np.round(val[1], 2)))
            alltext += '\n \n'
        self.ax1.text(1.2, 0.5, alltext, fontsize=fontsize,
                      transform=self.ax1.transAxes, color='b', fontweight='bold')
        self.ax1.text(1.2, 0.2, 'Budget: {}'.format(
            np.round(dd_budget, 3)), fontsize=fontsize, transform=self.ax1.transAxes, color='r', fontweight='bold')
        return zlim_median

    def plotNvisits(self):
        """
        Method to plot Nvisits vs redshift for median m5 (references)

        """
        season = 1

        zstep = 0.01
        zmin = 0.1
        zmax = self.zmax
        z = np.arange(zmin, zmax+zstep, zstep)

        lstyle = dict(
            zip([1, 2, 3, 4], ['solid', 'dotted', 'dashed', 'dashdot']))

        io = -1
        r = []
        rb = []
        cad = []

        for fieldName, cadence in self.fields_cad.items():
            io += 1
            lsa = lstyle[io+1]
            lla, = self.ax2.plot(z, self.nvisits_ref[fieldName][season](
                z), color='k', label='cadence={}'.format(cadence), ls=lsa)
            rb.append([lla])
            days = 'days$^{-1}$'
            cad.append('cadence: {} {}'.format(cadence, days))
            # self.ax2.plot(z, self.nvisits_ref[fieldName][season](z), color='k', label='sum')

            for b in 'grizy':
                color = filtercolors[b]
                label = '{}'.format(b)
                myinterp = self.interp_ref(
                    self.df_visits_ref, 'Nvisits_{}'.format(b))
                if self.runtype == 'Nvisits_single':
                    func = self.nvisits_band_ref[fieldName][season][b]
                else:
                    func = self.nvisits_band[fieldName][season][b]

                if io == 0:
                    la, = self.ax2.plot(
                        z, func(z), color=color, label=label, ls=lsa)
                    r.append([la])
                else:
                    self.ax2.plot(z, func(z), color=color, ls=lsa)

        leg1 = self.ax2.legend([l[0] for l in r], 'grizy',
                               bbox_to_anchor=(-0.06, -0.01), frameon=False)
        self.ax2.set_xlabel('z')
        self.ax2.set_ylabel('Nvisits')
        self.ax2.grid()
        # self.ax2.legend()
        self.ax2.legend([l[0] for l in rb], cad, loc=2, frameon=False)
        self.ax2.add_artist(leg1)
        self.ax2.set_ylim(ymax=self.conf['Nvisits_night'])
        self.ax2.set_xlim(xmax=self.zmaxp)

    def plotNvisits_zlim(self, zlim=0.5):
        """
        Method to estimate Nvisits for zlim

        Parameters
        ---------------
        zlim: float, opt
         redshift limit (default: 0.5)
        """

        # fieldName = 'COSMOS'
        season = 1

        fontsize = 12
        ylims = self.ax2.get_ylim()
        xlims = self.ax2.get_xlim()

        yref = 0.9*ylims[1]
        scale = 0.1*ylims[1]

        it = -1

        colors = 'rk{}'.format(''.join([filtercolors[b] for b in 'grizy']))

        for fieldName in self.conf['Fields']:
            it += 1
            words = []
            words.append(fieldName.ljust(7))
            func = self.nvisits_ref[fieldName][season]
            nvisits = int(np.round(func(zlim)))
            seas_length = self.seasonLength_nvisits[fieldName](
                nvisits)/30.
            self.conf[fieldName]['season_length'] = np.round(
                np.min([seas_length, self.conf_orig[fieldName]['season_length']]), 1)

            self.ax2.plot(xlims, [nvisits]*2, color='red',
                          linestyle='solid', linewidth=0.1)
            words.append('{}'.format(nvisits))

            for io, b in enumerate('grizy'):
                key = 'nvisits_{}'.format(b)
                func = self.nvisits_band_ref[fieldName][season][b]
                nvisits_b = int(np.round(func(zlim)))
                words.append('{}'.format(nvisits_b))

            ytext = 1.0-0.1*it
            self.rainbow_text(1.12, ytext, words, colors,
                              ax=self.ax2, fontsize=fontsize, fontweight='bold')

        ttxt = 'lim'
        ttxt = 'complete'
        zl = 'z$_{'+ttxt+'}$'
        self.ax2.text(zlim-0.05, 0.95*yref,
                      '{} = {}'.format(zl, np.round(zlim, 2)), fontsize=fontsize)
        self.ax2.arrow(zlim, 0.9*yref, 0., -0.9*yref,
                       length_includes_head=True, color='r',
                       head_length=5, head_width=0.01)
        self.ax2.set_ylim(0,)

    def printVisits(self, dd_value):
        """
        Method to print(in a table) the number of visits and zlim
        per field and per night corresponding to a DD budget

        Parameters
        ----------
        dd_value: float
         DD budget

        """

        # get infos corresponding to this budget

        nVisits = self.nVisits_Fields(dd_value)

        nameConv = dict(zip(['season', 'all', 'zlim', 'zref', 'r', 'i', 'z', 'y'],
                            ['season', 'Nvisits', 'zlim', 'zref',
                             'Nvisits_r', 'Nvisits_i',
                             'Nvisits_z', 'Nvisits_y']))

        for key, vals in nVisits.items():
            fig, ax = plt.subplots()
            names = []
            vv = []
            names.append('season')
            io = -1
            for keyb, valb in vals.items():
                io += 1
                ro = [keyb]
                for keyc, valc in valb.items():
                    if keyc != 'zlim' and keyc != 'zref':
                        ro.append(int(np.round(valc)))
                    else:
                        ro.append(np.round(valc, 2))
                    if io == 0:
                        names.append(nameConv[keyc])
                vv.append(ro)

            tab = np.rec.fromrecords(vv, names=names)
            ll = '{} - {} - DD budget: {}'.format(
                self.conf['confName'], key, dd_value)
            ax.text(0.2, 0.8, ll)
            ax.table(cellText=tab, colLabels=names, loc='center')

            plt.axis('off')

    def dict_from_fields(self, fields):
        """
        Method to transform a set of tk widgets to dict

        Parameters
        ----------------
        fields: dict
          dict of tk widgets

        """

        resdict = {}
        fieldlist = []
        for keya, vala in fields.items():
            for keyb, valb in vala.items():
                if keyb == 'state':
                    if valb.get() == 1:
                        resdict[keya] = {}
                        fieldlist.append(keya)
                    else:
                        break
                else:
                    if keyb != 'check':
                        data = valb.get()
                        if keyb != 'seasons':
                            resdict[keya][keyb] = float(data)
                        else:

                            resdict[keya][keyb] = self.getSeasons(data)

        resdict['Fields'] = fieldlist

        return resdict

    def getSeasons(self, data):
        """
        Method to get the number of seasons from GUI

        Parameters
        ---------------
        data: str
          data to process

        Returns
        -----------
        list of seasons to consider

        """
        li = None
        if '-' in data:
            spl = data.split('-')
            li = list(range(int(spl[0]), int(spl[1])+1))
        else:
            spl = data.split(',')
            li = list(map(int, spl))

        return li

    def cadence_fields(self):
        """
        Method to estimate the list of fields with the same cadence
        and to select one field per cadence

        """
        # get the cadences involved
        cad_fields = {}
        for field in self.conf['Fields']:
            cadence = int(self.conf[field]['cadence'])
            if cadence not in cad_fields.keys():
                cad_fields[cadence] = []
            cad_fields[cadence].append(field)

        fields_cad = {}
        for key, vals in cad_fields.items():
            fields_cad[vals[0]] = key

        self.cad_fields = cad_fields
        self.fields_cad = fields_cad

    def popupmsg(self, msg):
        LARGE_FONT = ("Verdana", 12)
        NORM_FONT = ("Verdana", 10)
        SMALL_FONT = ("Verdana", 8)
        popup = tk.Tk()
        popup.wm_title("!")
        label = ttk.Label(popup, text=msg, font=NORM_FONT)
        label.pack(side="top", fill="x", pady=10)
        B1 = ttk.Button(popup, text="Okay", command=popup.destroy)
        B1.pack()
        popup.mainloop()

    def rainbow_text(self, x, y, strings, colors, orientation='horizontal',
                     ax=None, **kwargs):
        """
        Take a list of *strings* and *colors* and place them next to each
        other, with text strings[i] being shown in colors[i].

        Parameters
        ----------
        x, y : float
          Text position in data coordinates.
        strings : list of str
          The strings to draw.
        colors : list of color
          The colors to use.
        orientation : {'horizontal', 'vertical'}
        ax : Axes, optional
        The Axes to draw into. If None, the current axes will be used.
        **kwargs
        All other keyword arguments are passed to plt.text(), so you can
        set the font size, family, etc.
        """
        if ax is None:
            ax = plt.gca()
        t = ax.transAxes
        canvas = ax.figure.canvas

        assert orientation in ['horizontal', 'vertical']
        if orientation == 'vertical':
            kwargs.update(rotation=90, verticalalignment='bottom')

        for s, c in zip(strings, colors):
            text = ax.text(x, y, s + " ", color=c, transform=t, **kwargs)

            # Need to draw to update the text position.
            text.draw(canvas.get_renderer())
            ex = text.get_window_extent()
            if orientation == 'horizontal':
                t = transforms.offset_copy(
                    text.get_transform(), x=ex.width, units='dots')
            else:
                t = transforms.offset_copy(
                    text.get_transform(), y=ex.height, units='dots')


class GUI_Budget(DD_Budget):
    """
    class building a GUI to show
    the DD budget vs redshift limit
    inherits from DD_Budget

    Parameters
    ---------------
    file_visits_ref: npy file
     with the following columns
          fieldname: name of the DD field
          season: season number
          Nvisits: total number of visits
          cadence: candence
          Nvisits_g,r,i,z,y: number of visits in g,r,i,z,y bands
         file_visits: npy file
          with same infos as df_visits_ref
         runtype: str, opt
           type of SNR run to consider (default: Nvisits_single)
        dir_config: str,opt
          directory where config files are located
        configName: str
         configuration (yaml) file of the DDF

    """

    def __init__(self, file_visits_ref, file_visits,
                 runtype='Nvisits_single', dir_config='input/sn_studies',
                 configName='DD_init'):
        super().__init__(file_visits_ref=file_visits_ref, file_visits=file_visits,
                         runtype=runtype, dir_config=dir_config,
                         configName=configName)

        # build the GUI
        root = tk.Tk()
        root.title('DD optimisation for SN')

        # figure where results are displayed
        self.fig = plt.Figure(figsize=(15, 9), dpi=100)
        # self.fig.suptitle(self.configName.split('.')[0], fontsize=15)
        gs = self.fig.add_gridspec(2, 1)
        self.ax1 = self.fig.add_subplot(gs[0, 0])
        self.ax2 = self.fig.add_subplot(gs[1, 0])

        self.fig.subplots_adjust(right=0.7, bottom=0.25, top=0.95)

        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)
        self.toolbar.update()

        self.plotBudget_zlim()
        self.cadence_fields()
        self.plotNvisits()
        self.zmax = 0.97
        self.ax1.set_xlim(self.zminp, self.zmax)
        self.ax2.set_xlim(self.zminp, self.zmax)
        # common font
        helv36 = tkFont.Font(family='Helvetica', size=15, weight='bold')

        # building the GUI
        # frames
        button_frame = tk.Frame(master=root, bg="white")
        button_frame.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        # button_frame.place(relx=.9, rely=.5, anchor="c")
        button_frame.place(relx=.9, rely=.8, anchor="c")

        fields_frame = tk.Frame(master=root, bg="white")
        fields_frame.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        # button_frame.place(relx=.9, rely=.5, anchor="c")
        fields_frame.place(relx=.4, rely=.85, anchor="c")

        # entries
        ents = self.make_entries(button_frame, font=helv36)
        fields = self.make_fields(fields_frame, font=helv36)
        self.fields_tab = fields
        # buttons
        heightb = 3
        widthb = 6

        nvisits_button = tk.Button(
            button_frame, text="from budget",
            command=(lambda e=ents, b=fields: self.updateData(e, b)),
            bg='yellow', height=heightb, width=widthb, fg='blue', font=helv36)

        z_button = tk.Button(
            button_frame, text="from zlim",
            command=(lambda e=ents, b=fields: self.updateData_z(e, b)), bg='yellow',
            height=heightb, width=widthb, fg='red', font=helv36)

        nvisits_night_button = tk.Button(
            button_frame, text="from \n Nvisits/night",
            command=(lambda e=ents, b=fields: self.updateData_nvisits_night(e, b)),
            bg='yellow', height=heightb, width=widthb, fg='green', font=helv36)

        quit_button = tk.Button(button_frame, text="Quit",
                                command=root.quit, bg='yellow',
                                height=heightb, width=widthb, fg='black', font=helv36)

        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        button_frame.columnconfigure(2, weight=2)

        nvisits_button.grid(row=5, column=0, sticky=tk.W+tk.E)
        z_button.grid(row=5, column=1, sticky=tk.W+tk.E)
        nvisits_night_button.grid(row=6, column=0, sticky=tk.W+tk.E)
        quit_button.grid(row=6, column=1, sticky=tk.W+tk.E)

        # quit_button.pack(side=tk.BOTTOM)
        root.mainloop()

    def make_entries(self, frame, font):
        """
        Method to make entries available to the GUI

        Parameters
        ---------------
        frame: tk.Frame
          frame where entries will be located
        font: tkFont
          font used for the label

        Returns
        ----------
        entries: dict
          dict of entries

        """

        tk.Label(frame, text='Nvisits/night(max)', bg='white',
                 fg='black', font=font).grid(row=0)
        tk.Label(frame, text='Nvisits(10 yrs)', bg='white',
                 fg='black', font=font).grid(row=1)
        tk.Label(frame, text='zcomp', bg='white',
                 fg='red', font=font).grid(row=2)
        tk.Label(frame, text='DD budget', bg='white',
                 fg='blue', font=font).grid(row=3)
        tk.Label(frame, text='Nvisits/night', bg='white',
                 fg='green', font=font).grid(row=4)

        entries = {}

        for vv in ['Nvisits', 'zlim', 'DDbudget', 'Nvisits_night', 'Nvisits_night_max']:
            entries[vv] = tk.Entry(frame, width=10, font=font)

        entries['Nvisits'].insert(10, 2390000)
        entries['Nvisits_night_max'].insert(10, 300)
        entries['zlim'].insert(10, '0.7')
        entries['DDbudget'].insert(10, '0.05')
        entries['Nvisits_night'].insert(10, 50)

        entries['Nvisits_night_max'].grid(row=0, column=1)
        entries['Nvisits'].grid(row=1, column=1)
        entries['zlim'].grid(row=2, column=1)
        entries['DDbudget'].grid(row=3, column=1)
        entries['Nvisits_night'].grid(row=4, column=1)

        return entries

    def make_fields(self, frame, font):

        # frame.grid_configure(minsize=150)
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(20, weight=1)
        fields = {}
        tk.Label(frame, text='cadence', bg='white',
                 fg='black', font=font).grid(row=1, column=0, sticky=tk.W)
        tk.Label(frame, text='season length', bg='white',
                 fg='black', font=font).grid(row=2, column=0, sticky=tk.W)
        tk.Label(frame, text='seasons', bg='white',
                 fg='black', font=font).grid(row=3, column=0, sticky=tk.W)
        for icol, fieldname in enumerate(self.Fields):
            fields[fieldname] = {}
            fields[fieldname]['state'] = tk.IntVar()
            fields[fieldname]['check'] = tk.Checkbutton(
                frame, text=fieldname, font=font, bg='white',
                highlightthickness=0, bd=0, variable=fields[fieldname]['state'])
            fields[fieldname]['check'].select()
            fields[fieldname]['check'].grid(
                row=0, column=icol+1, sticky=tk.W, padx=0)

            fields[fieldname]['cadence'] = tk.Entry(frame, width=10, font=font)
            fields[fieldname]['season_length'] = tk.Entry(
                frame, width=10, font=font)
            fields[fieldname]['seasons'] = tk.Entry(frame, width=10, font=font)

            for ib, val in enumerate(['cadence', 'season_length', 'seasons']):
                fields[fieldname][val].insert(
                    10, str(self.conf_orig[fieldname][val]))
                fields[fieldname][val].grid(row=1+ib, column=icol+1)

            # fields[fieldname]['season_length'].insert(10, '180')
            # fields[fieldname]['seasons'].insert(10, '1,2')
            """
            fields[fieldname]['cadence'].grid(row=1+id, column=icol+1)
            fields[fieldname]['season_length'].grid(row=2, column=icol+1)
            fields[fieldname]['seasons'].grid(row=3, column=icol+1)
            """
        return fields

    def updateData(self, entries, fields):
        """
        Method to update the plot when buttons are clicked

        Parameters
        ---------------
        entries: dict
           dict of entries

        """

        config = self.dict_from_fields(fields)

        self.Fields = config['Fields']
        config['Nvisits'] = int(entries['Nvisits'].get())
        config['Nvisits_night'] = int(entries['Nvisits_night_max'].get())

        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.cadence_fields()
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        ddbudget = float(entries['DDbudget'].get())
        if ddbudget > 0:
            zlim = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
            # update season length here
            self.update_tab(config)
        # update canvas
        self.ax1.set_xlim(self.zminp, self.zmax)
        self.ax2.set_xlim(self.zminp, self.zmax)
        self.canvas.draw()

    def updateData_z(self, entries, fields):
        """
        Method to update the plot when buttons are clicked

        Parameters
        ---------------
        entries: dict
           dict of entries

        """

        config = self.dict_from_fields(fields)

        self.Fields = config['Fields']
        config['Nvisits'] = int(entries['Nvisits'].get())
        config['Nvisits_night'] = int(entries['Nvisits_night_max'].get())

        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.cadence_fields()
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        zlim = float(entries['zlim'].get())
        if zlim > 0:
            # get the corresponding budget (median)
            ddbudget = self.interp_z_ddbudget(zlim)
            zlim = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
            # update season_length here
            self.update_tab(config)
        # update canvas
        self.ax1.set_xlim(self.zminp, self.zmax)
        self.ax2.set_xlim(self.zminp, self.zmax)
        self.canvas.draw()

    def updateData_nvisits_night(self, entries, fields):
        """
        Method to update the plot when buttons are clicked

        Parameters
        ---------------
        entries: dict
           dict of entries

        """

        config = self.dict_from_fields(fields)

        self.Fields = config['Fields']
        config['Nvisits'] = int(entries['Nvisits'].get())
        config['Nvisits_night'] = int(entries['Nvisits_night_max'].get())

        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        nvisits_night = float(entries['Nvisits_night'].get())
        # print('nvisits/night', nvisits_night)
        # zlim = {}
        if nvisits_night > 0:
            # for field in self.conf['Fields']:
            #    zlim[field] = self.z_ref[field][1](nvisits_night)
            # print('zlim', zlim)
            zlim = self.z_ref['COSMOS'][1](nvisits_night)
            # get the corresponding budget (median)
            ddbudget = self.interp_z_ddbudget(zlim)
            zlimb = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
            # update season length here
            self.update_tab(config)
        # update canvas
        self.ax1.set_xlim(self.zminp, self.zmax)
        self.ax2.set_xlim(self.zminp, self.zmax)
        self.canvas.draw()

    def update_tab(self, config, which='season_length'):
        """
        Method to update fields of the table according to config values

        Parameters
        ---------------
        config: dict
          config dict to use to update
        which: str, opt
          field to update (default: season_length)

        """
        for fieldname in config['Fields']:
            self.fields_tab[fieldname][which].delete(0, "end")
            self.fields_tab[fieldname][which].insert(
                10, str(config[fieldname][which]))

    def nSN(self, config, zlim):
        """
        Method to estimate the total number of SN up to z<zlim

        Parameters
        --------------
        config: dict
           configuration file of the fields
        zlim: float
          redshift limit for the calculation

        Returns
        -----------
        the total number of SN(z<zlim)

        """

        survey_area = 9.6  # 1 LSST foV
        zmin = 0.01
        account_for_edges = True

        n_SN = 0
        for fieldName in self.Fields:
            seasons = config[fieldName]['seasons']
            for seas in seasons:
                season_length = config[fieldName]['season_length']*30.
                zz, rate, err_rate, nsn, err_nsn = self.sn_rate(zmin=zmin,
                                                                zmax=zlim,
                                                                survey_area=survey_area,
                                                                account_for_edges=account_for_edges,
                                                                duration=season_length)
                n_SN += np.cumsum(nsn)[-1]

        return n_SN
