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
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from .wrapper import Mod_z


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
        self.zminp = 0.4
        self.zmaxp = 0.95
        self.colors = dict(zip(self.bands, ['c', 'g', 'y', 'r', 'm']))
        self.Fields = ['COSMOS', 'CDFS', 'XMM-LSS', 'ELAIS', 'ADFS1', 'ADFS2']

        self.df_visits_ref = Mod_z(
            '{}/{}'.format(dir_config, file_visits_ref)).nvisits
        self.df_visits = Mod_z('{}/{}'.format(dir_config, file_visits)).nvisits
        # self.dir_config = dir_config
        # loading the configuration file for this scenario
        name = '{}/{}.yaml'.format(dir_config, configName)
        config_orig = yaml.load(open(name), Loader=yaml.FullLoader)
        #config_orig = yaml.load(open(name))

        # little modif here: convert season param to lists

        self.conf_orig = config_orig
        config = copy.deepcopy(config_orig)

        for field in config['Fields']:
            config[field]['seasons'] = self.getSeasons(
                config[field]['seasons'])

        self.process(config)
        # self.gui()

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

        zr = np.arange(0.1, 0.95, 0.05)
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
                nvisits_season = nvisits_night*30 * \
                    self.conf[fieldname]['season_length'] / \
                    self.conf[fieldname]['cadence']

                zvals = self.z[fieldname][season](nvisits_night)
                
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
        
        return df_tot

    def summary_Nvisits_single(self):
        """
        Method to estimate a summary of the budget results regarding
        the number of visits and the redshift limits.

        """

        z = np.arange(0.1, 1.0, 0.01)
        #idx = (self.budget['DD_budget'] < 0.2)
        idx = (self.budget['DD_budget'] >= 0.)
        toplot = self.budget[idx]

        """
        print('here',toplot.columns,self.zcols)
        plt.plot(toplot['zref'],toplot['DD_budget'],'ko')
        plt.show()
        """
        
        medz = toplot[self.zcols].median(axis=1)  # median redshift limit
        medbud = toplot['DD_budget']
        medvisits = toplot['Nvisits_night']

        # get zlim_min and zlim_mac cols
        r = []

        # buds = np.arange(0.000001,0.10,0.0001)

        df_bud_z = pd.DataFrame()

        for io, col in enumerate(self.zcols):
            idx = toplot[col] >= 0.25
            idx &= toplot[col] <= 0.95
            # plt.plot(toplot[idx][col],toplot[idx]['DD_budget'],color='k')
            """
            interp_ddbudget = interpolate.interp1d(toplot[idx]['DD_budget'],
                                                        toplot[idx][col],
                                                        bounds_error=False,fill_value=0)
            interp_night = interpolate.interp1d(toplot[idx]['DD_budget'],
                                                        toplot[idx]['Nvisits_night'],
                                                        bounds_error=False,fill_value=0)
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

            # df_bud_z = pd.concat([df_bud_z, df], sort=False)
            df_bud_z = pd.concat([df_bud_z, df])

            # r.append((col,interp_ddbudget(0.04),interp_ddbudget(0.08)))

        df_buz_z = df_bud_z.sort_values(by=['budget'])

        idx = df_bud_z['budget'] > 0.

        """
        print(df_bud_z[idx][['budget','z']])
        plt.plot(df_bud_z[idx]['z'],df_bud_z[idx]['budget'])
        plt.show()
        """                   
        df_min = df_bud_z[idx].groupby(['z']).min().reset_index()
        df_max = df_bud_z[idx].groupby(['z']).max().reset_index()
        df_med = df_bud_z[idx].groupby(['z']).median().reset_index()
        
        """
        colors = dict(zip(['COSMOS','CDFS','XMM-LSS','ELAIS',
                      'ADFS1','ADFS2'],['k','b','r','g','m','orange']))
        for io,col in enumerate(self.zcols):
            idx = toplot[col]>=0.3
            idx &= toplot[col]<=0.85
            fieldname = col.split('_')[1]
            plt.plot(toplot[idx][col],toplot[idx]['DD_budget'],color=colors[fieldname])

        plt.plot(df_min['z'],df_min['budget'])
        plt.plot(df_max['z'],df_max['budget'])
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
        #self.zmax = 0.90
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
        z = np.arange(zminval, 0.9, 0.05)

        self.ax1.set_ylim(ymax=0.20)

        self.ax1.set_xlim([zminval+0.01, self.zmax])
        
        self.ax1.set_ylim([self.interp_z_ddbudget(zminval), np.min(
            [100., self.interp_z_ddbudget(self.zmax)])])
        zb = np.arange(zminval, self.zmax, 0.01)
        self.ax1.fill_between(zb, self.interpmin(
            zb), self.interpmax(zb), color='yellow', alpha=0.5)
        self.ax1.plot(self.medz, self.medbud, color='k')
        self.ax1.set_ylabel(r'DD budget')
        self.ax1.set_xlabel(r'z$_{lim}$')
        self.ax1.grid()
        print('hello',self.ax1.get_ylim())
    def plotBudget_zlim_budget(self, dd_budget):
        """
        Method to plot zlimits (min, median, max) depending on the dd_budget

        Parameters
        ---------------
        dd_budget: float
          dd budget

        """
        fontsize = 15
        zlim_median, zlim_min, zlim_max = self.zlim_Nvisits_single(
            dd_budget)
        self.ax1.plot(self.ax1.get_xlim(), [dd_budget]*2, color='r', ls='--')
        alltext = ''
        for ip, val in enumerate([('min', zlim_min), ('median', zlim_median), ('max', zlim_max)]):
            self.ax1.arrow(val[1], dd_budget, 0., -dd_budget,
                           length_includes_head=True, color='b',
                           head_length=0.005, head_width=0.01)
            # self.ax1.text(0.35, 0.05-0.01*ip,
            #              '$z_{lim}^{'+val[0]+'}$='+str(np.round(val[1], 2)), fontsize=fontsize)
            alltext += '$z_{lim}^{'+val[0]+'}$='+str(np.round(val[1], 2))
            alltext += '\n '
        self.ax1.text(1.1, 0.2, alltext, fontsize=fontsize,transform=self.ax1.transAxes)
        self.ax1.text(1.1, 0.02, 'Budget: {}'.format(
            np.round(dd_budget, 3)), fontsize=fontsize,transform=self.ax1.transAxes)
        return zlim_median

    def plotNvisits(self):

        fieldName = 'COSMOS'
        season = 1

        zmin = 0.1
        zmax = 0.9
        z = np.arange(zmin, zmax, 0.01)

        self.ax2.plot(z, self.nvisits_ref[fieldName][season](
            z), color='k', label='sum')

        for b in 'grizy':
            myinterp = self.interp_ref(
                self.df_visits_ref, 'Nvisits_{}'.format(b))
            if self.runtype == 'Nvisits_single':
                self.ax2.plot(z, self.nvisits_band_ref[fieldName][season][b](
                    z), color=filtercolors[b], label='{}'.format(b))
            else:
                self.ax2.plot(z, self.nvisits_band[fieldName][season][b](
                    z), color=filtercolors[b], label='{}'.format(b))
        self.ax2.set_xlabel('z')
        self.ax2.set_ylabel('Nvisits')
        self.ax2.grid()
        self.ax2.legend()
        #self.ax2.legend(bbox_to_anchor=(1.2, -0.1),ncol=1,fontsize=12,frameon=False,loc='lower right')
        
    def plotNvisits_zlim(self, zlim=0.5):

        fieldName = 'COSMOS'
        season = 1

        fontsize = 15
        ylims = self.ax2.get_ylim()
        xlims = self.ax2.get_xlim()
        nvisits = int(np.round(self.nvisits_ref[fieldName][season](zlim)))
        yref = 0.9*ylims[1]
        scale = 0.1*ylims[1]
        self.ax2.text(1.1, 0.9, 'Nvisits={}'.format(
            nvisits), fontsize=fontsize,transform=self.ax2.transAxes)
        self.ax2.plot(xlims, [nvisits]*2, color='red', linestyle='--')
        for io, b in enumerate('grizy'):
            key = 'nvisits_{}'.format(b)
            nvisits_b = int(
                np.round(self.nvisits_band_ref[fieldName][season][b](zlim)))
            #ytext = 0.8*ylims[1]-scale*io
            ytext = 0.8-0.1*io
            self.ax2.text(1.1, ytext,
                          'Nvisits - ${}$ ={}'.format(b, nvisits_b),
                          fontsize=fontsize,
                          color=self.colors[b],transform=self.ax2.transAxes)

        zl = 'z$_{lim}$'
        self.ax2.text(zlim-0.05, 0.95*yref,
                      '{} = {}'.format(zl, np.round(zlim, 2)), fontsize=fontsize)
        self.ax2.arrow(zlim, 0.9*yref, 0., -0.9*yref,
                       length_includes_head=True, color='r',
                       head_length=5, head_width=0.01)
        self.ax2.set_ylim(0,)

    def plot_budget_zlim_deprecated(self, dd_budget=-1):
        """
        Plot to display DD budget results as a function of the redshift limit
        if dd_budget>0: an estimation of zlim (min, median, max) corresponding to dd_budget
        is displayed

        Parameters
        ----------
        dd_budget: float, opt
         DD budget (default: -1)

        """

        zminval = 0.1
        z = np.arange(zminval, 0.9, 0.05)
        if dd_budget > 0.:
            zlim_median, zlim_min, zlim_max = self.zlim_Nvisits_single(
                dd_budget)

        # Now plot

        fig1, ax1 = plt.subplots(figsize=(8, 6))
        # fig1.suptitle(self.conf['confName'])
        # ax2.set_title('{} - season {}'.format(fieldName,season))
        ax1.set_ylim(ymax=0.10)

        ax1.set_xlim([zminval+0.01, self.zmax])
        ax1.set_ylim([self.interp_z_ddbudget(zminval), np.min(
            [0.10, self.interp_z_ddbudget(self.zmax)])])
        zb = np.arange(zminval, self.zmax, 0.01)
        ax1.fill_between(zb, self.interpmin(
            zb), self.interpmax(zb), color='yellow', alpha=0.5)
        ax1.plot(self.medz, self.medbud, color='k')
        ax1.set_ylabel(r'DD budget')
        ax1.set_xlabel(r'z$_{lim}$')
        ax1.grid()

        if dd_budget > 0.:
            # draw arrows corresponding to zlim (min, max, median) values
            ax1.plot(ax1.get_xlim(), [dd_budget]*2, color='r', ls='--')
            for ip, val in enumerate([('min', zlim_min), ('max', zlim_max), ('median', zlim_median)]):
                ax1.arrow(val[1], dd_budget, 0., -dd_budget,
                          length_includes_head=True, color='b',
                          head_length=0.005, head_width=0.01)
                ax1.text(0.35, 0.04-0.005*ip,
                         '$z_{lim}^{'+val[0]+'}$='+str(np.round(val[1], 2)))

    def plot_budget_visits_deprecated(self, fieldName, season, dd_budget=-1):
        """
        Plot to display DD budget results
        The plot has two parts:
        - left side: DD budget vs zlim for the field considered
        - right side: Number of visits vs redshift limit for the field considered

        Parameters
        ----------
        dd_budget: float
        DD budget
        fieldName: str
        name of the field to display
        season:
        name of the season to display
        dd_budget: float, opt
        DD budget (default: -1)
        """

        zmin = 0.1
        zmax = 0.9
        z = np.arange(zmin, zmax, 0.05)
        if dd_budget >= 0.0000001:
            # get the number of visits for the fields

            nVisits = self.nVisits_Fields(dd_budget)
            nvisits_choice = nVisits[fieldName][season]['all']
            Nvisits_band = nVisits
        # Now plot

        interp_bud_z = interpolate.interp1d(
            self.budget['DD_budget'].values,
            self.budget['z_{}_{}'.format(fieldName, season)].values,
            bounds_error=False, fill_value=0.0)

        interp_z_bud = interpolate.interp1d(self.budget['z_{}_{}'.format(fieldName, season)].values,
                                            self.budget['DD_budget'.format(
                                                fieldName, season)].values,
                                            bounds_error=False, fill_value=0.0)

        interp_z_visits_field = interpolate.interp1d(self.budget['z_{}_{}'.format(fieldName, season)].values,
                                                     self.budget['Nvisits_night_{}_{}'.format(
                                                         fieldName, season)].values,
                                                     bounds_error=False, fill_value=0.0)
        interp_z_visits = interpolate.interp1d(self.budget['zref'],
                                               self.budget['Nvisits_night_{}_{}'.format(
                                                   fieldName, season)].values,
                                               bounds_error=False, fill_value=0.0)
        interp_visits_z = interpolate.interp1d(self.budget['Nvisits_night_{}_{}'.format(fieldName, season)].values,
                                               self.budget['zref'],
                                               bounds_error=False, fill_value=0.0)

        zb = np.arange(zmin, zmax, 0.01)
        bud = interp_z_bud(zb)
        ymax = np.min([0.1, np.max(bud)])
        zmax = zb[np.argmax(bud)]

        fig, axs = plt.subplots(
            1, 2, gridspec_kw={'hspace': 0, 'wspace': 0}, figsize=(16, 8))
        (ax1, ax2) = axs

        fig.suptitle(
            '{} - {} - season {}'.format(self.conf['confName'], fieldName, season))

        # left-handed plot: DDbudget vs zlim

        ax1.set_ylim(ymax=ymax)

        ax1.set_xlim([zmin+0.01, zmax])
        # ax1.set_ylim([interp_z_bud(zmin), np.min(
        #    [0.10, ymax])])

        budvals = interp_z_bud(zb)
        ax1.plot(zb, budvals, color='k')

        ax1.set_ylabel(r'DD budget')
        ax1.set_xlabel(r'z$_{lim}$')
        ax1.grid()

        axa = ax1.twinx()

        axa.plot(zb, budvals, color='k')

        zlims = interp_bud_z(ax1.get_ylim())
        Nvisitslim = interp_z_visits_field(zlims)

        axa.set_ylim(Nvisitslim)

        zticks = interp_bud_z(ax1.get_yticks())

        Nvisits_ticks = interp_z_visits_field(zticks)

        axa.set_yticks(Nvisits_ticks)

        if dd_budget > 0.:
            ax1.plot(ax1.get_xlim(), [dd_budget]*2, color='r', ls='--')
            zlimit = interp_bud_z(dd_budget)
            print('zlimit', zlimit, dd_budget)
            ax1.arrow(zlimit, dd_budget, 0., -dd_budget,
                      length_includes_head=True, color='b',
                      head_length=0.005, head_width=0.01)
            """
            ax1.text(0.35, 0.05,
                     '$z_{lim}$='+str(np.round(zlimit, 2)), fontsize=15)
            """
            ax1.text(1.05*zlimit, 0.5*dd_budget,
                     '$z_{lim}$='+str(np.round(zlimit, 2)), fontsize=15)
        # right-hand side plot: nvisits vs zlim
        # adjust axes
        ax2.set_ylim(axa.get_ylim())
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_yticks(axa.get_yticks())
        # make the twin to have Nvisits on the right
        ax3 = ax2.twinx()
        ax3.set_ylim(axa.get_ylim())
        ax3.set_xlim(ax1.get_xlim())
        ax3.plot(zb, interp_z_visits(zb),
                 color='k', label='all')

        for b in 'grizy':
            myinterp = self.interp_ref(
                self.df_visits_ref, 'Nvisits_{}'.format(b))
            if self.runtype == 'Nvisits_single':
                ax3.plot(z, self.nvisits_band_ref[fieldName][season][b](
                    z), color=filtercolors[b], label='{}'.format(b))
            else:
                ax3.plot(z, self.nvisits_band[fieldName][season][b](
                    z), color=filtercolors[b], label='{}'.format(b))
        ax2.grid()
        ax3.legend()
        ax2.yaxis.set_ticklabels([])
        axa.yaxis.set_ticklabels([])
        ax3.set_ylabel(r'Nvisits/night')
        ax2.set_xlabel(r'z')
        ax3.set_yticks(ax2.get_yticks())
        ax3.set_yticklabels(np.round(ax3.get_yticks()).astype(int))

        if dd_budget > 0.:

            ax2.plot(ax2.get_xlim(), [
                     nVisits[fieldName][season]['all']]*2, color='r', ls='--')
            nvisits_choice_calc = 0
            zName = 'zref'
            ax2.arrow(nVisits[fieldName][season][zName], nvisits_choice, 0.0, ax3.get_ylim()[0]-nvisits_choice,
                      length_includes_head=True, color='b',
                      head_length=1., head_width=0.01)

            for io, band in enumerate('grizy'):
                nvisits_band = int(np.round(nVisits[fieldName][season][band]))
                if nvisits_band > 0:

                    ymax = ax2.get_ylim()[1]
                    ax2.arrow(nVisits[fieldName][season][zName], nVisits[fieldName][season][band],
                              ax2.get_xlim()[
                        1]-nVisits[fieldName][season][zName], 0.,
                        length_includes_head=True, color='b',
                        head_width=1., head_length=0.01)
                    """
                    ax2.text(0.35, 0.9*nvisits_choice-0.1*io*nvisits_choice,
                             '$N_{visits}$-'+band+' = '+str(nvisits_band), fontsize=15)
                    """
                    ax2.text(0.3, 0.8*ymax-0.07*io*ymax,
                             '$N_{visits}$-'+band+' = '+str(nvisits_band), fontsize=15)

                    nvisits_choice_calc += nvisits_band

            ax2.text(0.3, 0.9*ymax,
                     '$N_{visits}$ - sum = '+str(int(nvisits_choice_calc)), fontsize=15)
            """
            ax2.text(0.35, 1.1*nvisits_choice,
                     '$N_{visits}$ - sum = '+str(int(nvisits_choice_calc)), fontsize=15)
            """

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

        self.fig.subplots_adjust(right=0.8, bottom=0.25, top=0.95)

        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)
        self.toolbar.update()

        self.plotBudget_zlim()
        self.plotNvisits()
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

        nvisits_button.grid(row=4, column=0, sticky=tk.W+tk.E)
        z_button.grid(row=4, column=1, sticky=tk.W+tk.E)
        nvisits_night_button.grid(row=5, column=0, sticky=tk.W+tk.E)
        quit_button.grid(row=5, column=1, sticky=tk.W+tk.E)

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

        # Nvisits : 2390000
        tk.Label(frame, text='Nvisits(10 yrs)', bg='white',
                 fg='black', font=font).grid(row=0)
        tk.Label(frame, text='zlim', bg='white',
                 fg='red', font=font).grid(row=1)
        tk.Label(frame, text='DD budget', bg='white',
                 fg='blue', font=font).grid(row=2)
        tk.Label(frame, text='Nvisits/night', bg='white',
                 fg='green', font=font).grid(row=3)

        entries = {}

        for vv in ['Nvisits', 'zlim', 'DDbudget', 'Nvisits_night']:
            entries[vv] = tk.Entry(frame, width=10, font=font)

        entries['Nvisits'].insert(10, 2390000)
        entries['zlim'].insert(10, '0.7')
        entries['DDbudget'].insert(10, '0.05')
        entries['Nvisits_night'].insert(10, 50)

        entries['Nvisits'].grid(row=0, column=1)
        entries['zlim'].grid(row=1, column=1)
        entries['DDbudget'].grid(row=2, column=1)
        entries['Nvisits_night'].grid(row=3, column=1)

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
        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        ddbudget = float(entries['DDbudget'].get())
        if ddbudget > 0:
            zlim = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
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
        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        zlim = float(entries['zlim'].get())
        if zlim > 0:
            # get the corresponding budget (median)
            ddbudget = self.interp_z_ddbudget(zlim)
            zlim = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
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
        # reset axes
        self.ax1.cla()
        self.ax2.cla()

        self.process(config)
        # plot references
        self.plotBudget_zlim()
        self.plotNvisits()

        # plot results of entries
        nvisits_night = float(entries['Nvisits_night'].get())
        print('nvisits/night', nvisits_night)
        if nvisits_night > 0:
            zlim = self.z_ref['COSMOS'][1](nvisits_night)
            print('zlim', zlim)
            # get the corresponding budget (median)
            ddbudget = self.interp_z_ddbudget(zlim)
            zlimb = self.plotBudget_zlim_budget(ddbudget)
            self.plotNvisits_zlim(zlim)
        # update canvas
        self.ax1.set_xlim(self.zminp, self.zmax)
        self.ax2.set_xlim(self.zminp, self.zmax)
        self.canvas.draw()
