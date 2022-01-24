import pandas as pd
import numpy as np
from scipy.interpolate import interp1d


class Mod_z:
    """
    Method to modify values of input file

    """

    def __init__(self, fName):

        # load the file

        tab = np.load(fName, allow_pickle=True)

        tabdf = pd.DataFrame.from_records(tab)

        tabdf = tabdf.replace('ADFS1', 'ADFS')

        tabmod = tabdf.groupby(['cadence']).apply(
            lambda x: self.mod(x))

        #tabmod = tabdf
        ll = []
        for b in 'grizy':
            ll.append('Nvisits_{}'.format(b))
        tabmod['Nvisits'] = tabmod[ll].sum(axis=1)

        self.nvisits = tabmod

    def mod(self, grp):
        """
        Method to modify a group


        Parameters
        --------------
        grp : pandas df group

        Returns
        -----------
        modified grp group
        """
        for band in 'r':
            what = 'Nvisits_{}'.format(band)
            idx = grp[what] > 1.e-21
            idx &= grp[what] < 1.
            grp.loc[idx, what] = 1.

            if band == 'g' or band == 'r':
                Index_label = grp[grp[what] < 1.e-10].index.tolist()
                Index_label_p = grp[grp[what] > 1.e-10].index.tolist()
                #print('index', Index_label, Index_label_p)
                grp.loc[Index_label, what] = grp.loc[Index_label_p[-1]][what]
            #print('io', grp[what])

        grp['Nvisits_g'] = 2
        #grp['Nvisits_r'] = 4
        return grp


class DD_Budget:
    """
    class to estimate the DD budget for a survey configuration

    Parameters
    ---------------
    fDir: str
      location dir of (nvisits, zlim) file
    fName: str
      (nvisits, zlim) filename
    cadence: int, opt
      cadence of observation (default: 1 day)
    Nvisits_nonDD: int,opt
      total number of non-DD visits (default: 2122176)
    slDir: str, opt
      location dir of (seasonlength vs nvisits) file (default: input)
    slName: str, opt
      (seasonlength vs nvisits) (default: seasonlength_nvisits.py)
    """

    def __init__(self, fDir, fName,
                 cadence=1, Nvisits_nonDD=2122176,
                 slDir='input', slName='seasonlength_nvisits.npy'):

        self.cadence = cadence
        self.Nvisits_nonDD = Nvisits_nonDD
        self.zlim_visit, self.visit_zlim = self.load(fDir, fName)

        self.nvisits_seasonlength = self.load_seasonlength_nvisits(
            slDir, slName)

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

        # tab = np.load(fName, allow_pickle=True)
        tab = Mod_z('{}/{}'.format(fDir, fName)).nvisits
        res = {}
        resb = {}
        # print(tab)
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

    def __call__(self, cosmodf):

        # print(cosmodf.columns)
        cosmodf['budget'] = cosmodf.apply(
            lambda x: self.budget(x), axis=1)

        cosmodf['zcomp_dd_unique'] = cosmodf.apply(
            lambda x: np.mean(x['zcomp_dd']), axis=1)
        cosmodf['zcomp_ultra_unique'] = cosmodf.apply(
            lambda x: np.mean(x['zcomp_ultra']), axis=1)
        cosmodf['nseasons_ultra_unique'] = cosmodf.apply(
            lambda x: np.mean(x['nseasons_ultra']), axis=1)
        cosmodf['nseasons_dd_unique'] = cosmodf.apply(
            lambda x: np.mean(x['nseasons_dd']), axis=1)
        cosmodf['nddf_ultra'] = cosmodf.apply(
            lambda x:  len(x['ddf_ultra']), axis=1)
        cosmodf['nddf_dd'] = cosmodf.apply(
            lambda x:  len(x['ddf_dd']), axis=1)
        cosmodf['nddf'] = cosmodf['nddf_ultra']+cosmodf['nddf_dd']

        #print('finally', cosmodf)

        return cosmodf

    def budget(self, grp):
        """
        Method to estimate the budget corresponding to a survey summarized in grp

        Parameters
        ---------------
        grp: pandas group
          with infos related to the survey
        """

        Nvisits_tot = 0
        Nvisits_tot += self.Nvisits_field(grp, 'ultra')
        Nvisits_tot += self.Nvisits_field(grp, 'dd')

        budget = Nvisits_tot/(Nvisits_tot+self.Nvisits_nonDD)

        # return pd.DataFrame({'budget': [budget]})
        return budget

    def Nvisits_field(self, grp, suffix='ultra'):

        Nvisits_tot = 0
        ddf = grp['ddf_{}'.format(suffix)]
        zcomp = grp['zcomp_{}'.format(suffix)]
        nseasons = grp['nseasons_{}'.format(suffix)]
        nddf = len(ddf)
        #print('helli', ddf, zcomp)
        for i in range(nddf):
            # get the number of visits per night corresponding to zcomp
            nvisits_night = self.visit_zlim[self.cadence](zcomp[i])
            season_length = self.nvisits_seasonlength[ddf[i]](
                nvisits_night)
            season_length = np.min([season_length, 180.])
            Nvisits_all_seas = self.Nvisits_seasons(
                nvisits_night, season_length, nseasons[i])
            # print(ddf[i], zcomp[i],
            #      nseasons[i], nvisits_night, season_length, Nvisits_all_seas)
            Nvisits_tot += Nvisits_all_seas

        return Nvisits_tot

    def Nvisits_seasons(self, nvisits_night, season_length, nseasons):
        """
        Total number of visits over a season

        Parameters
        ---------------
        nvisits_night: int
          number of visits per obs night
        season_length: float
          season length (in days)
        nseasons: int
          number of seasons of obs.
        """

        N_nights = nseasons*season_length/self.cadence

        return N_nights*nvisits_night
