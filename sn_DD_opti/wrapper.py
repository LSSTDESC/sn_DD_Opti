import pandas as pd
import numpy as np


class Mod_z:
    """
    Method to modify values of input file

    """

    def __init__(self, fName):

        # load the file

        tab = np.load(fName, allow_pickle=True)
       
        tabdf = pd.DataFrame.from_records(tab)
        
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
        return grp
