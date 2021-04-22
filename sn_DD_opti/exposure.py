import yaml
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
import numpy.lib.recfunctions as rf
import copy


class DDF_Scenario:
    """

    """

    def __init__(self, config, fDir, fName, zfields=['COSMOS', 'XMM-LSS'], cadence=1):

        self.zlim_visit, self.visit_zlim = self.load(fDir, fName)

        # get the number of visits for the "z fixed" fields
        for field in zfields:
            zlim = config[field]['zlim']
            print('uuu', zlim)
            nvisits = self.visit_zlim[cadence](zlim)
            print('oooo', zlim, nvisits)
            config[field]['Nvisits'] = int(nvisits)

        print('estimated budget', self.nvisits_tot(config)/config['Nvisits'])

        configb = copy.deepcopy(config)
        self.res = self.budget_multiple_config(config, nADFS=1)
        self.resb = self.budget_multiple_config(config, nADFS=2)

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
          values: interp1d(Nvisits, z)

    """

        from wrapper import Mod_z
        #tab = np.load(fName, allow_pickle=True)
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
        return field['season_length']*field['Nseasons']*30./field['cadence']

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
        return field['Nvisits']*self.nights(field)

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
                bud += self.nvisits_field(config[field])

        return bud

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
        numpy array with three cols: Nvisits, zlim, budget
        """
        r = []
        print(config)
        for nv in range(10, 250, 10):
            for ff in fields:
                config[ff]['Nvisits'] = nv

            config['ADFS']['Nvisits'] *= nADFS
            # print(config)
            # print('estimated budget', nv, zlim_visit[1](
            #    nv), budget(config)/config['Nvisits'])
            r.append((nv, self.zlim_visit[1](
                nv), self.nvisits_tot(config)/config['Nvisits']))

        res = np.rec.fromrecords(r, names=['Nvisits', 'zlim', 'budget'])
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


parser = OptionParser()

parser.add_option('--config', type=str, default='config.yaml',
                  help='input config file [%default]')

opts, args = parser.parse_args()

config = yaml.load(open(opts.config))

print(config)

"""
# number of AGN visits
AGN = config['AGN']
Nvisits_AGN = int(AGN['Nvisits']*nights(AGN))

print('AGN: frac:', Nvisits_AGN/config['Nvisits'], 'Nvisits:', Nvisits_AGN)

# Number of DDF visits alloted to SN

Nvisits_DDF = int(config['frac_DDF']*config['Nvisits']-Nvisits_AGN)

print('Number of visits allowed for SN', Nvisits_DDF)
# this number of visits is equal to sum(Nvisits_per_night)*30*season_length*Nseasons/cadence

COSMOS = config['COSMOS']
DDF_not_COSMOS = config['DDF_not_COSMOS']

Nvisits_per_night = Nvisits_DDF/(nights(COSMOS)+nights(DDF_not_COSMOS))

print('Max allowed number of visits per night', int(Nvisits_per_night))

if 'Nvisits' in COSMOS.keys():
    Nvisits_new = COSMOS['Nvisits']*nights(COSMOS)
    Nvisits_new += DDF_not_COSMOS['Nvisits']*nights(DDF_not_COSMOS)
    print('DDF SN only (total):', Nvisits_new/config['Nvisits'])

r = []
for nvisits in range(10, 500):
    bud = budget(COSMOS, nvisits)+budget(DDF_not_COSMOS, nvisits)
    bud /= config['Nvisits']
    r.append((nvisits, bud))

res = np.rec.fromrecords(r, names=['nvisits', 'budget'])

plt.plot(res['nvisits'], res['budget'])

ninterp = interp1d(res['budget'], res['nvisits'])
for budget in [0.06, 0.1, 0.15]:
    print(budget, ninterp(budget))

plt.show()
"""

dd_scen = DDF_Scenario(config, 'visits_files',
                       'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_80.npy')

res = dd_scen.res
resb = dd_scen.resb

# res = budget_config(config)
fig, ax = plt.subplots()
color = dict(zip(range(1, 4), ['k', 'r', 'b']))
nseas_cosmos = config['COSMOS']['Nseasons']
nseas_xmm = config['XMM-LSS']['Nseasons']
nseas_adfs = config['ADFS']['Nseasons']
nseas_field = config['CDFS']['Nseasons']

legm = 'COSMOS$^{0.9}_'+str(nseas_cosmos) + \
    '$, XMM-LSS$^{0.9}_'+str(nseas_xmm)+'$'
ADFS = 'ADFS$_'+str(nseas_adfs)+'$'
field = 'field$_'+str(nseas_field)+'$'
fields = 'fields$_'+str(nseas_field)+'$'

legg = dict(
    zip(range(1, 4), [ADFS+'(1p)', ADFS+'(1p)+1 '+field, ADFS+'(1p)+2 '+fields]))
leggb = dict(
    zip(range(1, 4), [ADFS+'(2p)', ADFS+'(2p)+1 '+field, ADFS+'(2p)+2 '+fields]))
for nf in np.unique(res['Nf']):
    idx = res['Nf'] == nf
    sel = res[idx]
    idxb = resb['Nf'] == nf
    selb = resb[idxb]
    legt = '{}+{}'.format(legm, legg[nf])
    if nf > 1:
        legt += '/{}'.format(leggb[nf-1])
    ax.plot(sel['zlim'], sel['budget'], color=color[nf],
            label=legt)
    if nf == 3:
        legt = '{}+{}'.format(legm, leggb[nf])
        ax.plot(selb['zlim'], selb['budget'],
                color=color[nf], ls='--', label=legt)
ax.grid()

ax.set_xlim([0.65, 0.9])
ax.set_ylim([0.05, 0.15])
ax.set_xlabel('$z_{complete}$', fontsize=12)
ax.set_ylabel('Budget', fontsize=12)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)

ax.legend(fontsize=12, frameon=True)

plt.show()
