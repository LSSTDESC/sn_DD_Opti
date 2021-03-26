import yaml
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def nights(field):

    return field['season_length']*field['Nseasons']*field['Nfields']*30./field['cadence']


def budget(fields, nvisits):

    return nvisits*nights(fields)


parser = OptionParser()

parser.add_option('--config', type=str, default='config.yaml',
                  help='input config file [%default]')

opts, args = parser.parse_args()

config = yaml.load(open(opts.config))

print(config)

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
