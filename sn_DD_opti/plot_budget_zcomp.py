import pandas as pd
from __init__ import plt, linestyles
from scipy.interpolate import interp1d
import numpy as np

df = pd.read_csv('input/Nvisits_zcomp.csv')

df['season_length'] = 180.
df['nvisits_WFD'] = 2122176

df['nvisits_DD'] = df['nvisits']*df['season_length']
df['budget'] = df['nvisits_DD']/(df['nvisits_DD']+df['nvisits_WFD'])

fig, ax = plt.subplots(figsize=(12, 9))
ax.plot(df['zcomp'], df['budget'], color='r', lw=2)
interp = interp1d(df['zcomp'], df['budget'], bounds_error=False, fill_value=0.)
zref = 0.80
bud = interp(zref)

ax.plot([zref]*2, [0., bud], color='k', ls='dotted')
ax.plot([0., zref], [bud]*2, color='k', ls='dotted')

ax.text(0.67, bud+0.0005, 'budget = {}'.format(np.round(bud, 3)))
ax.grid()
ax.set_xlim([0.60, 0.9])
ax.set_ylim([0., 0.02])
ax.set_xlabel(r'$z_{complete}$')
ax.set_ylabel(r'budget/field/season')
plt.show()