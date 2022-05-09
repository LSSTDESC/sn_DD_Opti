from version import __version__
import matplotlib.pyplot as plt
from collections import OrderedDict

filtercolors = dict(zip('ugrizy', ['b', 'c', 'g', 'y', 'r', 'm']))
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['font.size'] = 22
plt.rcParams['figure.titlesize'] = 22
plt.rcParams['legend.fontsize'] = 22
plt.rcParams['lines.linewidth'] = 3
#plt.rcParams['font.family'] = 'Liberation Serif'
plt.rcParams['font.family'] = 'Arial'
#plt.rcParams['font.family'] = 'Verdana'
#plt.rcParams['font.weight'] = 'bold'

#plt.rcParams['mathtext.default'] = 'regular'
#plt.rcParams['font.serif'] = 'Times'

linestyles = OrderedDict(
    [('solid',               (0, ())),
     ('loosely dotted',      (0, (1, 10))),
     ('dotted',              (0, (1, 5))),
     ('densely dotted',      (0, (1, 1))),

     ('loosely dashed',      (0, (5, 10))),
     ('dashed',              (0, (5, 5))),
     ('densely dashed',      (0, (5, 1))),

     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
     ('dashdotted',          (0, (3, 5, 1, 5))),
     ('densely dashdotted',  (0, (3, 1, 1, 1))),

     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
