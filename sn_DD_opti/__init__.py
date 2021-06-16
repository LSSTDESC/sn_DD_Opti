from .version import __version__
import matplotlib.pyplot as plt

filtercolors = dict(zip('ugrizy', ['b', 'c', 'g', 'y', 'r', 'm']))
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 22
plt.rcParams['font.size'] = 22
plt.rcParams['figure.titlesize'] = 22
plt.rcParams['legend.fontsize'] = 22

#plt.rcParams['font.family'] = 'Liberation Serif'
plt.rcParams['font.family'] = 'Arial'
#plt.rcParams['font.family'] = 'Verdana'
#plt.rcParams['font.weight'] = 'bold'

#plt.rcParams['mathtext.default'] = 'regular'
#plt.rcParams['font.serif'] = 'Times'
