import pandas as pd
from __init__ import plt, linestyles
from scipy.interpolate import interp1d
import numpy as np
import os
from optparse import OptionParser
from wrapper import Mod_z

lw = 3


class Show_Visits:
    """
    class to display nvisits vs redshift

    Parameters
    ---------------
    files_visits: str
      file name containing the number of visits vs redshift
      for a set of cadences
    cadence: float,opt
      cadence chosen (default: 1 day-1)
    nvisits_max: int, opt
      max number of vivits in display (default: 300)
    zmin: float, opt
      min z value (default: 0.1)
    zmax: float, opt
       max z value (default: 0.85)
    dir_files: str, opt
       location dir of the file

    """

    def __init__(self, ax, file_visits, cadence=1., nvisits_max=300, zmin=0.5, zmax=0.85, dir_files='input'):

        self.ax = ax
        self.cadence = cadence
        self.nvisits_max = nvisits_max
        self.zmin = zmin
        self.zmax = zmax
        # modify the visits
        # the idea here is to keep the number of visits constant for z>zlim where
        # zlim is the limit of the use of the band
        data = Mod_z('{}/{}'.format(dir_files, file_visits)).nvisits

        # select data for this cadence
        idx = np.abs(data['cadence']-self.cadence) < 1.e-5
        idx &= data['z'] >= self.zmin
        sel = data[idx]

        # prepare interpolators for plot
        self.zmin = np.min(sel['z'])
        #self.zmax = np.max(sel['z'])
        self.dictvisits = {}
        self.dictvisits['nvisits'] = self.interp_z(sel['z'], sel['Nvisits'])

        self.bands = 'grizy'
        self.colors = dict(zip(self.bands, ['c', 'g', 'y', 'r', 'm']))
        for b in self.bands:
            self.dictvisits['nvisits_{}'.format(b)] = self.interp_z(
                sel['z'], sel['Nvisits_{}'.format(b)])

        self.z_nvisits = self.interp_z(sel['Nvisits'], sel['z'])

    def interp_z(self, x, y):
        """
        Method to perform a 1d interpolator

        Parameters
        ---------------
        x: array(float)
           x values
        y: array (float)
           y values

        Returns
        ----------
        interp1d(x,y)

        """
        return interp1d(x, y, bounds_error=False, fill_value=0.)

    def plotNvisits(self):
        """
        Method to plot the number of visits vs redshift (background)

        """

        zstep = 0.005
        zvals = np.arange(self.zmin, self.zmax+zstep, zstep)

        lsb = dict(zip(['sum', 'g', 'r', 'i', 'z', 'y'], [
                   'solid', 'dotted', 'dashdot', (0, (3, 5, 1, 5, 1, 5)), (0, (5, 1)), (0, (3, 1, 1, 1, 1, 1))]))

        self.ax.plot(zvals, self.dictvisits['nvisits']
                     (zvals), color='k', label='sum', lw=3, ls=lsb['sum'])

        for io, b in enumerate(self.bands):
            key = 'nvisits_{}'.format(b)
            self.ax.plot(zvals, self.dictvisits[key](
                zvals), color=self.colors[b], label='${}$-band'.format(b), lw=3, ls=lsb[b])
        self.ax.grid()
        # self.ax.set_xlabel('$z_{\mathrm{complete}}$')
        self.ax.set_ylabel(r'$\mathrm{N_{visits}}$')
        # self.ax.legend()
        # self.ax.legend(bbox_to_anchor=(1.2, -0.1), ncol=1,
        #               fontsize=12,frameon=False, loc='lower right')
        self.ax.legend(bbox_to_anchor=(0.5, 1.35),
                       loc='upper center', ncol=3, frameon=False)
        self.ax.set_ylim(0,)
        if self.nvisits_max > 0:
            self.ax.set_ylim(ymax=self.nvisits_max)
        self.ax.set_xlim(xmax=self.zmax)

        self.ax.set_xticklabels([])

    def plotzlim(self, z=0.6):
        """
        Method to write the number of visits corresponding to a given redshift

        Parameters
        ---------------
        z : float
          redshift considered

        """
        fontsize = 18

        ylims = self.ax.get_ylim()
        nvisits = int(np.round(self.dictvisits['nvisits'](z)))
        yref = 0.85*ylims[1]
        scale = 0.12*ylims[1]
        nvstr = 'N$_{\mathrm{{visits}}}$'
        zstr = '$z_{\mathrm{complete}}$'
        # self.ax.text(0.55, yref, '{}= {}'.format(nvstr,
        #                                         nvisits))
        self.ax.text(0.60, 140, '{}= {}'.format(nvstr,
                                                nvisits))
        coeffa = dict(zip(self.bands, [0.52, 0.52, 0.52, 0.62, 0.62]))
        coeffb = dict(zip(self.bands, [0, 1, 2, 0, 1]))
        for io, b in enumerate(self.bands):
            key = 'nvisits_{}'.format(b)
            nvisits_b = int(np.round(self.dictvisits[key](z)))
            ff = 'N$_{\mathrm{visits}}^'+'{}'.format(b)+'$'
            self.ax.text(coeffa[b], 0.99*ylims[1]-scale*(coeffb[b]+1),
                         '{} ={}'.format(ff, nvisits_b), color=self.colors[b])

        # self.ax.text(0.9*z, np.min([1.5*nvisits, 280]),
        #             '{} = {}'.format(zstr, np.round(z, 2)))
        self.ax.text(0.93*z, np.min([1.75*nvisits, 280]),
                     '{} = {}'.format(zstr, np.round(z, 2)))

        # self.ax.arrow(z, np.min([1.4*nvisits, 270]), 0., np.max([-1.4*nvisits, -270]),
        #              length_includes_head=True, color='r',
        #              head_length=5, head_width=0.01)
        self.ax.plot([z]*2, [0., 300], color='k', ls='dotted', lw=lw)
        self.ax.plot([0., z], [nvisits]*2, color='k', ls='dotted', lw=lw)
        self.ax.set_ylim(0,)

    def plotnvisits(self, nvisits):
        """
        Method to draw the filter allocation and the redshift limit 
        corresponding to a total number of visits

        Parameters
        --------------
        nvisits: int
          total number of visits

        """
        # get the redshift from the total number of visits

        zlim = self.z_nvisits(nvisits)
        # self.plotzlim(np.round(zlim, 2))
        self.plotzlim(zlim)

        # self.ax.plot(self.ax.get_xlim(), [nvisits]*2,
        #             color='r', linestyle='--')


def plot_budget_zcomp(ax, zref=0.80):

    df = pd.read_csv('input/Nvisits_zcomp.csv')

    df['season_length'] = 180.
    df['nvisits_WFD'] = 2122176

    df['nvisits_DD'] = df['nvisits']*df['season_length']
    df['budget'] = df['nvisits_DD']/(df['nvisits_DD']+df['nvisits_WFD'])
    df['budget_per'] = 100.*df['nvisits_DD'] / \
        (df['nvisits_DD']+df['nvisits_WFD'])

    ww = 'budget'
    ww = 'budget_per'
    ax.plot(df['zcomp'], df[ww], color='r', lw=lw)
    interp = interp1d(df['zcomp'], df[ww], bounds_error=False, fill_value=0.)
    bud = interp(zref)

    ax.plot([zref]*2, [0., 2.5], color='k', ls='dotted', lw=lw)
    ax.plot([0., zref], [bud]*2, color='k', ls='dotted', lw=lw)

    #ax.text(0.67, bud+0.0005, 'budget = {}'.format(np.round(bud, 3)))
    ax.text(0.60, bud+0.05, 'budget = {} %'.format(np.round(bud, 1)))
    ax.grid()
    ax.set_xlim([0.50, 0.95])
    #ax.set_ylim([0., 0.025])
    ax.set_ylim([0., 2.5])
    ax.set_xlabel(r'$z_{complete}$')
    # ax.set_ylabel(r'budget/field/season')
    ax.set_ylabel(r'budget/field/season [%]')


def check_get_file(webPath, fDir, fName):
    """
    Function checking if a file is available
    If not, grab it from a web server

    Parameters
    ---------------
    webPath: str
       web path name
    fDir: str
      location dir of the file
    fName: str
      name of the file

    """

    if os.path.isfile('{}/{}'.format(fDir, fName)):
        return

    path = '{}/{}'.format(webPath, fName)
    cmd = 'wget --no-clobber --no-verbose {} --directory-prefix {}'.format(
        path, fDir)
    os.system(cmd)


def check_grab(visitsDir, fileList, webPath='https://me.lsst.eu/gris/DD_design'):

    # check if dir exist
    # if not, create it
    if not os.path.exists(visitsDir):
        os.mkdir(visitsDir)

    # check if files exist
    # if not grab it from a server

    for fi in fileList:
        check_get_file(webPath, visitsDir, fi)


parser = OptionParser()

parser.add_option("--show", type="str", default='visits',
                  help="GUI to visualize - visits or budget[%default]")
parser.add_option("--cadence", type="float", default=1.0,
                  help="cadence - for show Visits only [%default]")
parser.add_option("--Nvisits_max", type=int, default=300,
                  help="Max number of visits for display - for show Visits only [%default]")
parser.add_option("--zmax", type=float, default=0.95,
                  help="zmax for display - for show Visits only [%default]")
parser.add_option("--Ny_max", type=int, default=20,
                  help="y-band max number of visits [%default]")
parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")


opts, args = parser.parse_args()

Nvisits_z_file = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
    opts.Ny_max)
Nvisits_z_fields_file = 'Nvisits_z_fields_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
    opts.Ny_max)

visitsDir = opts.visitsDir

check_grab(visitsDir, [Nvisits_z_file, Nvisits_z_fields_file])


fig, ax = plt.subplots(figsize=(9, 16), nrows=2)

myvisits = Show_Visits(ax[0], Nvisits_z_file,
                       cadence=opts.cadence,
                       nvisits_max=opts.Nvisits_max,
                       zmax=opts.zmax,
                       dir_files=visitsDir)

myvisits.plotNvisits()
myvisits.plotnvisits(131)
ax[0].set_xlim([0.5, 0.95])
plot_budget_zcomp(ax[1])

plt.subplots_adjust(hspace=0.04)
plt.show()
