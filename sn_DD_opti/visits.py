import matplotlib
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from .wrapper import Mod_z
import tkinter as tk
from tkinter import font as tkFont
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import numpy as np
matplotlib.use('tkagg')

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['figure.titlesize'] = 15
plt.rcParams['legend.fontsize'] = 15
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.family'] = 'serif'


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

    def __init__(self, file_visits, cadence=1., nvisits_max=300, zmin=0.5, zmax=0.85, dir_files='input'):

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

        self.ax.plot(zvals, self.dictvisits['nvisits']
                     (zvals), color='k', label='sum', lw=2)

        for io, b in enumerate(self.bands):
            key = 'nvisits_{}'.format(b)
            self.ax.plot(zvals, self.dictvisits[key](
                zvals), color=self.colors[b], label='${}$-band'.format(b), lw=2)
        self.ax.grid()
        self.ax.set_xlabel('$z_{complete}$')
        self.ax.set_ylabel('$N_{visits}$')
        # self.ax.legend()
        # self.ax.legend(bbox_to_anchor=(1.2, -0.1), ncol=1,
        #               fontsize=12,frameon=False, loc='lower right')
        self.ax.legend(bbox_to_anchor=(-0.005, 1.1),
                       loc='upper left', ncol=6, frameon=False)
        self.ax.set_ylim(0,)
        if self.nvisits_max > 0:
            self.ax.set_ylim(ymax=self.nvisits_max)
        self.ax.set_xlim(xmax=self.zmax)

    def plotzlim(self, z=0.6):
        """
        Method to write the number of visits corresponding to a given redshift

        Parameters
        ---------------
        z : float
          redshift considered

        """
        fontsize = 15

        ylims = self.ax.get_ylim()
        nvisits = int(np.round(self.dictvisits['nvisits'](z)))
        yref = 0.9*ylims[1]
        scale = 0.1*ylims[1]
        nvstr = 'N_{visits}'
        zstr = 'z_{complete}'
        self.ax.text(0.6, yref, '${}$ = {}'.format(nvstr,
                                                   nvisits), fontsize=fontsize)
        for io, b in enumerate(self.bands):
            key = 'nvisits_{}'.format(b)
            nvisits_b = int(np.round(self.dictvisits[key](z)))
            self.ax.text(0.6, 0.8*ylims[1]-scale*io,
                         '${}^{}$ ={}'.format(nvstr, b, nvisits_b), fontsize=fontsize, color=self.colors[b])

        self.ax.text(0.95*z, 1.5*nvisits,
                     '${}$ = {}'.format(zstr, np.round(z, 2)), fontsize=fontsize)
        self.ax.arrow(z, 1.4*nvisits, 0., -1.4*nvisits,
                      length_includes_head=True, color='r',
                      head_length=5, head_width=0.01)
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

        self.ax.plot(self.ax.get_xlim(), [nvisits]*2,
                     color='r', linestyle='--')


class GUI_Visits(Show_Visits):
    """
    class building a GUI where plots may be shown and
    inherits from ShowVisits

    Parameters
    ---------------
    files_visits: str
      file name containing the number of visits vs redshift
      for a set of cadences
    cadence: float,opt
      cadence chosen (default: 1 day-1)
    nvisits_max: int, opt
      max number of visits for display (default: 300)
    zmin: float, opt
      min z value (default: 0.1)
    zmax: float, opt
       max z value (default: 0.85)
    dir_files: str, opt
       location dir of the file
    """

    def __init__(self, file_visits, cadence=1., nvisits_max=300, zmin=0.5, zmax=0.85, dir_files='input'):
        super().__init__(file_visits, cadence=cadence, nvisits_max=nvisits_max,
                         zmin=zmin, zmax=zmax, dir_files=dir_files)

        # build the GUI here
        root = tk.Tk()
        # figure where the plots will be drawn
        self.fig = plt.Figure(figsize=(15, 10), dpi=100)
        self.ax = self.fig.add_subplot(111)
        #leg = 'days$^{-1}$'
        leg = 'day'
        #self.fig.suptitle('cadence: {} {}'.format(int(self.cadence), leg))
        self.fig.subplots_adjust(right=0.8)
        self.ax.set_xlim(self.zmin, self.zmax)
        # define the figure canvas here
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH)
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)

        self.toolbar = NavigationToolbar2Tk(self.canvas, root)
        self.toolbar.update()
        # self.ax.cla()

        # plot the number of visits vs z
        self.plotNvisits()

        # common font
        helv36 = tkFont.Font(family='Helvetica', size=15, weight='bold')

        # building the GUI
        # frame
        button_frame = tk.Frame(master=root, bg="white")
        button_frame.pack(fill=tk.X, side=tk.BOTTOM, expand=False)
        button_frame.place(relx=.9, rely=.5, anchor="c")

        # entries
        ents = self.make_entries(button_frame, font=helv36)

        # buttons
        heightb = 3
        widthb = 6

        nvisits_button = tk.Button(
            button_frame, text="Nvisits \n from \n zlim", command=(lambda e=ents: self.updateData_z(e)),
            bg='yellow', height=heightb, width=widthb, fg='blue', font=helv36)

        z_button = tk.Button(button_frame, text="zlim \n from \n Nvisits", command=(
            lambda e=ents: self.updateData_nv(e)), bg='yellow', height=heightb, width=widthb, fg='red', font=helv36)

        quit_button = tk.Button(button_frame, text="Quit",
                                command=root.quit, bg='yellow',
                                height=heightb, width=widthb, fg='black', font=helv36)

        button_frame.columnconfigure(0, weight=1)
        button_frame.columnconfigure(1, weight=1)
        button_frame.columnconfigure(2, weight=1)

        nvisits_button.grid(row=2, column=0, sticky=tk.W+tk.E)
        z_button.grid(row=2, column=1, sticky=tk.W+tk.E)
        quit_button.grid(row=2, column=2, sticky=tk.W+tk.E)

        root.mainloop()

    def make_entries(self, frame, font):
        """
        Method to define entries to the GUI

        Parameters
        ---------------
        frame: tk.Frame
          frame where entries will be located
        font: tk.Font
          font for entries

        Returns
        ----------
        entries: dict
          dict of tk.Entry

        """
        tk.Label(frame, text='zcomplete', bg='white',
                 fg='blue', font=font).grid(row=0)
        tk.Label(frame, text='Nvisits', bg='white',
                 fg='red', font=font).grid(row=1)

        entries = {}
        entries['zlim'] = tk.Entry(frame, width=5, font=font)
        entries['Nvisits'] = tk.Entry(frame, width=5, font=font)
        # entries['zlim'].pack(ipady=3)
        entries['zlim'].insert(10, "0.7")
        entries['Nvisits'] .insert(10, "50")
        entries['zlim'].grid(row=0, column=1)
        entries['Nvisits'] .grid(row=1, column=1)

        return entries

    def updateData_z(self, entries):
        """
        Method to update the figure according to request made on entries
        The number of visits will be plotted here

        Parameters
        ---------------
        entries: dict of tk.Entry

        """

        # reset axes
        self.ax.cla()
        # plot Nvisits vs z
        self.plotNvisits()

        # get the redshift and plot if redshift>0
        z = float(entries['zlim'].get())
        if z > 0:
            self.plotzlim(z=z)

        # update canvas
        self.ax.set_xlim(self.zmin, self.zmax)
        self.canvas.draw()

    def updateData_nv(self, entries):
        """
        Method to update the figure according to request made on entries
        zlim and filter allocation will be plotted here.

        Parameters
        ---------------
        entries: dict of tk.Entry

        """

        # reset axes
        self.ax.cla()

        # plot Nvisits vs z
        self.plotNvisits()

        # get total number of visits and plot if nv>0
        nv = float(entries['Nvisits'].get())
        if nv > 0:
            self.plotnvisits(nvisits=nv)

        # update canvas
        self.ax.set_xlim(self.zmin, self.zmax)
        self.canvas.draw()
