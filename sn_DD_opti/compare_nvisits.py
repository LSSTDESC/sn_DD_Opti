import numpy as np
from optparse import OptionParser
from scipy.interpolate import interp1d
from sn_DD_opti import plt


def load_interp(fDir, prefix, nvisits,cadence=1):

    fName = '{}/{}'.format(fDir,fileName(prefix,nvisits))

    fi = np.load(fName)

    print(fi.dtype)

    idx = np.abs(fi['cadence']-cadence)<1.e-5
    sel = fi[idx]

    interp = interp1d(sel['z'],sel['Nvisits'],bounds_error=False, fill_value=0.)

    return interp
    

def fileName(prefix, nvisits):

    return '{}_{}.npy'.format(prefix,nvisits)


parser = OptionParser()

parser.add_option("--visitsDir", type=str, default='visits_files',
                  help="directory where visits files are located[%default]")
parser.add_option("--prefix", type=str, default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny',
                  help="prefix for the files [%default]")


opts, args = parser.parse_args()

visitsDir = opts.visitsDir
prefix = opts.prefix


Ny_vals = [20,30,40,60,80]

zvals = np.arange(0.5,0.95,0.01)

interp = {}

for Ny in Ny_vals:
    interp[Ny] = load_interp(visitsDir,prefix,Ny)

fig, ax = plt.subplots()

for key, vals in interp.items():
    print(vals(zvals))
    norm = interp[20](zvals)
    ax.plot(zvals,vals(zvals)/norm,label='N$_{}$={}'.format('y',key))

ax.legend()

plt.show()
