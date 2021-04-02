import os
from optparse import OptionParser
from sn_DD_opti.budget import GUI_Budget
from sn_DD_opti.visits import GUI_Visits


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

parser.add_option("--show", type="str", default='Visits',
                  help="GUI to visualize - visits or budget[%default]")
parser.add_option("--cadence", type="float", default=1.0,
                  help="cadence - for show Visits only [%default]")
parser.add_option("--Nvisits_max", type=int, default=300,
                  help="Max number of visits for display - for show Visits only [%default]")
parser.add_option("--zmax", type=float, default=0.95,
                  help="zmax for display - for show Visits only [%default]")
parser.add_option("--Ny_max", type=int, default=20,
                  help="y-band max number of visits [%default]")

opts, args = parser.parse_args()

Nvisits_z_file = 'Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
    opts.Ny_max)
Nvisits_z_fields_file = 'Nvisits_z_fields_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_Ny_{}.npy'.format(
    opts.Ny_max)

visitsDir = 'visits_files'

check_grab(visitsDir, [Nvisits_z_file, Nvisits_z_fields_file])

if opts.show == 'visits':
    myvisits = GUI_Visits(Nvisits_z_file,
                          cadence=opts.cadence,
                          nvisits_max=opts.Nvisits_max,
                          zmax=opts.zmax,
                          dir_files=visitsDir)

if opts.show == 'budget':
    mybud = GUI_Budget(Nvisits_z_file,
                       Nvisits_z_fields_file,
                       runtype='Nvisits_single',
                       dir_config='input',
                       dir_files=visitsDir)
