import os
from optparse import OptionParser
from sn_DD_opti.budget import GUI_Budget
from sn_DD_opti.showvisits import GUI_Visits

parser = OptionParser()

parser.add_option("--show", type="str", default='Visits',
                  help="GUI to visualize - Visits or Budget[%default]")
parser.add_option("--cadence", type="float", default=3.0,
                  help="Cadence - for show Visits only [%default]")
parser.add_option("--Nvisits_z_file", type="str", default='Nvisits_z_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_selb.npy',
                  help="File with Nvisits vs redshift - median field [%default]")
parser.add_option("--Nvisits_z_field_file", type="str", default='Nvisits_z_fields_-2.0_0.2_error_model_ebvofMW_0.0_nvisits_selb.npy',
                  help="File with Nvisits vs redshift for each field [%default]")

opts, args = parser.parse_args()

nvisits_cadence = opts.Nvisits_z_file
nvisits_cadence_season = opts.Nvisits_z_field_file


if opts.show == 'Visits':
    myvisits = GUI_Visits(nvisits_cadence, cadence=opts.cadence)

if opts.show == 'Budget':
    mybud = GUI_Budget(nvisits_cadence,
                       nvisits_cadence_season,
                       runtype='Nvisits_single', dir_config='input')
