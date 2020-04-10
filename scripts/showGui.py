import os
from optparse import OptionParser
from sn_DD_opti.budget import DD_Budget
from sn_DD_opti.showvisits import ShowVisits

parser = OptionParser()

parser.add_option("--show", type="str", default='Visits',
                  help="GUI to visualize - Visits or Budget[%default]")
parser.add_option("--cadence", type="float", default=3.0,
                  help="Cadence - for show_Visits only [%default]")

opts, args = parser.parse_args()

nvisits_cadence = 'Nvisits_cadence_Nvisits_median_m5_filter.npy'
nvisits_cadence_season = 'Nvisits_cadence_Nvisits_median_m5_field_filter_season.npy'


if opts.show == 'Visits':
    myvisits = ShowVisits(nvisits_cadence, cadence=opts.cadence)

if opts.show == 'Budget':
    mybud = DD_Budget(nvisits_cadence,
                      nvisits_cadence_season,
                      runtype='Nvisits_single', dir_config='input')
