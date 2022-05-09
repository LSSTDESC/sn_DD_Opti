#!/bin/bash

#figs 20 and 21
###python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly_zspectro_x1c_1sigma_wa_sigp_uy --nspectro 20000tenyears --var_to_plot detfom --legy "DETF FoM [95% C.L.]"

###python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly_zspectro_x1c_1sigma_wa_sigp_uy --nspectro pfscurrent --var_to_plot detfom --legy "DETF FoM [95% C.L.]"



###python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly_zspectro_x1c_1sigma_sigp_uy --nspectro 20000tenyears

#python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly_zspectro_x1c_1sigma_sigp_uy,cosmo_files_yearly_zspectro_x1c_1sigma_wa_sigp_uy --nspectro pfscurrent 
python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly_zspectro_x1c_1sigma_sigp_uy,cosmo_files_yearly_zspectro_x1c_1sigma_wa_sigp_uy --nspectro 20000tenyears

# fig 19
#python sn_DD_opti/plot_survey_time.py --cosmoDir cosmo_files_yearly,cosmo_files_yearly --nspectro nosyste --var_to_plot nsn_ultra,nsn_dd --legy '$N_{SN}^{ultra-deep}$,$N_{SN}^{deep}$'

#fig 18.
#python sn_DD_opti/plot_survey_nspectro.py --cosmoDir cosmo_files_yearly --nspectro nosyste
