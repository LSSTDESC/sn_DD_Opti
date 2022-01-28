#!/bin/bash

#Early Deep Rolling

python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot nsn_DD --Ny 40
python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot sigma_w --Ny 40
#python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_0.60_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot nsn_DD --Ny 40
#python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_0.60_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot sigma_w --Ny 40
#python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_0.65_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot sigma_w --Ny 40
#python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_2_2_mini_0.65_Ny_40 --xaxis nseasons_ultra_unique --yaxis zcomp_ultra_unique --var_to_plot nsn_DD --Ny 40

#Deep Rolling 10 Years
python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_0.80_0.80_2_2_Ny_40 --var_to_plot nsn_DD --Ny 40
python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_deep_rolling_0.80_0.80_2_2_Ny_40 --var_to_plot sigma_w --Ny 40

#universal 
python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_universal_10_Ny_40 --Ny 40 --var_to_plot nsn_DD
python sn_DD_opti/deep_rolling_contour.py --cosmoDir cosmo_files_10years --cosmoFile cosmoSN_universal_10_Ny_40 --Ny 40 --var_to_plot sigma_w
