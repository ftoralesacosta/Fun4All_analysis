rm wtables_ch_jet_analysis_theta-energy_resolution
make wtables_ch_jet_analysis_theta-energy_resolution
./wtables_ch_jet_analysis_theta-energy_resolution $1 $2 $3 $4 1.4 1p4T_combined.root
scp ../output/No_Missing_const_output_mom_res_1p4T_combinedsigma_eta_5_p_6_.root fernando@169.229.123.248:/Users/fernando/Documents/Fun4All_analysis/analysis/notebook/
