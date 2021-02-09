rm wtables_ch_jet_analysis_theta-energy_resolution
make wtables_ch_jet_analysis_theta-energy_resolution
./wtables_ch_jet_analysis_theta-energy_resolution 2 1 1 $1 combined_electrons+jets.root
scp ../output/No_Missing_const_output_mom_res_combined_electrons+jetssigma_eta_5_p_6_.root fernando@169.229.123.248:~/Documents/Fun4All_analysis/output/
