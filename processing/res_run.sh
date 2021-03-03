#!/bin/bash
#Argument Map:
#A = 1 -> Widths from table will be used
#= 2 -> Widths from table won't be used
#B = 1 -> Table will be updated
#= 2 -> Table won't be updated
#C = 1 -> Run code and quit
#= 2 -> Run code and show plots
#D = 0 -> Don't Apply any N_Missing Cut
#1 -> Apply N_Missing < 1
#= 2 -> apply N_Missing >= 1 Selection
#E = B Field, either 1.4 or 3.0 T at this time. DONT TRY ANOTHER B FIELD

#get B Field from input and parse for filenames
if (($(echo "$5 == 3.0" |bc -l) ));
then
  bfield="3T"
elif (($(echo "$5 == 1.4" |bc -l) ));
then
  bfield="1p4T" #the "p" stands for "point". i.e. "1 point 4 tesla"
fi

#Start doing stuff
rm wtables_ch_jet_analysis_theta-energy_resolution
make wtables_ch_jet_analysis_theta-energy_resolution
./wtables_ch_jet_analysis_theta-energy_resolution $1 $2 $3 $4 $5 ${bfield}_combined.root

#File name depends on if/how missing jet constituent cut is applied
if (( ${4:-9} == 0 )); 
then
  prefix="NoCuts"
elif (( ${4:-9} == 1 ));
then
  prefix="No_Missing_const"
elif (( ${4:-9} == 2 ));
then
  prefix="Missing_const"
elif ((${4} == 9));
then echo "Incorrect 4th argument, or not enough arguments given.\n"
fi

#scp it to mac mini
scp ../output/${prefix}_output_mom_res_${bfield}_combinedsigma_eta_5_p_6_.root fernando@136.152.45.216:/Users/fernando/Documents/Fun4All_analysis/analysis/notebook/
