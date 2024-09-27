#! bin/bash

# This script can be used to write .txt files containing arguments for the create_samples.sh executable for htcondor.
# It can also be used to split the total number of events into more manageable sizes to utilize parallelization.
# The arguments are eta, energy, number of events, and path to the input file

# Set input variables
cap="zpos"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Trajectories/Rk_Rx2/0_PU/"
particles="singlemuon"
producer="flatEGun"
position="hgcalCenter"
tevts=500 # Number of total events
fevts=100 # Number of events per file
end=$((tevts/fevts))

gf="files.txt"
rm $gf
touch $gf


for xi in 1 2 5 10 15 20 100
do
# Create txt file
    #f="files_"$xi".txt"
    #rm $f
    #touch $f
    outdir="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Trajectories/Rk_Rx2/x"$xi"/0_PU/"
    for eta in 17
    do
        for energy in 10 100
        do
            for i in $(seq 1 $end)
            do
                #echo "$eta, $energy, $fevts, $i, $root, $outdir, $xi" >> $f
                echo "$eta, $energy, $fevts, $i, $root, $outdir, $xi" >> $gf
            done
        done
    done
done
