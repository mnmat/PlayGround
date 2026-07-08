#! bin/bash

# This script can be used to write .txt files containing arguments for the create_samples.sh executable for htcondor.
# It can also be used to split the total number of events into more manageable sizes to utilize parallelization.
# The arguments are eta, energy, number of events, and path to the input file

# Set input variables
cap="zpos"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Trajectories/RK_Rx/0_PU/"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Geant4ePropagator/Geant4ePropagator/0_PU/"
#root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/DebugLowEfficiency/200_PU/"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.2/Purity//0_PU"
root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.2/Scan/0_PU"
particles="singlemuon"
generator="flatPtGun"
position="hgcalCenter"
tevts=10000 # Number of total events
fevts=1000 # Number of events per file
end=$((tevts/fevts))

gf="files.txt"
rm $gf
touch $gf


for xi in 5
do
# Create txt file
    #f="files_"$xi".txt"
    #rm $f
    #touch $f
    #outdir="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Geant4ePropagator/RungeKutta/0_PU/Smoothing_G4ePropagator/x"$xi"/0_PU/"
    outdir=$root"Scan/x"$xi
    outdir=$root
    for eta in 17 19 21 23 25 27
    do
        for energy in 3 5 10 15 20 50 100
        do
            for i in $(seq 1 $end)
            do
                #echo "$eta, $energy, $fevts, $i, $root, $outdir, $xi" >> $f
                echo "$eta,$energy,$fevts,$i,$root,$outdir,$xi,$generator" >> $gf
            done
        done
    done
done
