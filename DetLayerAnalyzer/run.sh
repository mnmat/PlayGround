#!/bin/bash

eta="$1"
en="$2"
evts="$3"
idx="$4"
root="$5"
outdir="$6"
rescale="$7"
generator="$8"

#eta=17
#en=100
#evts=1
#idx=1
#root="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Geant4ePropagator/"
#outdir="/eos/cms/store/group/dpg_hgcal/comm_hgcal/mmatthew/PatternRecognitionByKalmanFilter/CMSSW_14_1_0_pre2/KFv0.1/Geant4ePropagator"
#rescale=1
#generator="flatPtGun"

cmsRun Config.py $eta $en $evts $idx $root $outdir $rescale $generator
