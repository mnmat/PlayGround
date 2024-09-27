#!/bin/bash

eta="$1"
en="$2"
evts="$3"
idx="$4"
root="$5"
outdir="$6"
rescale="$7"

cmsRun Config.py $eta $en $evts $idx $root $outdir $rescale
