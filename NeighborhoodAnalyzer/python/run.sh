#!/bin/bash

for energy in 10 20 50 100
do
    for idx in 1 2 3 4 5 6 7 8 9 10
    do
        cmsRun Config.py $energy $idx
    done 
done
