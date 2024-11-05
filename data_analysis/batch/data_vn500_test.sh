#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/

#cmsenv
eval `scramv1 runtime -sh`
#cmssw-el7
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/JetCollectivity/data_analysis/batch
echo PWD: $PWD
echo $HOSTNAME
../new_default_data_vn ./list_test.txt 0 1
