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
#../bin/new_default_data_vn.exe ./Run3_tree_list/2023/list_25/list_job$1 0 1
#../bin/new_default_data_vn.exe ./list_test.txt 0 1
