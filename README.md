### Batch jobs
Modify the output file name and path in new_default_data_vn.cc

mkdir condor_log

Make data file list that contains all paths to all files
e.g
```
find /eos/cms/store/group/comm_luminosity/gkrintir/highpT_jets/JetMET*/AK8PFJet500* -name '*root'> TreeList_Run3_2024.list

```
Split the file list:
```
mkdir Run3_tree_list
mkdir Run3_tree_list/list_25
cd Run3_tree_list/list_25
split -l25 -d -a 3 ../TreeList_Run3_2024.list list_job
```

Modify "batch/data_vn500.sh" to reflect the excutable and input file list that will be used. 
```
chmod +x data_vn500.sh
```

Modify "batch/OnOff.py" accordingly (e.g excutable name, file list name, line_count_frac=line_count/X, X needs to be the same as split -lX). This python script will divide jobs, give args etc. 

Submit condor jobs:
```Linux
python3 OnOff.py
```
