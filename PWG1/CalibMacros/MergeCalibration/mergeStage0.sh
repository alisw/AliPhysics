#
# Merging Stage 0
#
# Merging in separate process for each detector
# log files - merge<detector>.log and syswatch<detector>.log stored also
# per detector
#
# 
# arguments of the merge macros:
# 
# 1 - input list of the files
# 2 - name of output root file with merged calibration
# 3 - name filter - what to merge
# 4 - reject mask
# 5 - flag single key for Obj array
#
#
# merge T0
#
aliroot -b -q $ALICE_ROOT/PWG1/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"CalibObjectsT0.root\",\"T0Calib\",\"AliESDfriends\",kTRUE\) >> mergeT0.log
cp syswatch.log syswatchT0.log
#
# merge TOF
#
aliroot -b -q $ALICE_ROOT/PWG1/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"CalibObjectsTOF.root\",\"TOFHistos\",\"AliESDfriends\",kTRUE\) >> mergeTOF.log
cp syswatch.log syswatchTOF.log
#
# merge TPC
#
aliroot -b -q $ALICE_ROOT/PWG1/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"CalibObjectsTPC.root\",\"TPCCalib\",\"AliESDfriends\",kFALSE\) >> mergeTPC.log
cp syswatch.log syswatchTPC.log
#
# merge TRD
#
aliroot -b -q $ALICE_ROOT/PWG1/CalibMacros/MergeCalibration/mergeCustom.C\(\"calib.list\",\"CalibObjectsTRD.root\",\"TRDcalib\",\"AliESDfriends\",kFALSE\)  >> mergeTRD.log
cp syswatch.log syswatchTRD.log
