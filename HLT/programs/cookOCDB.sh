#!/bin/bash
SOURCE=/cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB
TARGET=/home/drohr/HLT/HCDB_new3
FUTURE_RUN=500000
SOURCE_RUN=255283

#Download default CDB entries for future run
aliroot -l -b -q $ALICE_SOURCE/HLT/programs/downloadCDB.C"($FUTURE_RUN,\"local://$SOURCE\",\"local://$TARGET/tmp\",\"*/*/*\")"

#Extend run ranges of all these objects to infinity
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$TARGET/tmp ocdbTarget=$TARGET cdbEntries=*/*/* sourceRun=$FUTURE_RUN

#Clearn up temporary copy of default OCDB
rm -Rf $TARGET/tmp

#Remove objects not wanted in HCDB
rm -Rf $TARGET/TPC/Calib/CorrectionMaps*

#Create HLT HCDB configuration objects that differ from OCDB
aliroot -l -q -b $ALICE_SOURCE/HLT/exa/makeComponentConfigurationObject.C"(\"HLT/ConfigTPC/TPCHWClusterDecoder\", \"\", \"local://$TARGET\")"
aliroot -l -q -b $ALICE_SOURCE/HLT/exa/makeComponentConfigurationObject.C"(\"HLT/ConfigTPC/TPCHWClusterFinder\",
    \"-debug-level 0 -do-mc 0 -deconvolute-time 1 -deconvolute-pad 1 -flow-control 0 -single-pad-suppression 0 -bypass-merger 0 -cluster-lower-limit 10 -single-sequence-limit 0 -use-timebin-window 1 -merger-distance 4 -charge-fluctuation 0 -rcu2-data 1\",
    \"local://$TARGET\")"
#ATTENTION: the -rcu2-data flag is NOT set for the HCDB!!!!!! This is set by the chain configuration not by the HCDB!!!!!! Thus, the setting is not applied running aliroot!!!!!!

#Update TPC CalibDB
aliroot -b << EOF
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$TARGET");
  man->SetRun($SOURCE_RUN)
  AliCDBEntry* entry = man->Get("TPC/Calib/RecoParam");
  TObjArray* array = dynamic_cast<TObjArray*>(entry->GetObject());
  AliTPCRecoParam* param;
  TIter iter(array);
  while (param = dynamic_cast<AliTPCRecoParam*>(iter())) \
  { \
    param->SetUseCorrectionMap(kFALSE); \
    param->SetAccountDistortions(kFALSE); \
    param->SetUseTOFCorrection(kFALSE); \
    param->SetUseComposedCorrection(kTRUE); \
  }
  man->Put(entry);
EOF

#Fetch some special objects from one reference run and extend them, because default object does not exist
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="GRP/GRP/Data" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="GRP/CTP/Config" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="GRP/CTP/CTPtiming" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="GRP/CTP/Scalers" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="TPC/Calib/HighVoltage" sourceRun=$SOURCE_RUN

#Fetch some special objects from one reference run and extend them, because we want the latest object not the default object
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="TPC/Calib/AltroConfig" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="TPC/Calib/Temperature" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="ITS/Calib/BadChannelsSSD" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="ITS/Calib/SPDDead" sourceRun=$SOURCE_RUN

#Fetch some more special objects needed such that offline reco does not crash (used to compare offline reco runs with hlt recraw local)
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="MUON/Calib/HV" sourceRun=$SOURCE_RUN
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="MUON/Calib/OccupancyMap" sourceRun=$SOURCE_RUN

#Fetch CTP Aliases to make offline reco work with trigger class filter
$ALICE_SOURCE/HLT/programs/extendHLTOCDB.sh ocdbSource=$SOURCE ocdbTarget=$TARGET cdbEntries="GRP/CTP/Aliases" sourceRun=$SOURCE_RUN

rm -f $ALICE_SOURCE/OCDB/HLT/ConfigTPC/TPCFastTransform/*
aliroot -l -q -b $ALICE_SOURCE/HLT/TPCLib/macros/makeTPCFastTransformOCDBObject.C"(\"local://$TARGET\", $SOURCE_RUN, $SOURCE_RUN)"
SRCFILE=`ls $ALICE_SOURCE/OCDB/HLT/ConfigTPC/TPCFastTransform | tail -n 1`
aliroot -l -q -b $ALICE_SOURCE/HLT/programs/adjustOCDBObject.C"(\"$ALICE_SOURCE/OCDB/HLT/ConfigTPC/TPCFastTransform/$SRCFILE\", \"local://$TARGET\", 0)"
rm -f $ALICE_SOURCE/OCDB/HLT/ConfigTPC/TPCFastTransform/*
