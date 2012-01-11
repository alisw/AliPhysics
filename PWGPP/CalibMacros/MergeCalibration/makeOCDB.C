/*
  macro to extract the OCDb entries

  Example:
  .L $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/makeOCDB.C
  makeOCDB("105160");

*/

void makeOCDB(TString runNumberString, TString  ocdbStorage="")
{
  gROOT->Macro("LoadLibraries.C");
  gROOT->LoadMacro("ConfigCalibTrain.C");
 
  // detector tasks
  
  gROOT->LoadMacro("makeOCDBTPC.C");
  
  AliLog::SetClassDebugLevel("AliESDEvent",19);
  Int_t runNumber = runNumberString.Atoi();
  printf("runNumber from runCalibTrain = %d\n",runNumber);
  ConfigCalibTrain(runNumber, "raw://");
  
  // Steering Tasks  -set output storage
  // DefaultStorage set already before - in ConfigCalibTrain.C
  // 
  AliCDBManager::Instance()->SetSpecificStorage("*/*/*",ocdbStorage.Data());
  
  
  // Detector Tasks
  makeOCDBTPC(runNumber,runNumber, "CalibObjects.root" , ocdbStorage);
  
  
  return;
}
