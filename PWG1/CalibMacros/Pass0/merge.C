/*
  merge output calib objects on Alien
  using AliFileMerger functionality

  Directory with runCalibTrain output: outputDir
  pattern: AliESDfriends_v1.root 
  Output file name: CalibObjects.root

  Example: 
  .L $ALICE_ROOT/ANALYSIS/CalibMacros/MergeCalibration/merge.C
  merge("alien:///alice/cern.ch/user/j/jotwinow/CalibTrain/output","AliESDfriends_v1.root");
*/

void merge(const char* outputDir, const char* pattern) {
  //
  // load libraries
  //
  gROOT->Macro("LoadLibraries.C");
  //
  TH1::AddDirectory(0);

  //
  AliFileMerger merger;
  merger.AddReject("esdFriend"); // do not merge 

  // local 
  // merger.IterTXT("calib.list","CalibObjects.root",kFALSE);
	 
  // alien 
  merger.IterAlien(outputDir, "CalibObjects.root", pattern);

return;
}
