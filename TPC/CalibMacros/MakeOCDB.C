void MakeOCDB(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity()){
  //
  //
  //
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  gSystem->AddIncludePath("-I$ALICE_ROOT/STEER");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC");

  gROOT->LoadMacro("$ALICE_ROOT/TPC/CalibMacros/CalibTimeVdrift.C+");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/CalibMacros/CalibTimeGain.C+");
  ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  // Make gain calibration
  //
  TFile fcalib("CalibObjectsTrain1.root");
  AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGain");    
  CalibTimeGain(gain, ocdbStorage.Data(),startRun,endRun,kTRUE);
  //
  // Make vdrift calibration
  //
  CalibTimeVdriftGlobal();


}
