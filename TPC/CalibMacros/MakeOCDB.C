void MakeOCDB(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity(),TString inputFile="CalibObjectsTrain1.root"){
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
  TFile fcalib(inputFile.Data());
  AliTPCcalibTimeGain * gain = ( AliTPCcalibTimeGain *)fcalib.Get("calibTimeGain");    
  if (gain==0) return;
  CalibTimeGain(gain, ocdbStorage.Data(),startRun,endRun,kTRUE);
  //
  // Make vdrift calibration
  //
  CalibTimeVdriftGlobal(inputFile.Data(),startRun,AliCDBRunRange::Infinity());
  //
  // Make calibration plot
  //
  Int_t run=endRun;
  ocdbStorage="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/TimeDrift",ocdbStorage.Data());  
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("TPC/Calib/TimeDrift",run);
  TObjArray * arr = (TObjArray*)entry->GetObject();
  TObjArray *picArray = new TObjArray;
  MakeDefaultPlots(arr,picArray);
  TFile fdrift("vdrift.root","recreate");
  picArray->Write("drift Plot");
  fdrift.Close();  
}
