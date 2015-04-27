// Macro to extract calibration summary information
// ConfigOCDB.C  macro has to be present in working directory
//
void CalibSummary(Int_t irun, const char* ocdbStorage){
  //
  //
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/TPC/macros/ConfigOCDB.C");
  ConfigOCDB(irun,ocdbStorage);  
  AliTPCcalibSummary *calibSummary = new AliTPCcalibSummary;
  calibSummary->ProcessRun(irun);
  delete calibSummary;
}
  
