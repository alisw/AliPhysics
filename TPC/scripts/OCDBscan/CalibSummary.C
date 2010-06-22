// Macro to extract calibration summary information
// ConfigOCDB.C  macro has to be present in working directory
//
void CalibSummary(Int_t irun){
  //
  //
  //
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");  
  gROOT->LoadMacro("ConfigOCDB.C");
  ConfigOCDB(irun);  
  AliTPCcalibSummary *calibSummary = new AliTPCcalibSummary;
  calibSummary->ProcessRun(irun);
  delete calibSummary;
   
}
  
