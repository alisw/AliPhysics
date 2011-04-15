/*
    Make OCDB entries for the calibration using tracks:
    The specific storage is set in the STEERING macro.


    void makeOCDBTPC(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity(),TString inputFile="CalibObjects.root",TString ocdbStorage)
    Example:

    .L makeOCDBTPC.C 
    makeOCDBTPC(0, AliCDBRunRange::Infinity(), "CalibObjects.root","");
 

*/



void makeOCDBTPC(Int_t startRun, Int_t endRun=AliCDBRunRange::Infinity(),TString inputFile="CalibObjects.root", TString ocdbStorage=""){
  //
  //
  //
  gSystem->Load("libSTEER");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTAT");
  gSystem->Load("libTPCcalib");
  gROOT->LoadMacro("CalibTimeVdrift.C");
  gROOT->LoadMacro("CalibTimeGain.C");
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";
  //
  // Make gain calibration
  //
  TFile fcalib(inputFile.Data());
  //
  //
  //
  CalibTimeGain(inputFile.Data(), startRun,endRun,ocdbStorage);
  //
  // Make vdrift calibration
  //
  CalibTimeVdrift(inputFile.Data(),startRun,AliCDBRunRange::Infinity(),ocdbStorage);
}
