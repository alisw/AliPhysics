/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_ROOT/PWG1/CalibMacros/Pass0/makeOCDB.C
  makeOCDB("105160");

*/

//__________________________________________________________________

void makeOCDB(TString runNumberString, TString ocdbStorage = "")
{
  makeOCDB("CalibObjects.root", "ALL", runNumberString, ocdbStorage);
}

//___________________________________________________________________

void makeOCDB(const Char_t *filename, TString component, TString runNumberString, TString ocdbStorage = "")
{
  //
  // extract TPC OCDB entries
  //
  gROOT->Macro("LoadLibraries.C");
  gROOT->LoadMacro("ConfigCalibTrain.C");

  // switch off log info
  AliLog::SetClassDebugLevel("AliESDEvent",0);

  // config GRP
  Int_t runNumber = runNumberString.Atoi();
  printf("runNumber from runCalibTrain = %d\n",runNumber);

  /* configCalibTrain only if needed */
  if (component == "TPC" || component == "TRD" || component == "ALL")
    ConfigCalibTrain(runNumber, "raw://");
  else
    AliCDBManager::Instance()->SetDefaultStorage("raw://");
  
  // Steering Tasks - set output storage
  // DefaultStorage set already before - in ConfigCalibTrain.C
//ocdbStorage+="?se=ALICE::CERN::SE";

  AliCDBManager::Instance()->SetSpecificStorage("*/*/*",ocdbStorage.Data());

  // set OCDB storage
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";

  
  /* makeOCDB for selected component */
  if (component == "TPC" || component == "ALL")
    makeOCDB_TPC(filename, runNumber, ocdbStorage);
  if (component == "TOF" || component == "ALL")
    makeOCDB_TOF(filename, runNumber, ocdbStorage);
  if (component == "T0" || component == "ALL")
    makeOCDB_T0(filename, runNumber, ocdbStorage);
  if (component == "TRD" || component == "ALL")
    makeOCDB_TRD(filename, runNumber, ocdbStorage);
  if (component == "MeanVertex" || component == "ALL")
    makeOCDB_MeanVertex(filename, runNumber, ocdbStorage);
  
  gSystem->Exec(Form("touch %s_ocdb_done", component.Data()));
  return;
}

//___________________________________________________________________

void makeOCDB_TPC(const Char_t *filename, Int_t runNumber, TString ocdbStorage)
{

  // TPC part
  TFile fcalib(filename);
  AliTPCPreprocessorOffline proces;

  // switch on parameter validation
  proces.SetTimeGainRange(0.5,3.0);
  proces.SwitchOnValidation();

  // Make timegain calibration
  //proces.CalibTimeGain(filename, runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  proces.CalibTimeGain(filename, runNumber,runNumber,ocdbStorage);

  // Make vdrift calibration
  //proces.CalibTimeVdrift(filename,runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  proces.CalibTimeVdrift(filename,runNumber,runNumber,ocdbStorage);

}

//___________________________________________________________________

void makeOCDB_TOF(const Char_t *filename, Int_t runNumber, TString ocdbStorage)
{
  AliTOFAnalysisTaskCalibPass0 calibTask;
  Printf("Calibrating TOF");
  calibTask.ProcessOutput(filename, ocdbStorage);
}

//___________________________________________________________________

void makeOCDB_T0(const Char_t *filename, Int_t runNumber, TString ocdbStorage)
{
  // T0 part
  AliT0PreprocessorOffline procesT0;
  // Make  calibration of channels offset
  procesT0.Process(filename,runNumber, runNumber, ocdbStorage);
}

//___________________________________________________________________

void makeOCDB_TRD(const Char_t *filename, Int_t runNumber, TString ocdbStorage)
{
   //TRD part
   AliTRDPreprocessorOffline procestrd;
   procestrd.SetLinearFitForVdrift(kTRUE);
   procestrd.SetMinStatsVdriftT0PH(600*10);
   procestrd.SetMinStatsVdriftLinear(50);
   procestrd.SetMinStatsGain(600);
   procestrd.Init(filename);
   Int_t versionVdriftUsed = procestrd.GetVersionVdriftUsed();
   Int_t subversionVdriftUsed = procestrd.GetSubVersionVdriftUsed();
   Int_t versionGainUsed = procestrd.GetVersionGainUsed();
   Int_t subversionGainUsed = procestrd.GetSubVersionGainUsed();
   Int_t versionExBUsed = procestrd.GetVersionExBUsed();
   Int_t subversionExBUsed = procestrd.GetSubVersionExBUsed();
   printf("version and subversion vdrift %d and %d\n",versionVdriftUsed,subversionVdriftUsed);
   printf("version and subversion gain %d and %d\n",versionGainUsed,subversionGainUsed);
   printf("version and subversion exb %d and %d\n",versionExBUsed,subversionExBUsed);
   procestrd.Process(filename,runNumber,runNumber,ocdbStorage);
   Int_t trdstatus = procestrd.GetStatus();

}

//___________________________________________________________________

void makeOCDB_MeanVertex(const Char_t *filename, Int_t runNumber, TString ocdbStorage)
{
  //Mean Vertex
  AliMeanVertexPreprocessorOffline procesMeanVtx;
  procesMeanVtx.ProcessOutput(filename, ocdbStorage, runNumber);
}

//___________________________________________________________________

