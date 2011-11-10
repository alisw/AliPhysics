/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_ROOT/PWG1/CalibMacros/Pass0/makeOCDB.C
  makeOCDB("105160");

*/

void makeOCDB(TString runNumberString, TString  ocdbStorage="")
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
  ConfigCalibTrain(runNumber, "raw://");

  // Steering Tasks - set output storage
  // DefaultStorage set already before - in ConfigCalibTrain.C
//ocdbStorage+="?se=ALICE::CERN::SE";

  AliCDBManager::Instance()->SetSpecificStorage("*/*/*",ocdbStorage.Data());

  // set OCDB storage
  if (ocdbStorage.Length()==0) ocdbStorage+="local://"+gSystem->GetFromPipe("pwd")+"/OCDB";

  // TPC part
  TFile fcalib("CalibObjects.root");
  AliTPCPreprocessorOffline proces;

  // switch on parameter validation
  proces.SetTimeGainRange(0.5,3.0);
  proces.SwitchOnValidation();

  // Make timegain calibration
  //proces.CalibTimeGain("CalibObjects.root", runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  proces.CalibTimeGain("CalibObjects.root", runNumber,runNumber,ocdbStorage);

  // Make vdrift calibration
  //proces.CalibTimeVdrift("CalibObjects.root",runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  proces.CalibTimeVdrift("CalibObjects.root",runNumber,runNumber,ocdbStorage);
  //
  // TOF part
  //
  AliTOFAnalysisTaskCalibPass0 calibTask;
  Printf("Calibrating TOF");
  calibTask.ProcessOutput("CalibObjects.root", ocdbStorage);
//
//

// T0 part
  AliT0PreprocessorOffline procesT0;
  // Make  calibration of channels offset
   procesT0.Process("CalibObjects.root",runNumber, runNumber, ocdbStorage);



   //TRD part
   AliTRDPreprocessorOffline procestrd;
   procestrd.SetLinearFitForVdrift(kTRUE);
   procestrd.SetMinStatsVdriftT0PH(600*10);
   procestrd.SetMinStatsVdriftLinear(50);
   procestrd.SetMinStatsGain(600);
   procestrd.SetLimitValidateNoData(40);
   procestrd.SetLimitValidateBadCalib(40);
   procestrd.Init("CalibObjects.root");
   Int_t versionVdriftUsed = procestrd.GetVersionVdriftUsed();
   Int_t subversionVdriftUsed = procestrd.GetSubVersionVdriftUsed();
   Int_t versionGainUsed = procestrd.GetVersionGainUsed();
   Int_t subversionGainUsed = procestrd.GetSubVersionGainUsed();
   Int_t versionExBUsed = procestrd.GetVersionExBUsed();
   Int_t subversionExBUsed = procestrd.GetSubVersionExBUsed();
   printf("version and subversion vdrift %d and %d\n",versionVdriftUsed,subversionVdriftUsed);
   printf("version and subversion gain %d and %d\n",versionGainUsed,subversionGainUsed);
   printf("version and subversion exb %d and %d\n",versionExBUsed,subversionExBUsed);
   procestrd.Process("CalibObjects.root",runNumber,runNumber,ocdbStorage);
   Int_t trdstatus = procestrd.GetStatus();
  
  
  //Mean Vertex
  AliMeanVertexPreprocessorOffline procesMeanVtx;
  procesMeanVtx.ProcessOutput("CalibObjects.root", ocdbStorage, runNumber);
	
  return;
}
