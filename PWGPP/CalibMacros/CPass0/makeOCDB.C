/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_ROOT/PWGPP/CalibMacros/CPass0/makeOCDB.C
  makeOCDB("105160");

*/

void PrintDetectorStatus();



void makeOCDB(TString runNumberString, TString  ocdbStorage="")
{
  //
  // extract TPC OCDB entries
  //
  gROOT->Macro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/CalibMacros/CPass0/ConfigCalibTrain.C");

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
  AliTPCPreprocessorOffline procesTPC;

  // switch on parameter validation
  procesTPC.SetTimeGainRange(0.5,3.0);
  procesTPC.SwitchOnValidation();

  // Make timegain calibration
  //proces.CalibTimeGain("CalibObjects.root", runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  procesTPC.CalibTimeGain("CalibObjects.root", runNumber,runNumber,ocdbStorage);

  // Make vdrift calibration
  //proces.CalibTimeVdrift("CalibObjects.root",runNumber,AliCDBRunRange::Infinity(),ocdbStorage);
  procesTPC.CalibTimeVdrift("CalibObjects.root",runNumber,runNumber,ocdbStorage);
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
   procesT0.setDArun(179000);
   procesT0.Process("CalibObjects.root",runNumber, runNumber, ocdbStorage);



   //TRD part
   AliTRDPreprocessorOffline procestrd;
   procestrd.SetLinearFitForVdrift(kTRUE);
   procestrd.SetMinStatsVdriftT0PH(600*10);
   procestrd.SetMinStatsVdriftLinear(50);
   procestrd.SetMinStatsGain(600);
   procestrd.SetLimitValidateNoData(60);
   procestrd.SetLimitValidateBadCalib(60);
   procestrd.SetAlternativeDriftVelocityFit(kTRUE);
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
  
  
   //Mean Vertex
   AliMeanVertexPreprocessorOffline procesMeanVtx;
   procesMeanVtx.ProcessOutput("CalibObjects.root", ocdbStorage, runNumber);

   //
   // Print calibration status into the stdout
   //
   Int_t trdStatus = procestrd.GetStatus();
   Int_t tofStatus = calibTask.GetStatus();
   Int_t tpcStatus = ((procesTPC.ValidateTimeDrift() || procesTPC.ValidateTimeGain())==kFALSE);
   //
   printf("\n\n\n\n");
   printf("CPass0 calibration status\n");
   printf("TRD calibration status=%d\n",trdStatus);
   printf("TOF calibration status=%d\n",tofStatus);
   printf("TPC calibration status=%d\n",tpcStatus);
   PrintDetectorStatus();
   return;
}




void PrintDetectorStatus(){
  //
  // GetStatus for the detector which did not implement GetStatus function: 
  //
  // AliTOFAnalysisTaskCalibPass0 calibTask;            // GetStatus implemented
  // AliTRDPreprocessorOffline procestrd;               // GetStatus implemented 
  // AliTPCPreprocessorOffline proces;                  // GetStatus not implemented (next release)
  //                                                    // Logical or of the Validation function used instead
  //
  // AliMeanVertexPreprocessorOffline procesMeanVtx;    // GetStatus not implemented 
  // AliT0PreprocessorOffline procesT0;                 // GetStatus not implemented 
  // SDDcalib                                           // Not automatic update - Not needed 
  //
  // CODE TO BE WRITTEN HERE FOR DETECTOT WITHOUT STATUS 
}
