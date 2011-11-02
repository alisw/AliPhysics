/*
  macro to extract the OCDB entries

  input: CalibObjects.root
  ouput: TimeGain and TimeVdrift calibration objects for TPC and TRD

  Example:
  .L $ALICE_ROOT/PWG1/CalibMacros/Pass0/makeOCDB.C
  makeOCDB("105160");

*/

const AliTRDCalDet *GetCalDetGain(Int_t runNumber, Int_t version, Int_t subversion);
const AliTRDCalDet *GetCalDetVdrift(Int_t runNumber, Int_t version, Int_t subversion);

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
  procestrd.Init("CalibObjects.root");
  Int_t versionVdriftUsed = procestrd.GetVersionVdriftUsed();
  Int_t subversionVdriftUsed = procestrd.GetSubVersionVdriftUsed();
  Int_t versionGainUsed = procestrd.GetVersionGainUsed();
  Int_t subversionGainUsed = procestrd.GetSubVersionGainUsed();
  if((versionVdriftUsed != 0) && (versionGainUsed != 0)) {
    
    AliTRDCalDet *caldetVdrift =GetCalDetVdrift(runNumber,versionVdriftUsed,subversionVdriftUsed);
    procestrd.SetCalDetVdrift(caldetVdrift);
    AliTRDCalDet *caldetGain =GetCalDetGain(runNumber,versionGainUsed,subversionGainUsed);
    procestrd.SetCalDetGain(caldetGain);
    
    if(caldetVdrift && caldetGain) {
      
      procestrd.SetMinStatsVdriftT0PH(600*10);
      procestrd.SetMinStatsVdriftLinear(50);
      procestrd.SetMinStatsGain(600);
       
      procestrd.CalibVdriftT0("CalibObjects.root",runNumber,runNumber,ocdbStorage);
      procestrd.CalibGain("CalibObjects.root",runNumber,runNumber,ocdbStorage);
      procestrd.CalibChamberStatus(runNumber,runNumber,ocdbStorage);
    }
  }
  
  //Mean Vertex
  AliMeanVertexPreprocessorOffline procesMeanVtx;
  procesMeanVtx.ProcessOutput("CalibObjects.root", ocdbStorage, runNumber);
	
  return;
}

const AliTRDCalDet *GetCalDetVdrift(Int_t runNumber, Int_t version, Int_t subversion){
  //
  // Get Cal Det used during reconstruction for vdrift
  //


  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberVdrift",runNumber, version, subversion);
  if(!entry) {
    printf("Found no entry\n");
    return 0x0;
  }
  const AliCDBId id = entry->GetId();
  version = id.GetVersion();
  subversion = id.GetSubVersion();
  //printf("Found version %d and subversion %d for vdrift\n",version,subversion);
  const AliTRDCalDet* calDet = (AliTRDCalDet *)entry->GetObject();

  return calDet;

}
const AliTRDCalDet *GetCalDetGain(Int_t runNumber, Int_t version, Int_t subversion){
  //
  // Get Cal Det used during reconstruction for vdrift
  //


  AliCDBEntry *entry = AliCDBManager::Instance()->Get("TRD/Calib/ChamberGainFactor",runNumber, version, subversion);
  if(!entry) {
    printf("Found no entry\n");
    return 0x0;
  }
  const AliCDBId id = entry->GetId();
  version = id.GetVersion();
  subversion = id.GetSubVersion();
  //printf("Found version %d and subversion %d for vdrift\n",version,subversion);
  const AliTRDCalDet* calDet = (AliTRDCalDet *)entry->GetObject();

  return calDet;

}
