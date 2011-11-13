/*

 This macros setup the TPC calibration task AddTaskTPCCalib
 for Pass0.
 - the run number is required to config TPC OCDB
 
 The following calibration components are added to the AliTPCAnalysisTaskcalib task:
 1. AliTPCcalibCalib - redo reconstruction with current calibration
 2. AliTPCcalibTimeGain - TPC time dependent gain calibration
 3. AliTPCcalibTime - TPC time dependent drift time calibration

*/

// function to set TPC OCDB parameters
void ConfigOCDB(Int_t crun);

//_____________________________________________________________________________
AliAnalysisTask  *AddTaskTPCCalib(Int_t runNumber)
{
  //
  // add calibration task
  //
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskTPCCalib", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskTPCCalib", "This task requires an input event handler");
    return NULL;
  }  

  // set TPC OCDB parameters
  ConfigOCDB(runNumber);

  // setup task
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  SetupCalibTaskTrain1(task1);
  mgr->AddTask(task1);

  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("TPCCalib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");  

  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  return task1;
}

//_____________________________________________________________________________
void AddCalibCalib(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // calibCalib is a prefilter 
  // The current OCDB entries transformation are applied on cluster, tracks are refitted
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  calibCalib->SetDebugLevel(0);
  calibCalib->SetStreamLevel(0);
  calibCalib->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibCalib);
}

//_____________________________________________________________________________
void AddCalibTimeGain(TObject* task, Bool_t isCosmic = kFALSE, char * name = "calibTimeGain"){
  //
  //  Responsible: Alexander Kalweit
  //  Description: Time Gain calibration
  //

  // Set run time ranges (time stamps)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  if(!entry) { 
    ::Error("AddCalibTimeGain","Cannot get AliCDBEntry");
    return;
  }
  const AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if(!grpData) { 
    ::Error("AddCalibTimeGain","Cannot get AliGRPObject");
    return;
  }
  time_t sTime = grpData->GetTimeStart(); 
  time_t eTime = grpData->GetTimeEnd(); 
  TTimeStamp startRunTime(sTime);
  TTimeStamp stopRunTime(eTime);

  UInt_t year;
  startRunTime.GetDate(kTRUE,0,&year);
  TTimeStamp startTime(year,1,1,0,0,0);
  TTimeStamp stopTime(year,12,31,23,59,59);

  // 
  // setup calibration component
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain(name,"calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 10*60);
  calibTimeGain->SetIsCosmic(isCosmic);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(kTRUE);
  calibTimeGain->SetDebugLevel(0);
  calibTimeGain->SetStreamLevel(0);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibTimeGain->SetLowerTrunc(0.02);
  calibTimeGain->SetUpperTrunc(0.6);

  myTask->AddJob(calibTimeGain);

  AliTPCcalibGainMult *calibGainMult = new AliTPCcalibGainMult("calibGainMult","calibGainMult");
  calibGainMult->SetUseMax(kTRUE);
  calibGainMult->SetDebugLevel(0);
  calibGainMult->SetStreamLevel(0);
  calibGainMult->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibGainMult->SetLowerTrunc(0.02);
  calibGainMult->SetUpperTrunc(0.6);

  myTask->AddJob(calibGainMult);

}

//_____________________________________________________________________________
void AddCalibTime(TObject* task){
  //
  // Responsible: Dag Larsen
  // Description: Time V drift calibration
  //

  // Set run time ranges (time stamps)
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  if(!entry) { 
    ::Error("AddCalibTime","Cannot get AliCDBEntry");
    return;
  }
  const AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if(!grpData) { 
    ::Error("AddCalibTime","Cannot get AliGRPObject");
    return;
  }
  time_t sTime = grpData->GetTimeStart(); 
  time_t eTime = grpData->GetTimeEnd(); 

  TTimeStamp startRunTime(sTime);
  TTimeStamp stopRunTime(eTime);

  UInt_t year;
  startRunTime.GetDate(kTRUE,0,&year);
  TTimeStamp startTime(year,1,1,0,0,0);
  TTimeStamp stopTime(year,12,31,23,59,59);

  // 
  // setup calibration component
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 10*60, 2);
  calibTime->SetDebugLevel(0);
  calibTime->SetStreamLevel(0);
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 

  // max 200 tracks per event
  calibTime->SetCutTracks(200);

  myTask->AddJob(calibTime);
}

//_____________________________________________________________________________
void SetupCalibTaskTrain1(TObject* task){
  //
  // Setup tasks for calibration train
  //
  AddCalibCalib(task);
  AddCalibTimeGain(task);
  AddCalibTime(task);
}

//_____________________________________________________________________________
void ConfigOCDB(Int_t run){
  //
  // Configure TPC OCDB
  //
  printf("SETUP OCBD for TPC\n");
  printf("SETUP OCBD for TPC\n");
  printf("SETUP OCBD for TPC Run =%d\n", run);
  //
  //
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  param->ReadGeoMatrices();
  //
  AliMagF* magF= TGeoGlobalMagField::Instance()->GetField();
  printf("\n\nSET EXB FIELD\t\n\n");
  AliTPCcalibDB::Instance()->SetExBField(magF);
  //
  //
  //
  AliTPCTransform *transform     = AliTPCcalibDB::Instance()->GetTransform() ;
  //
  //AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  //
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("TPC/Calib/RecoParam");
  if (!entry){
    printf("TPC reco param not available");
  }
  TObjArray * array = (TObjArray*)entry->GetObject();
  if (!array){
    printf("TPC reco param not available");
  }
  // 0 - Low Flux (pp), 1- High Flux (Pb-Pb)
  AliTPCRecoParam * tpcRecoParam = (AliTPCRecoParam*)array->At(1);

  transform->SetCurrentRecoParam(tpcRecoParam);
  tpcRecoParam->SetUseGainCorrectionTime(0);
  tpcRecoParam->SetUseRPHICorrection(kFALSE); 
  tpcRecoParam->SetUseTOFCorrection(kFALSE);
  //
  tpcRecoParam->SetUseDriftCorrectionTime(0);
  tpcRecoParam->SetUseDriftCorrectionGY(0);
  //
  tpcRecoParam->SetUseRadialCorrection(kFALSE);
  tpcRecoParam->SetUseQuadrantAlignment(kFALSE);
  //
  tpcRecoParam->SetUseSectorAlignment(kFALSE);
  tpcRecoParam->SetUseFieldCorrection(kFALSE);
  tpcRecoParam->SetUseExBCorrection(kFALSE);
  //
  tpcRecoParam->SetUseMultiplicityCorrectionDedx(kFALSE);
  tpcRecoParam->SetUseAlignmentTime(kFALSE);
  tpcRecoParam->SetUseComposedCorrection(kTRUE);

  AliTPCcalibDB::Instance()->SetRun(run); 
}
