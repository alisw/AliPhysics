//=============================================================================
//
// *** AddTaskTPCCalib
//
// This macros setup the TPC calibration task
//
//=============================================================================




Int_t debugLevel  = 2;
Int_t streamLevel = 20;
TTimeStamp startTime(2010,2,1,0,0,0);
TTimeStamp stopTime(2010,12,31,0,0,0);
char * prefix = "/V6/";

void ConfigOCDB(Int_t crun);

AliAnalysisTask  *AddTaskTPCCalib(Int_t runNumber)
{
  gSystem->Load("libTPCcalib");
  // pointer to the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTPCCalib", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }  
  ConfigOCDB(runNumber);
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  //
  SetupCalibTaskTrain1(task1);
  mgr->AddTask(task1);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput1 =mgr->CreateContainer("TPCCalib",TObjArray::Class(), AliAnalysisManager::kOutputContainer, "AliESDfriends_v1.root");  
  mgr->ConnectInput(task1,0,cinput1);
  mgr->ConnectOutput(task1,0,coutput1);
  return task1;
}




void AddCalibCalib(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // calibCalib is a prefilter 
  // The current OCDB entries transformation are applied on cluster, tracks are refitted
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  calibCalib->SetDebugLevel(debugLevel);
  calibCalib->SetStreamLevel(streamLevel);
  calibCalib->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibCalib);

}
void AddCalibTimeGain(TObject* task, Bool_t isCosmic = kFALSE, char * name = "calibTimeGain"){
  //
  //  Responsible: Alexander Kalweit
  //  Description:
  //  Parameters to set
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain(name,"calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 30*60);
  //calibTimeGain->SetLowMemoryConsumption(kTRUE);
  //calibTimeGain->SetMIP(25.);
  calibTimeGain->SetIsCosmic(isCosmic);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(kFALSE);
  calibTimeGain->SetDebugLevel(debugLevel);
  calibTimeGain->SetStreamLevel(streamLevel);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  myTask->AddJob(calibTimeGain);
}

void AddCalibTime(TObject* task){
  //
  // Responsible: Dag Larsen
  // Description:
  //
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 20*60);
  calibTime->SetDebugLevel(debugLevel);
  calibTime->SetStreamLevel(streamLevel);
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibTime);
}

void AddCalibLaser(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibLaser *calibLaser = new AliTPCcalibLaser("laserTPC","laserTPC");
  calibLaser->SetDebugLevel(debugLevel);
  calibLaser->SetStreamLevel(streamLevel);
  calibLaser->SetTriggerMask(-1,-1,kFALSE);        //accept everything
  myTask->AddJob(calibLaser);
}


void AddCalibAlign(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibAlign *calibAlign = new AliTPCcalibAlign("alignTPC","Alignment of the TPC sectors");
  calibAlign->SetDebugLevel(debugLevel);
  calibAlign->SetStreamLevel(streamLevel);
  calibAlign->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  myTask->AddJob(calibAlign);
}

void AddCalibCosmic(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // Histogram residuals and pulls of the track parameters in bins of track parameters
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  calibCosmic->SetDebugLevel(debugLevel);
  calibCosmic->SetStreamLevel(streamLevel);
  calibCosmic->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  myTask->AddJob(calibCosmic);
}




void SetupCalibTaskTrain1(TObject* task){
  //
  //
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AddCalibCalib(task);
  AddCalibTimeGain(task);
  AddCalibTimeGain(task,kTRUE,"calibTimeGainCosmic"); // 2nd task for cosmic runs
  AddCalibTime(task);
  AddCalibLaser(task);
  AddCalibAlign(task);
  AddCalibCosmic(task);
  //
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  myTask->SetDebugOuputhPath(path.Data());

}



void ConfigOCDB(Int_t run){
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
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  //
  transform->SetCurrentRecoParam(tpcRecoParam);
  tpcRecoParam->SetUseGainCorrectionTime(0);
  tpcRecoParam->SetUseRPHICorrection(kTRUE); 
  tpcRecoParam->SetUseTOFCorrection(kFALSE);
  //
  tpcRecoParam->SetUseDriftCorrectionTime(0);
  tpcRecoParam->SetUseDriftCorrectionGY(0);
  //
  tpcRecoParam->SetUseRadialCorrection(kFALSE);
  tpcRecoParam->SetUseQuadrantAlignment(kFALSE);
  //
  tpcRecoParam->SetUseSectorAlignment(kTRUE);
  tpcRecoParam->SetUseGainCorrectionTime(kFALSE);
  tpcRecoParam->SetUseFieldCorrection(kFALSE);
  tpcRecoParam->SetUseExBCorrection(kTRUE);
  AliTPCcalibDB::Instance()->SetRun(run); 
}




