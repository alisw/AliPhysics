/// \file AddTaskTPCCalib.C
/// \brief This macros setup the TPC calibration task

Int_t debugLevel  = 2;
Int_t streamLevel = 20;
TTimeStamp startTime(2009,8,7,0,0,0);
TTimeStamp stopTime(2009,12,31,0,0,0);
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
  /// Responsible: Marian Ivanov
  /// Description:
  /// calibCalib is a prefilter
  /// The current OCDB entries transformation are applied on cluster, tracks are refitted

  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCalib *calibCalib = new AliTPCcalibCalib("calibTPC","calibTPC");
  calibCalib->SetDebugLevel(debugLevel);
  calibCalib->SetStreamLevel(streamLevel);
  calibCalib->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibCalib);

}
void AddCalibTimeGain(TObject* task){
  ///  Responsible: Alexander Kalweit
  ///  Description:
  ///  Parameters to set

  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain("calibTimeGain","calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 30*60);
  //calibTimeGain->SetLowMemoryConsumption(kTRUE);
  //calibTimeGain->SetMIP(25.);
  calibTimeGain->SetIsCosmic(kFALSE);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(kFALSE);
  calibTimeGain->SetDebugLevel(debugLevel);
  calibTimeGain->SetStreamLevel(streamLevel);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  myTask->AddJob(calibTimeGain);
}

void AddCalibTime(TObject* task){
  /// Responsible: Dag Larsen
  /// Description:

  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTime *calibTime = new AliTPCcalibTime("calibTime","calibTime",  startTime.GetSec(), stopTime.GetSec(), 20*60);
  calibTime->SetDebugLevel(debugLevel);
  calibTime->SetStreamLevel(streamLevel);
  calibTime->SetTriggerMask(-1,-1,kFALSE);        //accept everything 
  myTask->AddJob(calibTime);
}


void SetupCalibTaskTrain1(TObject* task){
  ///

  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  //AddCalibCalib(task);
  AddCalibTimeGain(task);
  AddCalibTime(task);
  //
  TString path=gSystem->pwd();
  path+=prefix;
  gSystem->mkdir(path);
  myTask->SetDebugOuputhPath(path.Data());

}



void ConfigOCDB(Int_t run){
  ///

  printf("SETUP OCBD for TPC\n");
  printf("SETUP OCBD for TPC\n");
  printf("SETUP OCBD for TPC\n");
  //
  //
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  param->ReadGeoMatrices();
  //  
 
  AliTPCTransform *transform     = AliTPCcalibDB::Instance()->GetTransform() ;
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  //
  transform->SetCurrentRecoParam(tpcRecoParam);
  tpcRecoParam->SetUseGainCorrectionTime(0);
  tpcRecoParam->SetUseRPHICorrection(kTRUE); 
  tpcRecoParam->SetUseTOFCorrection(kFALSE);
  //
  tpcRecoParam->SetUseDriftCorrectionTime(1);
  tpcRecoParam->SetUseDriftCorrectionGY(1);
  //
  tpcRecoParam->SetUseRadialCorrection(kFALSE);
  tpcRecoParam->SetUseQuadrantAlignment(kTRUE);
  //
  tpcRecoParam->SetUseSectorAlignment(kFALSE);
  tpcRecoParam->SetUseGainCorrectionTime(kFALSE);
  tpcRecoParam->SetUseFieldCorrection(kFALSE);
  tpcRecoParam->SetUseExBCorrection(kTRUE);
  AliTPCcalibDB::Instance()->SetRun(run);
}




