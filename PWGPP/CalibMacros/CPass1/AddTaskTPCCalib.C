/*

 This macros setup the TPC calibration task AddTaskTPCCalib
 for CPass1.
 - the run number is required to config TPC OCDB
 
 The following calibration components are added to the AliTPCAnalysisTaskcalib task:
 1. AliTPCcalibCalib - redo reconstruction with current calibration
 2. AliTPCcalibTimeGain - TPC time dependent gain calibration
 3. AliTPCcalibTime - TPC time dependent drift time calibration

*/

// function to set TPC OCDB parameters
void ConfigOCDB(Int_t crun);
void SetupCalibTaskTrain1(TObject* task);
void SetupCalibTaskTrainAlign(TObject* task);
void SetupCalibTaskTrainCluster(TObject* task);

Int_t debugLevel=0;
Int_t streamLevel=0;
Float_t lowtrunc = 0.02;
Float_t hightrunc = 0.6;

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

  // setup task TPCCalib
  TString outputFileName=mgr->GetCommonFileName();
  AliTPCAnalysisTaskcalib *task1=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  SetupCalibTaskTrain1(task1);
  mgr->AddTask(task1);
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  for (Int_t i=0; i<task1->GetJobs()->GetEntries(); i++) {
    if (task1->GetJobs()->At(i)) {
      AliAnalysisDataContainer* coutput = mgr->CreateContainer(task1->GetJobs()->At(i)->GetName(),
                                                               AliTPCcalibBase::Class(), 
                                                               AliAnalysisManager::kOutputContainer, 
                                                               "CalibObjects.root:TPCCalib"); 
      mgr->ConnectOutput(task1,i,coutput);
    }
  }
  mgr->ConnectInput(task1,0,cinput1);
  //
  // setup task TPCAlign
  AliTPCAnalysisTaskcalib *taskAlign=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  SetupCalibTaskTrainAlign(taskAlign);
  mgr->AddTask(taskAlign);
  AliAnalysisDataContainer *cinput2 = mgr->GetCommonInputContainer();
  if (!cinput2) cinput2 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  for (Int_t i=0; i<taskAlign->GetJobs()->GetEntries(); i++) {
    if (taskAlign->GetJobs()->At(i)) {
      AliAnalysisDataContainer* coutput = mgr->CreateContainer(taskAlign->GetJobs()->At(i)->GetName(),
                                                               AliTPCcalibBase::Class(), 
                                                               AliAnalysisManager::kOutputContainer, 
                                                               "CalibObjects.root:TPCAlign"); 
      mgr->ConnectOutput(taskAlign,i,coutput);
    }
  }
  mgr->ConnectInput(taskAlign,0,cinput2);
  //
  // setup task TPCCluster
  AliTPCAnalysisTaskcalib *taskCluster=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
  SetupCalibTaskTrainCluster(taskCluster);
  mgr->AddTask(taskCluster);
  AliAnalysisDataContainer *cinput3 = mgr->GetCommonInputContainer();
  if (!cinput3) cinput3 = mgr->CreateContainer("cchain",TChain::Class(), 
                                      AliAnalysisManager::kInputContainer);
  for (Int_t i=0; i<taskCluster->GetJobs()->GetEntries(); i++) {
    if (taskCluster->GetJobs()->At(i)) {
      AliAnalysisDataContainer* coutput = mgr->CreateContainer(taskCluster->GetJobs()->At(i)->GetName(),
                                                               AliTPCcalibBase::Class(), 
                                                               AliAnalysisManager::kOutputContainer, 
                                                               "CalibObjects.root:TPCCluster"); 
      mgr->ConnectOutput(taskCluster,i,coutput);
    }
  }
  mgr->ConnectInput(taskCluster,0,cinput3);
  //

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
void AddCalibTimeGain(TObject* task, Bool_t isCosmic = kFALSE, const char * name = "calibTimeGain"){
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

  Bool_t useQmax = (grpData->GetBeamType()).Contains("Pb-Pb") || (grpData->GetBeamType()).Contains("A-A");

  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibTimeGain *calibTimeGain = new AliTPCcalibTimeGain(name,"calibTimeGain", startTime.GetSec(), stopTime.GetSec(), 10*60);
  calibTimeGain->SetIsCosmic(isCosmic);
  calibTimeGain->SetUseCookAnalytical(kTRUE);
  calibTimeGain->SetUseMax(useQmax);
  calibTimeGain->SetDebugLevel(0);
  calibTimeGain->SetStreamLevel(0);
  calibTimeGain->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibTimeGain->SetLowerTrunc(lowtrunc);
  calibTimeGain->SetUpperTrunc(hightrunc);

  myTask->AddJob(calibTimeGain);

  AliTPCcalibGainMult *calibGainMult = new AliTPCcalibGainMult("calibGainMult","calibGainMult");
  calibGainMult->SetUseMax(useQmax);
  calibGainMult->SetDebugLevel(0);
  calibGainMult->SetStreamLevel(0);
  calibGainMult->SetTriggerMask(-1,-1,kTRUE);        //reject laser
  calibGainMult->SetLowerTrunc(lowtrunc);
  calibGainMult->SetUpperTrunc(hightrunc);

  myTask->AddJob(calibGainMult);

  // ===| get reco param                     |==================================
  AliTPCTransform *transform = AliTPCcalibDB::Instance()->GetTransform() ;
  const AliTPCRecoParam *recoParam = transform->GetCurrentRecoParam();
  const Int_t spec = recoParam->GetEventSpecie();

  Float_t minTPCsignalN = spec&AliRecoParam::kLowMult?100.:90.;

  // ===| settings via environment variables |==================================
  //
  // ---| minimum number of PID clusters     |----------------------------------
  const TString sMinTPCsignalN(gSystem->Getenv("TPC_GainCalib_minSignalN"));
  if (!sMinTPCsignalN.IsNull()) {
    minTPCsignalN=sMinTPCsignalN.Atof();
    ::Info("AddTaskTPCCalib","Setting minium number of PID clusters for gain calibration from environment variable TPC_CPass0_GainCalib_minSignalN: %.0f", minTPCsignalN);
  }

  calibTimeGain->SetMinTPCsignalN(minTPCsignalN);
  calibGainMult->SetMinTPCsignalN(minTPCsignalN);
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

  // max 15000 tracks per event
  calibTime->SetCutTracks(15000);

  myTask->AddJob(calibTime);
}


void AddCalibTracks(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // Histogram residuals and pulls of the track parameters in bins of track parameters
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task; 
  AliTPCClusterParam * clusterParam = AliTPCcalibDB::Instance()->GetClusterParam();

   AliTPCcalibTracksCuts *cuts = new AliTPCcalibTracksCuts(30, 0.4, 5, 0.13, 0.018);
  //
  AliTPCcalibTracks *calibTracks =  new AliTPCcalibTracks("calibTracks", "Resolution calibration object for tracks", clusterParam, cuts);
  calibTracks->SetDebugLevel(debugLevel);
  calibTracks->SetStreamLevel(streamLevel);
  calibTracks->SetTriggerMask(-1,-1,kTRUE);       
  myTask->AddJob(calibTracks); 
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

void AddCalibAlignInterpolation(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task; 
  AliTPCcalibAlignInterpolation *calibAlign = new AliTPCcalibAlignInterpolation("alignTPCInterpolation","Space point distortion calibration using interpolation",0);
  calibAlign->SetDebugLevel(debugLevel);
  calibAlign->SetStreamLevel(streamLevel);
  calibAlign->SetStreamLevelTrack(0);
  calibAlign->SetSyswatchStep(10000);
  //  Uncomenting following 2 lines we can enable extendad version of infromation dump
  //  calibAlign->SetStreamLevelTrack(AliTPCcalibAlignInterpolation::kStremInterpolation);
  //  calibAlign->SetSyswatchStep(10);

  calibAlign->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  ::Info("AddCalibAlignInterpolation", "Trigger mask set to accept everything");
  calibAlign->Dump();
  myTask->AddJob(calibAlign);
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


void AddCalibCosmic(TObject* task){
  //
  // Responsible: Marian Ivanov
  // Description:
  // Histogram residuals and pulls of the track parameters in bins of track parameters
  // Dump cosmic tracks to the tree
  //
  AliTPCAnalysisTaskcalib* myTask = (AliTPCAnalysisTaskcalib*) task;
  AliTPCcalibCosmic *calibCosmic = new AliTPCcalibCosmic("cosmicTPC","cosmicTPC");
  calibCosmic->SetDebugLevel(debugLevel);
  calibCosmic->SetStreamLevel(1);
  calibCosmic->SetTriggerMask(-1,-1,kTRUE);        //accept everything
  myTask->AddJob(calibCosmic);
}



//_____________________________________________________________________________
void SetupCalibTaskTrain1(TObject* task){
  //
  // Setup tasks for calibration train
  //
  // AddCalibCalib(task); - disable refitting
  AddCalibTimeGain(task);
  AddCalibTime(task);
  //AddCalibCosmic(task);
  AddCalibAlignInterpolation(task);
}

void SetupCalibTaskTrainAlign(TObject* task){
  //
  // Setup tasks for calibration train
  //
  AddCalibAlign(task);
  AddCalibLaser(task);
}

void SetupCalibTaskTrainCluster(TObject* task){
  //
  // Setup tasks for calibration train
  //
  AddCalibTracks(task);
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
  AliMagF* magF= (AliMagF*)TGeoGlobalMagField::Instance()->GetField();
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
    ::Error("AddTaskTPCCalib","TPC reco param not available");
    return;
  }
  TObjArray * array = (TObjArray*)entry->GetObject();
  if (!array){
    ::Error("AddTaskTPCCalib","TPC reco param not available");
    return;
  }
  
  //get the beam type from OCDB to decide which type of reco param we need -
  //high or low flux
  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    ::Error("AddTaskTPCCalib","GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }
  // 0 - Low Flux (pp), 1- High Flux (Pb-Pb)
  Int_t fluxType=0;
  if (beamType.Contains("p-p")) {fluxType=0;}
  if (beamType.Contains("Pb-Pb") || beamType.Contains("A-A")) {fluxType=1;}
  AliTPCRecoParam * tpcRecoParam = (AliTPCRecoParam*)array->At(fluxType);
  ::Info("AddTaskTPCCalib","Beam type: %s, using fluxType=%i",beamType.Data(),fluxType);
  tpcRecoParam->Print();
  lowtrunc = tpcRecoParam->GetMinFraction();
  hightrunc = tpcRecoParam->GetMaxFraction();

  transform->SetCurrentRecoParam(tpcRecoParam);

  // ===| set up gain calibration type |========================================
  //
  // default is Combined Calibration + Residual QA in CPass1
  // will be overwritte by recCPass1.sh
  // NOTE: This must be consistent to the settings in mergeMakeOCDB.byComponent.perStage.sh (makeOCDB.C)
  //
  AliTPCPreprocessorOffline::EGainCalibType tpcGainCalibType=AliTPCPreprocessorOffline::kResidualGainQA;

  // --- check for overwrites from environment variable
  //
  const TString sGainTypeFromEnv(gSystem->Getenv("TPC_CPass1_GainCalibType"));
  if (!sGainTypeFromEnv.IsNull()) {
    const AliTPCPreprocessorOffline::EGainCalibType tpcGainCalibTypeEnv=AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString(sGainTypeFromEnv);
    if (tpcGainCalibTypeEnv==AliTPCPreprocessorOffline::kNGainCalibTypes) {
      ::Fatal("AddTaskTPCCalib","Could not set up gain calibration type from environment variable TPC_CPass1_GainCalibType: %s",sGainTypeFromEnv.Data());
    }

    ::Info("AddTaskTPCCalib","Setting gain calibration type from environment variable TPC_CPass1_GainCalibType: %d", tpcGainCalibTypeEnv);
    tpcGainCalibType=tpcGainCalibTypeEnv;
  }

  if (tpcGainCalibType==AliTPCPreprocessorOffline::kFullGainCalib) {
    tpcRecoParam->SetUseGainCorrectionTime(0);
    tpcRecoParam->SetUseRPHICorrection(kFALSE);
    tpcRecoParam->SetUseMultiplicityCorrectionDedx(kFALSE);
    tpcRecoParam->SetCorrectionHVandPTMode(1);
  }

  // TOF correction should always be on
  // tpcRecoParam->SetUseTOFCorrection(kFALSE);

  // ===| For the rest, in CPass1 use a default setting |=======================
  //
// tpcRecoParam->SetUseDriftCorrectionTime(0);
// tpcRecoParam->SetUseDriftCorrectionGY(0);
// //
// tpcRecoParam->SetUseRadialCorrection(kFALSE);
// tpcRecoParam->SetUseQuadrantAlignment(kFALSE);
// //
// tpcRecoParam->SetUseSectorAlignment(kFALSE);
// tpcRecoParam->SetUseFieldCorrection(kFALSE);
// tpcRecoParam->SetUseExBCorrection(kFALSE);
// //
// tpcRecoParam->SetUseAlignmentTime(kFALSE);
// tpcRecoParam->SetUseComposedCorrection(kTRUE);

  // ===| Initialise AliTPCcalibDB |============================================
  //
  AliTPCcalibDB::Instance()->SetRun(run);
}
