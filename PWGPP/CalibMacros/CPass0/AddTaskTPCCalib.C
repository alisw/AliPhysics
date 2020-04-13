/*

 This macros setup the TPC calibration task AddTaskTPCCalib
 for CPass0.
 - the run number is required to config TPC OCDB
 
 The following calibration components are added to the AliTPCAnalysisTaskcalib task:
 1. AliTPCcalibCalib - redo reconstruction with current calibration
 2. AliTPCcalibTimeGain - TPC time dependent gain calibration
 3. AliTPCcalibTime - TPC time dependent drift time calibration

*/

// function to set TPC OCDB parameters
int ConfigOCDB();

// functions to set up the trains
void SetupCalibTaskTrain1(TObject* task, const char* options="ALL");
void SetupCalibTaskTrainAlign(TObject* task, const char* options="ALL");
void SetupCalibTaskTrainCluster(TObject* task, const char* options="ALL");
AliAnalysisTask  *AddTaskTPCCalib(const char* options="ALL");

Int_t debugLevel=0;
Int_t streamLevel=0;
Float_t lowtrunc = 0.02;
Float_t hightrunc = 0.6;

//_____________________________________________________________________________
Bool_t isOptionSelected(const char* optionstr, const char* optionsstr, const char* allstr="ALL")
{
  //check if option is enabled and not disabled in options
  TString option=optionstr;
  TString options=optionsstr;
  TString all=allstr;
  TString notOption="-";
  notOption+=option;
  if (options.Contains(notOption)) return kFALSE;
  return (options.Contains(option) || options.Contains(all));
}

//_____________________________________________________________________________
AliAnalysisTask  *AddTaskTPCCalib(const char* options)
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
  if (ConfigOCDB())
  {
    return(NULL);
  }

  AliTPCAnalysisTaskcalib *task=NULL;
  
  // setup task TPCCalib
  if (isOptionSelected("TPCCalib",options))
  {
    ::Info("AddTaskTPCCalib", "Adding TPCCalib");
    TString outputFileName=mgr->GetCommonFileName();
    task=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
    SetupCalibTaskTrain1(task,options);
    mgr->AddTask(task);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
        AliAnalysisManager::kInputContainer);
    for (Int_t i=0; i<task->GetJobs()->GetEntries(); i++) {
      if (task->GetJobs()->At(i)) {
        AliAnalysisDataContainer* coutput = mgr->CreateContainer(task->GetJobs()->At(i)->GetName(),
            AliTPCcalibBase::Class(), 
            AliAnalysisManager::kOutputContainer, 
            "CalibObjects.root:TPCCalib"); 
        mgr->ConnectOutput(task,i,coutput);
      }
    }
    mgr->ConnectInput(task,0,cinput1);
  }
  //
  // setup task TPCAlign
  if (isOptionSelected("TPCAlign",options))
  {
    ::Info("AddTaskTPCCalib", "Adding TPCAlign");
    task=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
    SetupCalibTaskTrainAlign(task,options);
    mgr->AddTask(task);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
        AliAnalysisManager::kInputContainer);
    for (Int_t i=0; i<task->GetJobs()->GetEntries(); i++) {
      if (task->GetJobs()->At(i)) {
        AliAnalysisDataContainer* coutput = mgr->CreateContainer(task->GetJobs()->At(i)->GetName(),
            AliTPCcalibBase::Class(), 
            AliAnalysisManager::kOutputContainer, 
            "CalibObjects.root:TPCAlign"); 
        mgr->ConnectOutput(task,i,coutput);
      }
    }
    mgr->ConnectInput(task,0,cinput1);
  }
  //
  // setup task TPCCluster
  if (isOptionSelected("TPCCluster",options))
  {
    ::Info("AddTaskTPCCalib", "Adding TPCCluster");
    task=new AliTPCAnalysisTaskcalib("CalibObjectsTrain1");
    SetupCalibTaskTrainCluster(task,options);
    mgr->AddTask(task);
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    if (!cinput1) cinput1 = mgr->CreateContainer("cchain",TChain::Class(), 
        AliAnalysisManager::kInputContainer);
    for (Int_t i=0; i<task->GetJobs()->GetEntries(); i++) {
      if (task->GetJobs()->At(i)) {
        AliAnalysisDataContainer* coutput = mgr->CreateContainer(task->GetJobs()->At(i)->GetName(),
            AliTPCcalibBase::Class(), 
            AliAnalysisManager::kOutputContainer, 
            "CalibObjects.root:TPCCluster"); 
        mgr->ConnectOutput(task,i,coutput);
      }
    }
    mgr->ConnectInput(task,0,cinput1);
  }
  //

  return task;
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibCalib");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibTimeGain");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibTime");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibTracks");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibAlign");
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
  ::Info("AddTaskTPCCalib", "AddCalibAlignInterpolation");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibLaser");
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
  ::Info("AddTaskTPCCalib", "AddTaskTPCCalib:AddCalibCosmic");
}



//_____________________________________________________________________________
void SetupCalibTaskTrain1(TObject* task, const char* options){
  //
  // Setup tasks for calibration train
  //
  if (isOptionSelected(":CalibCalib",options)) AddCalibCalib(task);
  if (isOptionSelected(":CalibTimeGain",options)) AddCalibTimeGain(task);
  if (isOptionSelected(":CalibTimeDrift",options)) AddCalibTime(task);
  if (isOptionSelected(":CalibAlignInterpolation",options)) AddCalibAlignInterpolation(task);
}

void SetupCalibTaskTrainAlign(TObject* task, const char* options){
  //
  // Setup tasks for calibration train
  //
  //if (isOptionSelected(":CalibAlign",options)) AddCalibAlign(task);
  if (isOptionSelected(":CalibLaser",options)) AddCalibLaser(task);
  //AddCalibCosmic(task);
}

void SetupCalibTaskTrainCluster(TObject* task, const char* options){
  //
  // Setup tasks for calibration train
  //
  if (isOptionSelected(":CalibTracks",options)) AddCalibTracks(task);
}

//_____________________________________________________________________________
int ConfigOCDB(){
  //
  // Configure TPC OCDB
  //
  bool print_info = !(getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0);

  if (print_info) ::Info("AddTaskTPCCalib", "SETUP OCBD\n");
  //
  //
  AliTPCParam *param= AliTPCcalibDB::Instance()->GetParameters();
  param->ReadGeoMatrices();
  //
  AliMagF* magF= (AliMagF*)(TGeoGlobalMagField::Instance()->GetField());
  if (print_info) ::Info("AddTaskTPCCalib", "\n\nSET EXB FIELD\t\n\n");
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
    return(1);
  }
  TObjArray * array = (TObjArray*)entry->GetObject();
  if (!array){
    ::Error("AddTaskTPCCalib","TPC reco param not available");
    return(1);
  }
  
  //get the beam type from OCDB to decide which type of reco param we need -
  //high or low flux
  entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());  // new GRP entry
  
  if ((grpData->GetDetectorMask() & AliDAQ::kTPC) == 0)
  {
    ::Error("AddTaskTPCCalib", "TPC not in list of active detectors for this run");
    return(1);
  }

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
  ::Info("AddTaskTPCCalib", "Beam type: %s, using fluxType=%i",beamType.Data(),fluxType);
  if (print_info) tpcRecoParam->Print();
  lowtrunc = tpcRecoParam->GetMinFraction();
  hightrunc = tpcRecoParam->GetMaxFraction();

  transform->SetCurrentRecoParam(tpcRecoParam);

  // ===| set up gain calibration type |========================================
  //
  // default is Full Calibration in CPass0
  // will be overwritte by recCPass0.sh
  // NOTE: This must be consistent to the settings in mergeMakeOCDB.byComponent.perStage.sh (makeOCDB.C)
  //
  AliTPCPreprocessorOffline::EGainCalibType tpcGainCalibType=AliTPCPreprocessorOffline::kFullGainCalib;

  // --- check for overwrites from environment variable
  //
  const TString sGainTypeFromEnv(gSystem->Getenv("TPC_CPass0_GainCalibType"));
  if (!sGainTypeFromEnv.IsNull()) {
    const AliTPCPreprocessorOffline::EGainCalibType tpcGainCalibTypeEnv=AliTPCPreprocessorOffline::GetGainCalibrationTypeFromString(sGainTypeFromEnv);
    if (tpcGainCalibTypeEnv==AliTPCPreprocessorOffline::kNGainCalibTypes) {
      ::Fatal("AddTaskTPCCalib","Could not set up gain calibration type from environment variable TPC_CPass0_GainCalibType: %s",sGainTypeFromEnv.Data());
    }

    ::Info("AddTaskTPCCalib","Setting gain calibration type from environment variable TPC_CPass0_GainCalibType: %d", tpcGainCalibTypeEnv);
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

  // ===| Further specific overwrites for CPass0 |==============================
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
  tpcRecoParam->SetUseAlignmentTime(kFALSE);
  tpcRecoParam->SetUseComposedCorrection(kTRUE);
  
  return(0);
}
