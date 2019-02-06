//=============================================================================
//
// *** AddTaskTRDCalib
//
// This macros setup the TRD calibration task
//
//=============================================================================

AliAnalysisTask  *AddTaskTRDCalib(Int_t runNumber)
{
  gSystem->Load("libTRDcalib");
  // pointer to the analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskTRDCalib", "No analysis manager to connect to.");
    return NULL;
  }

  // check the input handler
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTask", "This task requires an input event handler");
    return NULL;
  }

  // check the presence of the detectors
  AliCDBEntry* entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject()); 
  if (!grpData) {printf("Failed to get GRP data for run %d",runNumber); return NULL;}
  Int_t activeDetectors = grpData->GetDetectorMask(); 
  TString detStr = AliDAQ::ListOfTriggeredDetectors(activeDetectors);
  TString type = grpData->GetBeamType();
  TString LHCperiod = grpData->GetLHCPeriod();
  Bool_t isLHC10 =  LHCperiod.Contains("LHC10");
  Bool_t isLHC11 =  LHCperiod.Contains("LHC11");
  Bool_t isLHC12 =  LHCperiod.Contains("LHC12");
  Bool_t isLHC13 =  LHCperiod.Contains("LHC13");
  printf("TRD add macro, LHCperiod:%s\n isLHC10:%d isLHC11:%d isLHC12:%d isLHC13:%d\n",LHCperiod.Data(),(Int_t)isLHC10,(Int_t)isLHC11,(Int_t)isLHC12,(Int_t)isLHC13);
  
  ////////////////////////////////////////////
  // Number of timebins
  ///////////////////////////////////////////
  AliTRDcalibDB *calib = AliTRDcalibDB::Instance();
  Int_t nbOfTimeBins = calib->GetNumberOfTimeBinsDCS();
  if(nbOfTimeBins < 0) nbOfTimeBins = 27;
  ////////////////////////////////////////////
  //
  /////////////////////////////////////////////
  Int_t versiongain, subversiongain, versionvdrift, subversionvdrift;


  /////////////////////////
  // The TRD calib Task
  /////////////////////////
  AliTRDCalibTask *calibTask = new AliTRDCalibTask();

  // Disabling TRD CPAss0 as per https://savannah.cern.ch/bugs/?88813
  //calibTask->SetMaxEvent(-1);
  
  if (strstr(type.Data(),"A-A")) {
    calibTask->SetMaxNbTracks(4000);
  }
  else calibTask->SetMaxNbTracks(999999999);
  if(!isLHC10) calibTask->SetRejectPileUpWithTOFOrITS(kTRUE);
  calibTask->SetHisto2d(kTRUE);
  calibTask->SetVector2d(kFALSE);
  calibTask->SetVdriftLinear(kTRUE);
  calibTask->SetExbAlt(kFALSE);
  calibTask->SetNz(0,0);
  calibTask->SetNrphi(0,0);
  calibTask->SetNz(0,1);
  calibTask->SetNrphi(0,1);
  calibTask->SetNz(0,2);
  calibTask->SetNrphi(0,2);
  calibTask->SetLow(0);
  calibTask->SetHigh(30);
  calibTask->SetFillZero(kFALSE);
  // now
  calibTask->AddSelectedTriggerClass("C0OB0-ABCE-NOPF-ALL");
  calibTask->AddSelectedTriggerClass("CTRDCO2-ABCE-NOPF-CENT");
  calibTask->AddSelectedTriggerClass("CTRDCO2-ABCE-NOPF-TRD");
  calibTask->AddSelectedTriggerClass("CTRDCO2-ABCE-NOPF-ALL");
  calibTask->SetReject(kTRUE);
  // before
  //calibTask->AddSelectedTriggerClass("CINT1B-ABCE-NOPF-ALL");
  //calibTask->AddSelectedTriggerClass("CINT1WU-B-NOPF-ALL");
  //calibTask->AddSelectedTriggerClass("CINT7WU-B-NOPF-ALL");
  //calibTask->AddSelectedTriggerClass("CINT7WU-I-NOPF-ALL");
  //calibTask->SetReject(kFALSE);
  //calibTask->SetDebug(2);
  calibTask->SetNbTimeBins(nbOfTimeBins);
  calibTask->SetNumberBinCharge(100);
  //calibTask->SetMaxEvent(10);
  //calibTask->SetThresholdP(1.0);
  calibTask->SetRequirePrimaryVertex(kTRUE);
  calibTask->SetMinNbOfContributors(1);
  if((!isLHC10) && (!isLHC11) && (!isLHC12) && (!isLHC13)) {
    printf("RunII \n");
    calibTask->SetMaxCluster(-1.0);
    calibTask->SetNbMaxCluster(0);
  }
  else {
    printf("RunI \n");
    calibTask->SetMaxCluster(100.0);
    calibTask->SetNbMaxCluster(2);
  }

  if ( detStr.Contains("ITSSPD") && (!detStr.Contains("ITSSDD") || !detStr.Contains("ITSSSD"))) calibTask->SetUseSPDVertex();

  //calibTask->SetLimitChargeIntegration(kTRUE);


  /////////////////////////////
  // Track cuts
  /////////////////////////////
  AliESDtrackCuts *trackCuts = new AliESDtrackCuts("trackcuts","trackcuts");
  trackCuts->SetMinNClustersTPC(50);
  trackCuts->SetMaxChi2PerClusterTPC(3.5);
  //trackCuts->SetMaxCovDiagonalElements(2,2,0.5,0.5,2);
  trackCuts->SetRequireTPCRefit(kTRUE);
  //  trackCuts->SetRequireITSRefit(kTRUE); //removed by BD after discussion with Raphaelle, was introduced to handle pile-up
  trackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  //trackCuts->SetMinNsigmaToVertex(10);
  trackCuts->SetRequireSigmaToVertex(kFALSE);
  trackCuts->SetAcceptKinkDaughters(kFALSE);
  trackCuts->SetMaxDCAToVertexZ(30.0);
  trackCuts->SetMaxDCAToVertexXY(3.0);
  trackCuts->SetDCAToVertex2D(kFALSE);

  calibTask->SetESDtrackCuts(trackCuts);

  mgr->AddTask(calibTask);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  //AliAnalysisDataContainer *cinput = mgr->GetCommonOutputContainer();

  if (!cinput) cinput = mgr->CreateContainer("cchain",TChain::Class(),
                                      AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutput =mgr->CreateContainer("TRDCalib",TList::Class(), AliAnalysisManager::kOutputContainer, "CalibObjects.root");


  mgr->ConnectInput(calibTask,0,cinput);
  mgr->ConnectOutput(calibTask,1,coutput);
  return calibTask;

}
