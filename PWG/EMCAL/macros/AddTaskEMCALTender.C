// $Id$

AliAnalysisTaskSE *AddTaskEMCALTender(
  const char*  p       = "lhc11h",
  Bool_t timeCut       = kFALSE,
  Bool_t distBC        = kTRUE, 
  Bool_t recalibClus   = kTRUE, 
  Bool_t recalibClusPos = kTRUE, 
  Bool_t nonLinearCorr = kTRUE, 
  Bool_t remExotic     = kTRUE,
  Bool_t fidRegion     = kFALSE,
  Bool_t calibEnergy   = kTRUE,
  Bool_t calibTime     = kTRUE,
  Bool_t remBC         = kTRUE,
  Bool_t reclusterize  = kFALSE,
  UInt_t clusterizer   = AliEMCALRecParam::kClusterizerNxN,
  Bool_t trackMatch    = kFALSE,
  Bool_t updateCellOnly= kFALSE,
  const char* pass     = 0)
{
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskEMCALTender", "No analysis manager to connect to.");
    return NULL;
  }

  AliVEventHandler *evhand = mgr->GetInputEventHandler();

  // Create the task and configure it.
  //===========================================================================

  AliAnalysisTaskSE *ana = 0;
  AliAnalysisDataContainer *coutput1 = 0;

  UInt_t nonLinFunct = AliEMCALRecoUtils::kBeamTestCorrected;

  TString period(p);
  period.ToLower();
  if (period == "lhc12a15e")
    nonLinFunct = AliEMCALRecoUtils::kPi0MCv3;
  else if (period == "lhc12a15a")
    nonLinFunct = AliEMCALRecoUtils::kPi0MCv2;

  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/ConfigEmcalTenderSupply.C");

  AliEMCALTenderSupply *EMCALSupply = ConfigEmcalTenderSupply(timeCut, distBC, recalibClus, recalibClusPos, nonLinearCorr, remExotic, 
							      fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, clusterizer, trackMatch, updateCellOnly);
  if (pass) 
    EMCALSupply->SetPass(pass);

  if (evhand->InheritsFrom("AliESDInputHandler")) {
    AliTender* alitender = new  AliTender("AliTender");
    alitender->AddSupply(EMCALSupply);
    alitender->SetDefaultCDBStorage("raw://"); 
    ana = alitender;

    coutput1 = mgr->CreateContainer("emcal_tender_event", 
				    AliESDEvent::Class(), 
				    AliAnalysisManager::kExchangeContainer, 
				    "default_tender");
  }
  else if (evhand->InheritsFrom("AliAODInputHandler")) {
    AliEmcalTenderTask* emcaltender = new  AliEmcalTenderTask("AliEmcalTenderTask");
    emcaltender->SetEMCALTenderSupply(EMCALSupply);
    ana = emcaltender;
    coutput1 = mgr->CreateContainer("emcal_tender_event", 
				    AliAODEvent::Class(), 
				    AliAnalysisManager::kExchangeContainer, 
				    "default_tender");
  }
  else {
    ::Error("AddTaskEMCALTender", "Input event handler not recognized, AOD/ESD expected. Returning...");
    return NULL;
  }

  mgr->AddTask(ana);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  mgr->ConnectInput(ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ana, 1, coutput1 );
   
  return ana;
}
