// $Id$

AliAnalysisTaskSE *AddTaskEMCALTender(
  Bool_t distBC         = kTRUE,   //distance to bad channel
  Bool_t recalibClus    = kTRUE,   //recalibrate cluster energy
  Bool_t recalcClusPos  = kTRUE,   //recalculate cluster position
  Bool_t nonLinearCorr  = kTRUE,   //apply non-linearity
  Bool_t remExoticCell  = kTRUE,   //remove exotic cells
  Bool_t remExoticClus  = kTRUE,   //remove exotic clusters
  Bool_t fidRegion      = kFALSE,  //apply fiducial cuts
  Bool_t calibEnergy    = kTRUE,   //calibrate energy
  Bool_t calibTime      = kTRUE,   //calibrate timing
  Bool_t remBC          = kTRUE,   //remove bad channels
  UInt_t nonLinFunct    = AliEMCALRecoUtils::kBeamTestCorrected,
  Bool_t reclusterize   = kTRUE,   //reclusterize
  Float_t seedthresh    = 0.100,   //seed threshold
  Float_t cellthresh    = 0.050,   //cell threshold
  UInt_t clusterizer    = AliEMCALRecParam::kClusterizerv2,
  Bool_t trackMatch     = kTRUE,   //track matching
  Bool_t updateCellOnly = kFALSE,  //only change if you run your own clusterizer task
  Float_t timeMin       = 100e-9,  //minimum time of physical signal in a cell/digit (s)
  Float_t timeMax       = 900e-9,  //maximum time of physical signal in a cell/digit (s)
  Float_t timeCut       = 900e-9,  //maximum time difference between the digits inside EMC cluster (s)
  const char *pass      = 0        //string defining pass (use none if figured out from path)
) 
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


  gROOT->LoadMacro("$ALICE_ROOT/PWG/EMCAL/macros/ConfigEmcalTenderSupply.C");

  AliEMCALTenderSupply *EMCALSupply = ConfigEmcalTenderSupply(distBC, recalibClus, recalcClusPos, nonLinearCorr, remExoticCell, remExoticClus, 
							      fidRegion, calibEnergy, calibTime, remBC, nonLinFunct, reclusterize, seedthresh, 
							      cellthresh, clusterizer, trackMatch, updateCellOnly, timeMin, timeMax, timeCut);

  if (pass) 
    EMCALSupply->SetPass(pass);

  if (evhand->InheritsFrom("AliESDInputHandler")) {
    AliTender* alitender = mgr->GetTopTasks()->FindObject("AliTender");
    if (!alitender)
      alitender = new  AliTender("AliTender");
    alitender->AddSupply(EMCALSupply);
    alitender->SetDefaultCDBStorage("raw://"); 
    ana = alitender;
    coutput1 = mgr->CreateContainer("emcal_tender_event", 
				    AliESDEvent::Class(), 
				    AliAnalysisManager::kExchangeContainer, 
				    "default_tender");
  } else if (evhand->InheritsFrom("AliAODInputHandler")) {
    AliEmcalTenderTask* emcaltender = mgr->GetTopTasks()->FindObject("AliEmcalTenderTask"); 
    if (!emcaltender)
      emcaltender = new  AliEmcalTenderTask("AliEmcalTenderTask");
    emcaltender->SetEMCALTenderSupply(EMCALSupply);
    ana = emcaltender;
    coutput1 = mgr->CreateContainer("emcal_tender_event", 
				    AliAODEvent::Class(), 
				    AliAnalysisManager::kExchangeContainer, 
				    "default_tender");
  } else {
    ::Error("AddTaskEMCALTender", "Input event handler not recognized, AOD/ESD expected. Returning...");
    return NULL;
  }

  mgr->AddTask(ana);

  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================

  mgr->ConnectInput(ana, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(ana, 1, coutput1 );

  ::Info("AddTaskEMCALTender", "Tender configured");
   
  return ana;
}
