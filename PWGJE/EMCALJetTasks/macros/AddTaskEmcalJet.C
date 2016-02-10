AliEmcalJetTask* AddTaskEmcalJet(
  const char *nTracks                        = "usedefault",
  const char *nClusters                      = "usedefault",
  const Int_t algo                           = 1,                // 1 = AKT, 0 = KT
  const Double_t radius                      = 0.4,
  const AliJetContainer::EJetType_t jetType  = AliJetContainer::kFullJet,
  const Double_t minTrPt                     = 0.15,
  const Double_t minClPt                     = 0.30,
  const Double_t ghostArea                   = 0.005,
  const AliJetContainer::ERecoScheme_t recoSch = AliJetContainer::pt_scheme,
  const char *tag                            = "Jet",
  const Double_t minJetPt                    = 0.,
  const Bool_t selectPhysPrim                = kFALSE,
  const Bool_t lockTask                      = kTRUE,
  const Int_t useExchangeCont                = 0,
  const Bool_t bFillGhosts                   = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJet", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJet", "This task requires an input event handler");
    return 0;
  }
  
  enum EDataType_t {
    kUnknown,
    kESD,
    kAOD
  };

  EDataType_t dataType = kUnknown;

  if (handler->InheritsFrom("AliESDInputHandler")) {
    dataType = kESD;
  }
  else if (handler->InheritsFrom("AliAODInputHandler")) {
    dataType = kAOD;
  }

  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString trackName(nTracks);
  TString clusName(nClusters);

  if (trackName == "usedefault") {
    if (dataType == kESD) {
      trackName = "Tracks";
    }
    else if (dataType == kAOD) {
      trackName = "tracks";
    }
    else {
      trackName = "";
    }
  }

  if (clusName == "usedefault") {
    if (dataType == kESD) {
      clusName = "CaloClusters";
    }
    else if (dataType == kAOD) {
      clusName = "caloClusters";
    }
    else {
      clusName = "";
    }
  }

  AliParticleContainer* partCont = 0;
  if (!trackName.IsNull()) {
    partCont = new AliParticleContainer(trackName);
    partCont->SelectPhysicalPrimaries(selectPhysPrim);
    partCont->SetParticlePtCut(minTrPt);
  }

  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  AliJetContainer::EJetAlgo_t jetAlgo;
  if (algo == 0) {
    jetAlgo = AliJetContainer::kt_algorithm;
  }
  else {
    jetAlgo = AliJetContainer::antikt_algorithm;
  }

  TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, recoSch, radius, partCont, clusCont, tag);

  Printf("Jet task name: %s", name.Data());
 
  AliEmcalJetTask* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask) return mgrTask;  

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name, useExchangeCont);
  jetTask->SetJetType(jetType);
  jetTask->SetJetAlgo(jetAlgo);
  jetTask->SetRecombScheme(recoSch);
  jetTask->SetRadius(radius);
  if (partCont) jetTask->AdoptParticleContainer(partCont);
  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetTracksName(trackName);
  jetTask->SetClusName(clusName);
  jetTask->SetJetsName(tag);
  jetTask->SetMinJetPt(minJetPt);
  jetTask->SetGhostArea(ghostArea);

  if (bFillGhosts) jetTask->SetFillGhost();
  if (lockTask) jetTask->SetLocked();

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(jetTask, 0, cinput);

  TObjArray* cnt = mgr->GetContainers();

  if (useExchangeCont > 0) {
    if (!trackName.IsNull()) {
      AliAnalysisDataContainer* trackCont = static_cast<AliAnalysisDataContainer*>(cnt->FindObject(trackName));
      if (trackCont) {
        mgr->ConnectInput(jetTask, 1, trackCont);
      }
      else {
        ::Error("AddTaskEmcalJet", "Could not find input container '%s'!", nTracks);
      }
    }
  }

  if (useExchangeCont > 1) {
    if (!clusName.IsNull()) {
      AliAnalysisDataContainer* clusCont = static_cast<AliAnalysisDataContainer*>(cnt->FindObject(clusName));
      if (clusCont) {
        mgr->ConnectInput(jetTask, 2, clusCont);
      }
      else {
        ::Error("AddTaskEmcalJet", "Could not find input container '%s'!", nClusters);
      }
    }
  }

  return jetTask;
}
