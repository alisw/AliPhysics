/// \file AddTaskEmcalJetToyModel.C
/// \brief AddTask macro for jet finder used for toy models
///
/// Modification of commonly used AddTaskEmcalJet.C macro used for toy model jets analysis.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \author Leticia Conqueiro Mendez <leticia.cunqueiro.mendez@cern.ch>
/// \date Jul 22, 2016

AliEmcalJetTask* AddTaskEmcalJetToyModel(
  const char *nTracks                        = "usedefault",
  const char *nClusters                      = "usedefault",
  const AliJetContainer::EJetAlgo_t jetAlgo  = AliJetContainer::antikt_algorithm,
  const Double_t radius                      = 0.4,
  const AliJetContainer::EJetType_t jetType  = AliJetContainer::kFullJet,
  const Double_t minTrPt                     = 0.15,
  const Double_t minClPt                     = 0.30,
  const Double_t ghostArea                   = 0.005,
  const AliJetContainer::ERecoScheme_t reco  = AliJetContainer::pt_scheme,
  const char *tag                            = "Jet",
  const Double_t minJetPt                    = 0.,
  const Bool_t lockTask                      = kTRUE,
  const Bool_t bFillGhosts                   = kFALSE,
  const Double_t ptScale                     = 1.,
  const Bool_t kOnTheFly                     =kTRUE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalJetToyModel", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetToyModel", "This task requires an input event handler");
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

    if(kOnTheFly==kFALSE){
   AliTrackContainerToyModel* trackCont = new AliTrackContainerToyModel(trackName);
  trackCont->SetTrackScalePt(ptScale);
  trackCont->SetParticlePtCut(minTrPt);}

 if(kOnTheFly==kTRUE){
  AliMCParticleContainerToyModel* trackCont = new AliMCParticleContainerToyModel(trackName);
  trackCont->SetTrackScalePt(ptScale);
  trackCont->SetParticlePtCut(minTrPt);
  }



  AliClusterContainer* clusCont = 0;
  if (!clusName.IsNull()) {
    clusCont = new AliClusterContainer(clusName);
    clusCont->SetClusECut(0.);
    clusCont->SetClusPtCut(0.);
    clusCont->SetClusHadCorrEnergyCut(minClPt);
    clusCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }

  TString tagStr = TString::Format("%sPtScale%04d", tag, TMath::FloorNint(ptScale*1000));

  TString name = AliJetContainer::GenerateJetName(jetType, jetAlgo, reco, radius, trackCont, clusCont, tagStr);

  Printf("Jet task name: %s", name.Data());
 
  AliEmcalJetTask* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask) return mgrTask;  

  AliEmcalJetTask* jetTask = new AliEmcalJetTask(name);
  jetTask->SetJetType(jetType);
  jetTask->SetJetAlgo(jetAlgo);
  jetTask->SetRecombScheme(reco);
  jetTask->SetRadius(radius);
  if (trackCont) jetTask->AdoptParticleContainer(trackCont);
  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetJetsName(tagStr);
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

  return jetTask;
}
