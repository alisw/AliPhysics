/// \file AddTaskDmesonJets.C
/// \brief AddTask macro for the AliAnalysisTaskDmesonJets class.
///
/// AddTask macro for the AliAnalysisTaskDmesonJets class.
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 12, 2016

AliAnalysisTaskDmesonJets* AddTaskDmesonJets(
    const char *ntracks    = "usedefault",
    const char *nclusters  = "usedefault",
    const char *nMCpart    = "",
    Int_t       nMaxTrees  = 2,
    const char *suffix     = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskDmesonJets", "No analysis manager to connect to.");
    return 0;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalJetSpectraQA", "This task requires an input event handler");
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

  TString trackName(ntracks);
  TString clusName(nclusters);
  TString mcPartName(nMCpart);

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

  if (mcPartName == "usedefault") {
    mcPartName = "mcparticles"; // Always needs AliAODMCParticle objects
  }

  TString name("AliAnalysisTaskDmesonJets");
  if (strcmp(suffix, "") != 0) {
    name += TString::Format("_%s", suffix);
  }

  AliAnalysisTaskDmesonJets* jetTask = new AliAnalysisTaskDmesonJets(name, nMaxTrees);
  jetTask->SetVzRange(-10,10);

  if (!trackName.IsNull()) {
    AliHFTrackContainer* trackCont = new AliHFTrackContainer(trackName);
    jetTask->AdoptParticleContainer(trackCont);
  }

  if (!mcPartName.IsNull()) {
    AliMCParticleContainer* partCont = new AliHFAODMCParticleContainer(mcPartName);
    partCont->SetEtaLimits(-1.5, 1.5);
    jetTask->AdoptParticleContainer(partCont);
  }

  jetTask->AddClusterContainer(clusName);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
  
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer* cinput1 = mgr->GetCommonInputContainer();
  TString contname1(name);
  contname1 += "_histos";
  AliAnalysisDataContainer* coutput1 = mgr->CreateContainer(contname1.Data(), 
      TList::Class(), AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(jetTask, 0, cinput1);
  mgr->ConnectOutput(jetTask, 1, coutput1);

  for (Int_t i = 0; i < nMaxTrees; i++) {
    TString contname = TString::Format("%s_tree_%d", name.Data(), i);
    AliAnalysisDataContainer *coutput = mgr->CreateContainer(contname.Data(),
        TTree::Class(),AliAnalysisManager::kOutputContainer,
        Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(jetTask, 2+i, coutput);
  }
  return jetTask;
}
