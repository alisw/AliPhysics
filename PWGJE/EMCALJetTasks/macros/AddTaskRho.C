// $Id$

AliAnalysisTaskRho* AddTaskRho(
   const char    *nJets       = "Jets",
   const char    *nTracks     = "PicoTracks",
   const char    *nClusters   = "CaloClusters",  
   const char    *nRho        = "Rho",
   Double_t       jetradius   = 0.2,
   const char    *cutType     = "TPC",
   Double_t       jetareacut  = 0.01,
   Double_t       emcareacut  = 0,
   TF1           *sfunc       = 0,
   const UInt_t   exclJets    = 2,
   const Bool_t   histo       = kFALSE,
   const char    *taskname    = "Rho"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRho", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRho", "This task requires an input event handler");
    return NULL;
  }


  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s", taskname, nJets,cutType));
  AliAnalysisTaskRho* mgrTask = mgr->GetTask(name.Data());
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRho *rhotask = new AliAnalysisTaskRho(name, histo);
  rhotask->SetExcludeLeadJets(exclJets);
  rhotask->SetScaleFunction(sfunc);
  rhotask->SetOutRhoName(nRho);

  AliParticleContainer *trackCont = rhotask->AddParticleContainer(nTracks);
  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(nClusters);

  AliJetContainer *jetCont = rhotask->AddJetContainer(nJets,cutType,jetradius);
  if (jetCont) {
    jetCont->SetJetAreaCut(jetareacut);
    jetCont->SetAreaEmcCut(emcareacut);
    jetCont->SetJetPtCut(0);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							      TList::Class(),AliAnalysisManager::kOutputContainer,
							      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}
