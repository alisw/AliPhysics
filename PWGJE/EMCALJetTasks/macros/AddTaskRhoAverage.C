// $Id$

AliAnalysisTaskRhoAverage* AddTaskRhoAverage(
   const char    *nJets       = "Jets",
   const char    *nTracks     = "PicoTracks",
   const char    *nClusters   = "CaloClusters",  
   const char    *nRho        = "Rho",
   Double_t       jetradius   = 0.2,
   const char    *cutType     = "TPC",
   Double_t       jetareacut  = 0.01,
   Double_t       emcareacut  = 0,
   Double_t       trackptcut  = 0.15,
   Double_t       clusptcut   = 0.30,
   TF1           *sfunc       = 0,
   const UInt_t   exclPart    = 2,
   const UInt_t   rhotype     = 1,
   const Bool_t   histo       = kFALSE,
   const char    *taskname    = "RhoAverage"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoAverage", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRhoAverage", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s_%s", taskname, nTracks, nClusters, cutType));

  AliAnalysisTaskRhoAverage *rhotask = new AliAnalysisTaskRhoAverage(name, histo);
  rhotask->SetExcludeLeadPart(exclPart);
  rhotask->SetScaleFunction(sfunc);
  rhotask->SetOutRhoName(nRho);
  rhotask->SetRhoType(rhotype);

  AliParticleContainer *trackCont = rhotask->AddParticleContainer(nTracks);
  if (trackCont) trackCont->SetTrackPtCut(trackptcut);

  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(nClusters);
  if (clusterCont) clusterCont->SetClusPtCut(clusptcut);

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
