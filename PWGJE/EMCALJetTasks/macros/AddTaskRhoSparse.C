// $Id$

AliAnalysisTaskRhoSparse* AddTaskRhoSparse(
					   const char    *nJetsBkg    = "JetsBkg",
					   const char    *nJetsSig    = "JetsSig",
					   const char    *nTracks     = "PicoTracks",
					   const char    *nClusters   = "CaloClusters",  
					   const char    *nRho        = "Rho",
					   Double_t       jetradius   = 0.2,
					   const char    *cutType     = "TPC",
					   Double_t       jetareacut  = 0.01,
					   Double_t       jetptcut    = 0.0,
					   Double_t       emcareacut  = 0,
					   TF1           *sfunc       = 0x0,
					   const UInt_t   exclJets    = 2,
					   const Bool_t   histo       = kFALSE,
					   const char    *taskname    = "Rho",
					   const Bool_t   fRhoCMS      = kTRUE
					   )
{  

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoSparse", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRhoSparse", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s", taskname, nJetsBkg,cutType));
  AliAnalysisTaskRhoSparse* mgrTask = static_cast<AliAnalysisTaskRhoSparse *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRhoSparse *rhotask = new AliAnalysisTaskRhoSparse(name, histo);
  rhotask->SetHistoBins(1000,-0.1,9.9);
  rhotask->SetRhoCMS(fRhoCMS);
  rhotask->SetExcludeLeadJets(exclJets);
  rhotask->SetScaleFunction(sfunc);
  rhotask->SetOutRhoName(nRho);

  AliParticleContainer *trackCont = rhotask->AddParticleContainer(nTracks);
  AliClusterContainer *clusterCont = rhotask->AddClusterContainer(nClusters);

  AliJetContainer *bkgJetCont = rhotask->AddJetContainer(nJetsBkg,cutType,jetradius);
  if (bkgJetCont) {
    bkgJetCont->SetJetAreaCut(jetareacut);
    bkgJetCont->SetAreaEmcCut(emcareacut);
    bkgJetCont->SetJetPtCut(0.);
    bkgJetCont->ConnectParticleContainer(trackCont);
    bkgJetCont->ConnectClusterContainer(clusterCont);
  }

  AliJetContainer *sigJetCont = rhotask->AddJetContainer(nJetsSig,cutType,jetradius);
  if (sigJetCont) {
    sigJetCont->SetJetAreaCut(jetareacut);
    sigJetCont->SetAreaEmcCut(emcareacut);
    sigJetCont->SetJetPtCut(jetptcut);
    sigJetCont->ConnectParticleContainer(trackCont);
    sigJetCont->ConnectClusterContainer(clusterCont);
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  rhotask->SetSmallSystem(fRhoCMS);
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
