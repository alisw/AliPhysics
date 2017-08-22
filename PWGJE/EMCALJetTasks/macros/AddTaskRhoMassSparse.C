// $Id$

AliAnalysisTaskRhoMassSparse* AddTaskRhoMassSparse(
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
   TF1           *sfunc       = 0,
   const UInt_t   exclJets    = 2,
   const Bool_t   histo       = kFALSE,
   const char    *taskname    = "RhoMass",
   const Bool_t   fRhoCMS      = kTRUE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoMassSparse", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRhoMassSparse", "This task requires an input event handler");
    return NULL;
  }


  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s", taskname, nJetsBkg,cutType));
  AliAnalysisTaskRhoMassSparse* mgrTask = dynamic_cast<AliAnalysisTaskRhoMassSparse *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRhoMassSparse *rhomtask = new AliAnalysisTaskRhoMassSparse(name, histo);
  rhomtask->SetExcludeLeadJets(exclJets);
  rhomtask->SetScaleFunction(sfunc);
  rhomtask->SetOutRhoMassName(nRho);

  AliParticleContainer *trackCont = rhomtask->AddParticleContainer(nTracks);
  AliClusterContainer *clusterCont = rhomtask->AddClusterContainer(nClusters);

  AliJetContainer *bkgJetCont = rhomtask->AddJetContainer(nJetsBkg,cutType,jetradius);
  if (bkgJetCont) {
    bkgJetCont->SetJetAreaCut(jetareacut);
    bkgJetCont->SetAreaEmcCut(emcareacut);
    bkgJetCont->SetJetPtCut(jetptcut);
    bkgJetCont->ConnectParticleContainer(trackCont);
    bkgJetCont->ConnectClusterContainer(clusterCont);
  }

  AliJetContainer *sigJetCont = rhomtask->AddJetContainer(nJetsSig,cutType,jetradius);
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

  rhomtask->SetSmallSystem(fRhoCMS);
  mgr->AddTask(rhomtask);

  // Create containers for input/output
  mgr->ConnectInput(rhomtask, 0, mgr->GetCommonInputContainer());
  if (histo) {
    TString contname(name);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							      TList::Class(),AliAnalysisManager::kOutputContainer,
							      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhomtask, 1, coutput1);
  }

  return rhomtask;
}
