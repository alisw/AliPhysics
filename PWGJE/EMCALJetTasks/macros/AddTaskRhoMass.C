// $Id$

AliAnalysisTaskRhoMass* AddTaskRhoMass(
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
   const char    *taskname    = "RhoMass"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoMass", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskRhoMass", "This task requires an input event handler");
    return NULL;
  }


  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s", taskname, nJets,cutType));
  AliAnalysisTaskRhoMass* mgrTask = dynamic_cast<AliAnalysisTaskRhoMass*>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  AliAnalysisTaskRhoMass *rhomtask = new AliAnalysisTaskRhoMass(name, histo);
  rhomtask->SetExcludeLeadJets(exclJets);
  rhomtask->SetScaleFunction(sfunc);
  rhomtask->SetOutRhoMassName(nRho);

  AliParticleContainer *trackCont = rhomtask->AddParticleContainer(nTracks);
  AliClusterContainer *clusterCont = rhomtask->AddClusterContainer(nClusters);

  AliJetContainer *jetCont = rhomtask->AddJetContainer(nJets,cutType,jetradius);
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
