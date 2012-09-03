// $Id$

AliAnalysisTaskRhoBase* AddTaskRhoBase(
   const char    *nJets       = "Jets",
   const char    *nTracks     = "PicoTracks",
   const char    *nClusters   = "CaloClusters",  
   const char    *nRho        = "Rho",
   Double_t       jetradius   = 0.2,
   UInt_t         type        = AliAnalysisTaskEmcal::kTPC,
   Double_t       jetareacut  = 0.01,
   Double_t       emcareacut  = 0,
   Double_t       ptcut       = 0.15,
   TF1           *sfunc       = 0,
   TF1           *rfunc       = 0,
   const Bool_t   histo       = kFALSE,
   const char    *taskname    = "RhoBase"

)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskRhoBase", "No analysis manager to connect to.");
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

  AliAnalysisTaskRhoBase *rhotask = new AliAnalysisTaskRhoBase(taskname,histo);
  rhotask->SetAnaType(type);
  rhotask->SetRhoFunction(rfunc);
  rhotask->SetScaleFunction(sfunc);
  rhotask->SetJetsName(nJets);
  rhotask->SetTracksName(nTracks);
  rhotask->SetClusName(nClusters);
  rhotask->SetRhoName(nRho);
  rhotask->SetJetAreaCut(jetareacut);
  rhotask->SetAreaEmcCut(emcareacut);
  rhotask->SetPtCut(ptcut);
  rhotask->SetJetPtCut(0);
  rhotask->SetJetRadius(jetradius);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(rhotask);

  // Create containers for input/output
  mgr->ConnectInput(rhotask, 0, mgr->GetCommonInputContainer());

  if (histo) {
    TString contname(taskname);
    contname += "_histos";
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							      TList::Class(),AliAnalysisManager::kOutputContainer,
							      Form("%s", AliAnalysisManager::GetCommonFileName()));
    mgr->ConnectOutput(rhotask, 1, coutput1);
  }

  return rhotask;
}
