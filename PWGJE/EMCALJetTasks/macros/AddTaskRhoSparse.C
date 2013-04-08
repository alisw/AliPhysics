// $Id: AddTaskRho_pA.C 58584 2012-09-13 10:37:42Z loizides $

AliAnalysisTaskRhoSparse* AddTaskRhoSparse(
   const char    *nJetsBkg    = "JetsBkg",
   const char    *nJetsSig    = "JetsSig",
   const char    *nTracks     = "PicoTracks",
   const char    *nClusters   = "CaloClusters",  
   const char    *nRho        = "Rho",
   Double_t       jetradius   = 0.2,
   UInt_t         type        = AliAnalysisTaskEmcal::kTPC,
   Double_t       jetareacut  = 0.01,
   Double_t       jetptcut    = 0.0,
   Double_t       emcareacut  = 0,
   TF1           *sfunc       = 0,
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

  TString name(Form("%s_%s_", taskname, nJetsBkg));
  if (type == AliAnalysisTaskEmcal::kTPC) 
    name += "TPC";
  else if (type == AliAnalysisTaskEmcal::kEMCAL) 
    name += "EMCAL";
  else if (type == AliAnalysisTaskEmcal::kUser) 
    name += "USER";
  AliAnalysisTaskRhoSparse *rhotask = new AliAnalysisTaskRhoSparse(name, histo);
  rhotask->SetHistoBins(1000,-0.1,9.9);
  rhotask->SetAnaType(type);
  rhotask->SetScaleFunction(sfunc);
  rhotask->SetJetsName(nJetsBkg);
  rhotask->SetSigJetsName(nJetsSig);
  rhotask->SetTracksName(nTracks);
  rhotask->SetClusName(nClusters);
  rhotask->SetRhoName(nRho);
  rhotask->SetJetAreaCut(jetareacut);
  rhotask->SetAreaEmcCut(emcareacut);
  rhotask->SetJetPtCut(jetptcut);
  rhotask->SetJetRadius(jetradius);
  rhotask->SetExcludeLeadJets(exclJets);
  rhotask->SetRhoCMS(fRhoCMS);

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
