AliAnalysisTaskRho* AddTaskRho (
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
   const char    *suffix      = ""
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

  TString name(Form("AliAnalysisTaskRho_%s_%s", nJets,cutType));
  if (strcmp(suffix,"") != 0) {
    name += "_";
    name += suffix;
  }

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

AliAnalysisTaskRho* AddTaskRho(
   const char    *nJets       = "Jets",
   const char    *nTracks     = "PicoTracks",
   const char    *nClusters   = "CaloClusters",  
   const char    *nRho        = "Rho",
   Double_t       jetradius   = 0.2,
   const char    *cutType     = "TPC",
   Double_t       jetareacut  = 0.01,
   Double_t       emcareacut  = 0,
   const char    *sfuncPath   = "alien:///alice/cern.ch/user/s/saiola/LHC11h_ScaleFactorFunctions.root",
   const char    *sfuncName   = "LHC11h_HadCorr20_ClustersV2",
   const UInt_t   exclJets    = 2,
   const Bool_t   histo       = kFALSE,
   const char    *suffix      = ""
)
{
  TFile *file = TFile::Open(sfuncPath);
  if (!file || file->IsZombie()) {
    ::Error("AddTaskRho", "Could not open scale function file");
    return NULL;
  }

  TF1* sfunc = dynamic_cast<TF1*>(file->Get(sfuncName));

  if (sfunc) {
    ::Info("AddTaskRho", Form("Scale function %s loaded from file %s.", sfuncName, sfuncPath));
  }
  else {
    ::Error("AddTaskRho", Form("Scale function %s not found in file %s.", sfuncName, sfuncPath));
    return NULL;
  }

  AliAnalysisTaskRho *task = AddTaskRho(nJets, nTracks, nClusters, nRho, jetradius, cutType, jetareacut, emcareacut, sfunc, exclJets, histo, suffix);

  file->Close();
  delete file;
  file = 0;

  return task;
}
