AliEmcalClusterMaker* AddTaskEmcalClusterMaker(const UInt_t nonLinFunct   = AliEMCALRecoUtils::kBeamTestCorrected,
                                               const Bool_t remExClus     = kTRUE,
                                               const char *nClusters      = "usedefault",
                                               const char *outClusName    = "EmcCaloClusters",
                                               const Double_t emin        = 0.3,
                                               const Bool_t   histo       = kFALSE,
                                               const char *outputname     = "AnalysisResults.root"
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskHadCorr", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskHadCorr", "This task requires an input event handler");
    return NULL;
  }

  TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"

  TString nCells = "emcalCells";
  if (inputDataType == "ESD")
    nCells = "EMCALCells";

  if (nClusters==0 || strcmp(nClusters,"usedefault")==0) {
    if (inputDataType != "ESD")
      nClusters = "caloClusters";
    else 
      nClusters = "CaloClusters";
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("EmcalClusterMaker_%s", nClusters));
  if (strcmp(outClusName, "") != 0) name += Form("_%s", outClusName);
  AliEmcalClusterMaker *ecm = new AliEmcalClusterMaker(name, histo);
  ecm->SetOutClusName(outClusName);
  ecm->SetCaloCellsName(nCells);
  AliEMCALRecoUtils *ru = new AliEMCALRecoUtils;
  ru->SetNonLinearityFunction(nonLinFunct);
  if(remExClus) ru->SwitchOnRejectExoticCluster();
  ecm->SetRecoUtils(ru);
  AliClusterContainer *clusCont = ecm->AddClusterContainer(nClusters);
  clusCont->SetClusECut(emin);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(ecm);
    
  // Create containers for input/output
  mgr->ConnectInput (ecm, 0, mgr->GetCommonInputContainer());

  if (histo) {
    AliAnalysisDataContainer *coecm = mgr->CreateContainer(name,
							   TList::Class(),
							   AliAnalysisManager::kOutputContainer,
							   outputname);
    mgr->ConnectOutput(ecm,1,coecm);
  }
    
  return ecm;
}
