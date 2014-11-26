AliAnalysisTaskScale* AddTaskScale(
  const char *nTracks        = "Tracks",
  const char *nClusters      = "CaloClustersCorr",
  Double_t    trackptcut     = 0.150,
  Double_t    clusptcut      = 0.150,
  const char *taskname       = "Scale",
  const char *sfuncPath      = 0,
  const char *sfuncName      = 0
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskScale", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskScale", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(Form("%s_%s_%s_%d_%d", taskname, nTracks, nClusters, TMath::FloorNint(trackptcut*1000), TMath::FloorNint(clusptcut*1000)));
  AliAnalysisTaskScale *scaletask = new AliAnalysisTaskScale(name);
  AliParticleContainer *pcont = scaletask->AddParticleContainer(nTracks);
  if(pcont) {
    pcont->SetParticlePtCut(trackptcut);
    pcont->SetParticleEtaLimits(-0.7,0.7); // only accept tracks in the EMCal eta range
  }
  AliClusterContainer  *ccont = scaletask->AddClusterContainer(nClusters);
  if(ccont) ccont->SetClusPtCut(clusptcut);

  if (sfuncPath != 0 && sfuncName != 0) {
    TFile *file = TFile::Open(sfuncPath);
    if (file && !file->IsZombie()) {
      TF1* sfunc = dynamic_cast<TF1*>(file->Get(sfuncName));

      if (sfunc) {
        scaletask->SetScaleFunction(sfunc);
        ::Info("AddTaskScale", Form("Scale function %s loaded from file %s.", sfuncName, sfuncPath));
      }
      else {
        ::Error("AddTaskScale", Form("Scale function %s not found in file %s.", sfuncName, sfuncPath));
      }
      
      file->Close();
      delete file;
      file = 0;
    }
    else {
      ::Warning("AddTaskScale", "Could not open scale function file");
    }
  }

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(scaletask);

  // Create containers for input/output
  TString contname(name);
  contname += "_Histos";
  mgr->ConnectInput (scaletask, 0, mgr->GetCommonInputContainer() );
  AliAnalysisDataContainer *coscale = mgr->CreateContainer(contname,
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectOutput(scaletask,1,coscale);

  return scaletask;
}
