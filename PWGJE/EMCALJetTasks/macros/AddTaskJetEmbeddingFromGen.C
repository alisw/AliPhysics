// $Id$

AliJetEmbeddingFromGenTask* AddTaskJetEmbeddingFromGen(
  AliGenerator   *genGen,
  const char     *tracksName   = "GenParticles",
  const char     *taskName     = "JetEmbeddingFromGenTask",
  const Double_t  minPt        = 10,
  const Double_t  maxPt        = 10,
  const Double_t  minEta       = -0.9,
  const Double_t  maxEta       = 0.9,
  const Double_t  minPhi       = 0,
  const Double_t  maxPhi       = TMath::Pi() * 2,
  const Bool_t    copyArray    = kTRUE,
  const Bool_t    drawQA       = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskJetEmbedding", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskJetEmbedding", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  AliJetEmbeddingFromGenTask *jetEmb = new AliJetEmbeddingFromGenTask(taskName,drawQA);
  jetEmb->SetGen(genGen);
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetEtaRange(minEta, maxEta);
  jetEmb->SetPhiRange(minPhi, maxPhi);
  jetEmb->SetPtRange(minPt, maxPt);
  jetEmb->SetCopyArray(copyArray);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(jetEmb);
    
  // Create containers for input/output
  mgr->ConnectInput (jetEmb, 0, mgr->GetCommonInputContainer() );

  if (drawQA) {
    TString contName = taskName;
    contName += "_histos";
    AliAnalysisDataContainer *outc = mgr->CreateContainer(contName,
                                                          TList::Class(),
                                                          AliAnalysisManager::kOutputContainer,
                                                          "AnalysisResults.root");
    mgr->ConnectOutput(jetEmb, 1, outc);
  }


  return jetEmb;
}
