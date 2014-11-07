// $Id$

AliJetEmbeddingFromGenTask* AddTaskJetEmbeddingFromGen(
  Int_t           genType        = 0, //use Pythia as default
  Double_t        ptHardMin      = 50.,
  Double_t        ptHardMax      = 1000.,
  Double_t        ecms           = 2760.,
  const char     *tracksName     = "GenParticles",
  const char     *taskName       = "JetEmbeddingFromGenTask",
  const Double_t  minPt          = 10,
  const Double_t  maxPt          = 10,
  const Double_t  minEta         = -0.9,
  const Double_t  maxEta         = 0.9,
  const Double_t  minPhi         = 0,
  const Double_t  maxPhi         = TMath::Pi() * 2,
  const Bool_t    copyArray      = kTRUE,
  const Bool_t    drawQA         = kTRUE,
  const char     *partonInfoName = "PartonInfo"
)
{
  AliGenerator *genGen = NULL;
  if(genType==0) { //PYTHIA Perugia 2011
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenPythia.C");
    genGen = AddMCGenPythia(ecms, ptHardMin, ptHardMax, 2);
  }
  else if(genType==1 || genType==2) { //QPYTHIA and PYQUEN
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddMCGenQuench.C");
    genGen = AddMCGenQuench(ecms, ptHardMin, ptHardMax, genType);
  }
  if(!genGen)   {
    ::Error("AddTaskJetEmbeddingFromGenTask", "Generator does not exist");
    return NULL;
  }

  AliJetEmbeddingFromGenTask *task = AddTaskJetEmbeddingFromGen(genGen,tracksName,taskName,minPt,maxPt,minEta,maxEta,minPhi,maxPhi,copyArray,drawQA,partonInfoName);

  return task;

}

AliJetEmbeddingFromGenTask* AddTaskJetEmbeddingFromGen(
  AliGenerator   *genGen,
  const char     *tracksName     = "GenParticles",
  const char     *taskName       = "JetEmbeddingFromGenTask",
  const Double_t  minPt          = 10,
  const Double_t  maxPt          = 10,
  const Double_t  minEta         = -0.9,
  const Double_t  maxEta         = 0.9,
  const Double_t  minPhi         = 0,
  const Double_t  maxPhi         = TMath::Pi() * 2,
  const Bool_t    copyArray      = kTRUE,
  const Bool_t    drawQA         = kTRUE,
  const char     *partonInfoName = ""
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
  // if pythia is used as a generator, tell it not to print the event history to the screen
  if(genGen) {
    if(TString(genGen->IsA()->GetName()).EqualTo("AliGenPythia")) genGen->AliGenPythia::SetEventListRange(-10, -10);
    jetEmb->SetGen(genGen);
  }
  jetEmb->SetTracksName(tracksName);
  jetEmb->SetPartonInfoName(partonInfoName);
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
