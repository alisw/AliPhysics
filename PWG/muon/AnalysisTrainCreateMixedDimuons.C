//===========================================================================================

// Analysis Train driving an AliAnalysisTaskCreateMixedDimuons multiple-event task 
// creating mixed muon pairs
//
// Authors Alessandro De Falco and Antonio Uras, INFN Cagliari
// alessandro.de.falco@ca.infn.it  antonio.uras@ca.infn.it

//===========================================================================================

void AnalysisTrainCreateMixedDimuons(Char_t *nameTagDir = ".",
				     Char_t *nameOutFileAOD = "AliAODsMixedEvents.root",
				     Int_t   bufferSize = 2) {
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWGmuon");

  AliMultiEventInputHandler *inputHandler  = new AliMultiEventInputHandler(bufferSize);

  // setting pools definition...

  AliEventPoolMuon *pool = new AliEventPoolMuon("eventPool", "AOD");
  pool->SetTagDirectory(nameTagDir);
  pool->SetMultiplicityRange(1, 100, 100);         // min, max, step
  pool->SetNFWMuonRange(1, 10, 10);                // min, max, step
  pool->SetPrimaryVertexZRange(-100., -90., 1.);   // min, max, step
  pool->Init();
  
  // ... done!

  AliAnalysisManager *mgr  = new AliAnalysisManager("MixingAnalysisTest");
  mgr -> SetInputEventHandler(inputHandler);
  mgr -> SetEventPool(pool);

  inputHandler->SetEventPool(pool);
    
  AliAnalysisTaskCreateMixedDimuons *mixTask = new AliAnalysisTaskCreateMixedDimuons("AliAnalysisTaskCreateMixedDimuons");
  mixTask -> SetDebugLevel(10);
  mgr -> AddTask(mixTask);
  
  AliAnalysisDataContainer *cInputChain  = mgr->GetCommonInputContainer();
  
  AliAnalysisDataContainer *cOutputUserAODTree = mgr->CreateContainer("cOutputUserAODTree", TTree::Class(),    
								      AliAnalysisManager::kOutputContainer, 
								      nameOutFileAOD);
 
  mgr -> ConnectInput (mixTask, 0, cInputChain);
  mgr -> ConnectOutput(mixTask, 1, cOutputUserAODTree);

  mgr -> SetDebugLevel(10);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("mix");
  
}

//===========================================================================================
